"""
EinsteinToolkit data structures
"""

# -----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from __future__ import print_function

from yt.utilities.on_demand_imports import _h5py as h5py

import os
import re
import numpy as np
import weakref
import time

from yt.data_objects.grid_patch import \
    AMRGridPatch
from yt.geometry.grid_geometry_handler import \
    GridIndex
from yt.data_objects.static_output import \
    Dataset
from yt.utilities.lib.misc_utilities import get_box_grids_level
from yt.units.yt_array import YTArray
from yt.funcs import mylog as log_handler
from yt.funcs import setdefaultattr
from yt.visualization.plot_window import SlicePlot
from .fields import EinsteinToolkitFieldInfo
from .geometry import DatasetEntry, RefinementLevel, SymmetryTypes, variable_placeholder
from .io import format_filesize_output

global_parameters = 'Parameters and Global Attributes'

class EinsteinToolkitGrid(AMRGridPatch):
    _id_offset = 0

    def __init__(self, id, index, level, dimensions, component):
        AMRGridPatch.__init__(self, id, filename=index.index_filename,
                              index=index)
        self.Parent = []
        self.Children = []
        self.Level = level
        self.component = component
        self.start_index = None
        self.populate_default_field_parameters()
        self._set_default_field_parameters()

    def populate_default_field_parameters(self):
        self._default_field_parameters = {}
        
        self._default_field_parameters['center'] = self.ds.arr(np.zeros(3, dtype='float64'), 'cm')
        self._default_field_parameters['normal'] = self.ds.arr([0, 0, 1], 'code_length')

    def __repr__(self):
        return "EinsteinToolkitGrid_%04i_rf%02i (%s)" % (self.id, self.Level, self.ActiveDimensions)

    def _setup_dx(self):
        if len(self.Parent) > 0:
            if not hasattr(self.Parent[0], 'dds'):
                self.Parent[0]._setup_dx()
            self.dds = self.Parent[0].dds.d / self.ds.refine_by
        else:
            LE, RE = (
                self.index.grid_left_edge[self.id, :].d, self.index.grid_right_edge[self.id, :].d)
            self.dds = (RE - LE) / self.ActiveDimensions

        self.dds = self.dds.view(YTArray)
        self.dds.units = self.index.grid_left_edge.units

    def get_global_startindex(self):
        if self.start_index is not None:
            return self.start_index

        offset = self.LeftEdge - self.ds.domain_left_edge
        self.start_index = np.rint(
            offset / self.dds).astype(np.int64).ravel()

        return self.start_index


class EinsteinToolkitHierarchy(GridIndex):
    grid = EinsteinToolkitGrid

    def __init__(self, ds, dataset_type='EinsteinToolkit'):
        self.dataset_type = dataset_type
        self.dataset = weakref.proxy(ds)
        self.index_filename = ds.parameter_filename
        self.directory = os.path.dirname(self.index_filename)
        self.float_type = np.float64
        GridIndex.__init__(self, ds, dataset_type)

    def _detect_output_fields(self):
        self.field_list = [(self.dataset_type, basename) for basename in self.ds.dataset_basenames]

    def _count_grids(self):
        self.num_grids = self.ds.num_grids

    def _parse_index(self):
        self.grid_left_edge = self.ds.arr(
            self.ds.grid_left_edge, 'code_length')
        self.grid_right_edge = self.ds.arr(
            self.ds.grid_right_edge, 'code_length')
        self.grid_delta = self.ds.arr(self.ds.grid_delta, 'code_length')
        self.grid_dimensions = self.ds.grid_dimensions
        self.grid_particle_count = self.ds.grid_particle_count
        self.grid_levels = self.ds.grid_levels
        self.grids = np.empty(self.num_grids, dtype=object)
        self.max_level = self.grid_levels.max()

        for index in range(self.num_grids):
            self.grids[index] = self.grid(index, self, self.grid_levels[index, 0],
                                          self.grid_dimensions[index], self.dataset.grid_components[index])


    def reconstruct_children(self):
        mask = np.empty(len(self.grids), dtype=np.int32)
        for index, grid in enumerate(self.grids):
            get_box_grids_level(self.grid_left_edge[index, :], self.grid_right_edge[index, :],
                                self.grid_levels[index] + 1, self.grid_left_edge, self.grid_right_edge, self.grid_levels, mask)
            ids = np.where(mask.astype(bool))[0]
            for child_id in ids:
                overlap = (np.minimum(self.grid_right_edge[child_id].d, self.grid_right_edge[index].d) - np.maximum(self.grid_left_edge[child_id].d, self.grid_left_edge[index].d)).prod()/self.grid_delta[child_id].d.prod()
                # More than one cell of the potential child grid lies within the parent. To avoid rounding errors at boundaries.
                if overlap > 1.0:
                    grid.Children.append(self.grids[child_id])
                    self.grids[child_id].Parent.append(grid)

    def _populate_grid_objects(self):
        for grid in self.grids:
            grid._prepare_grid()
            grid._setup_dx()

        self.reconstruct_children()


class EinsteinToolkitDataset(Dataset):
    _index_class = EinsteinToolkitHierarchy
    _field_info_class = EinsteinToolkitFieldInfo
    _possible_symmetries = ['rotating180', 'rotating90', 'reflectx', 'reflecty', 'reflectz']

    def __init__(self, filename, iteration=-1, dataset_type='EinsteinToolkit', code_mass_solar=1.,
                 slice_plane=None, symmetry=None,
                 storage_filename=None,
                 units_override=None):
        self.fluid_types += ('EinsteinToolkit',)
        self.code_mass_solar = code_mass_solar
        self.slice_plane = slice_plane

        if iteration == -1:
            try:
                self.file_handle = h5py.File(filename, 'r')
            except:
                log_handler.error('Input file no longer valid.')
                self.file_handle = None
                return

            iterations = self.get_all_iterations()
            self.iteration = min(iterations)
            self.file_handle.close()
            self.file_handle = None
        else:
            self.iteration = iteration

        self.symmetries = []
        self.validate_symmetries(symmetry)

        # NOTE: This will be overwritten later when it is determined from the iteration.
        #      The next constructor call complains if this is not yet set.
        self.current_time = 0.0

        Dataset.__init__(self, filename, dataset_type,
                         units_override=units_override)
        self.storage_filename = storage_filename

    def validate_symmetries(self, symmetry):
        if isinstance(symmetry, list) or isinstance(symmetry, tuple):
            if 'rotating180' in symmetry and 'rotating90' in symmetry:
                log_handler.error('Invalid symmetry options. Cannot activate both rotating180 and rotating90.')
                raise ValueError('Invalid symmetry option.')
            if ('rotating180' in symmetry or 'rotating90' in symmetry) and ('reflectx' in symmetry or 'reflecty' in symmetry):
                log_handler.error('Invalid symmetry options. Cannot activate rotating180/90 symmetries at the same time as reflectx or reflecty.')
                raise ValueError('Invalid symmetry option.')

            for entry in symmetry:
                self.validate_symmetries(entry)
            return

        if symmetry is not None and symmetry not in self._possible_symmetries:
            log_handler.error('Invalid symmetry option %s. Can choose from: %s'%(symmetry, str(self._possible_symmetries)))
            raise ValueError('Invalid symmetry option.')

        if symmetry is 'rotating90':
            self.symmetries.append(SymmetryTypes.rotating90)
        elif symmetry is 'rotating180':
            self.symmetries.append(SymmetryTypes.rotating180)
        elif symmetry is 'reflectx':
            self.symmetries.append(SymmetryTypes.reflectx)
        elif symmetry is 'reflecty':
            self.symmetries.append(SymmetryTypes.reflecty)
        elif symmetry is 'reflectz':
            self.symmetries.append(SymmetryTypes.reflectz)

    def _set_code_unit_attributes(self):
        setdefaultattr(self, 'length_unit', self.quan(
            self.code_mass_solar, 'l_geom'))
        setdefaultattr(self, 'mass_unit', self.quan(
            self.code_mass_solar, 'm_geom'))
        setdefaultattr(self, 'time_unit', self.quan(
            self.code_mass_solar, 't_geom'))
        setdefaultattr(self, 'velocity_unit', self.quan(1.0, 'c'))

    def slice_plot(self, *args, **kwargs):
        if self.dimensionality != 2:
            log_handler.warning('dataset.slice_plot is only meant to be used with 2D data.')
            return None

        if self.slice_plane == 'xy':
            return SlicePlot(self, 'z', *args, **kwargs)
        elif self.slice_plane == 'xz':
            return SlicePlot(self, 'y', *args, **kwargs)
        elif self.slice_plane == 'yz':
            return SlicePlot(self, 'x', *args, **kwargs)

    def _parse_parameter_file(self):
        # Try to open file
        try:
            self.file_handle = h5py.File(self.parameter_filename, 'r')
        except:
            log_handler.error('Input file no longer valid.')
            self.file_handle = None
            return

        # Load string with all parameters
        self.parameter_string = self.file_handle[global_parameters]['All Parameters'][(
            )].decode('utf-8')

        self.unique_identifier = "{}-{}".format(
            self.parameter_filename, time.ctime())
        self.refine_by = self.get_refinement_factor()
        self.periodicity = 3 * (False,)

        self.determine_dimensionality()
        self.determine_slice_plane()
        self.parse_datasets()
        self.process_symmetries()
        self.initialize_grid_components()
        self.determine_domain_geometry()

        # Non-cosmological
        self.cosmological_simulation = 0
        self.current_redshift = 0
        self.omega_lambda = 0
        self.omega_matter = 0
        self.hubble_constant = 0

    def determine_dimensionality(self):
        # Pick a dataset at random (the first) from which we can extract the dimensionality.
        dataset_name = next((x for x in self.file_handle.keys() if x != global_parameters), None)
        if dataset_name is None:
            log_handler.error('No datasets found in file %s'%self.parameter_filename)
            return 

        dataset = self.file_handle[dataset_name]
        self.dimensionality = dataset.ndim

    def determine_slice_plane(self):
        if self.dimensionality == 2:
            if self.slice_plane is not None:
                return

            match = re.search('.*\.([xyz]{2})\.*.h5', self.parameter_filename)
            if match is None:
                log_handler.error(
                    'Cannot determine slicing plane from input file even though the data appears to be two dimensional. Specify slice_plane keyword argument in yt.load(), or use a file that matches the pattern ".*\.([xyz]{2})\.*".')
            else:
                self.slice_plane = match.groups()[0]

                if self.slice_plane not in ['xy', 'xz', 'yz']:
                    log_handler.error('Invalid slice plane: %s' %
                                     (self.slice_plane))
                    raise TypeError(
                        'File name incorrect for a two dimensional dataset.')

    #####################################################################
    # NOTE: The following function assumes that both the grid structure #
    #       and the parameters are identical for every variable.        #
    #                                                                   #
    #       This simplifies things greatly, and is _usually_ okay.      #
    #####################################################################
    def parse_datasets(self):
        dataset_regex = re.compile('(.*)\sit=\d.*')
        basename_set = set()

        for dataset_name in self.file_handle:
            if dataset_name == global_parameters:
                continue

            match = dataset_regex.search(dataset_name)
            if match is None:
                log_handler.warning('Dataset name %s not properly formatted. Skipping...'%dataset_name)
                continue

            basename_set.add(match.groups()[0])

        self.dataset_basenames = list(basename_set)
        if len(self.dataset_basenames) is 0:
            log_handler.error('No usable datasets found in file %s'%self.parameter_filename)
            return

        representative_dataset_basename = self.dataset_basenames[0]

        entry_set = set()
        have_current_time = False
        regex = re.compile('%s\sit=(\d+).*rl=(\d+).*'%(re.escape(representative_dataset_basename)))
        levels = set()

        for key in self.file_handle:
            generic_key = key.replace(representative_dataset_basename, variable_placeholder)
            match = regex.search(key)
            if match is None:
                continue

            iteration = int(match.groups()[0])
            level = int(match.groups()[1])
            if iteration != self.iteration:
                continue

            if not have_current_time:
                self.current_time = self.file_handle[key].attrs['time']
                have_current_time = True

            entry_set.add(DatasetEntry(self, key, generic_key, iteration, level))
            levels.add(level)

        if not have_current_time:
            log_handler.error('No datasets found for the specified iteration.')
            raise ValueError('No datasets found for the specified iteration.')

        entry_set = self.remove_duplicate_grids(entry_set)
        final_entry_set = set()
        for level in levels:
            rf = RefinementLevel(level, entry_set)
            if rf.require_removal():
                final_entry_set.update(rf.determine_base_set())
            else:
                final_entry_set.update([chunk.entry for chunk in rf.chunks])

        self.dataset_entries = list(final_entry_set)
        self.num_grids = len(self.dataset_entries)

    def process_symmetries(self):
        warning_string = 'File contains grid patches that are incompatible with the given symmetries. These will be ignored in the context of symmetry.'

        # Be sure to process rotating symmetries before reflections so that the full domain is filled in.
        # Can also take advantage of the fact that we can only have one rotation symmetry or the other.
        if SymmetryTypes.rotating90 in self.symmetries or SymmetryTypes.rotating180 in self.symmetries:
            symmetry = SymmetryTypes.rotating90 if SymmetryTypes.rotating90 in self.symmetries else SymmetryTypes.rotating180

            for index in range(self.num_grids):
                entry = self.dataset_entries[index]
                epsneg = -entry.delta/self.refine_by
                if (   (symmetry is SymmetryTypes.rotating90 and (entry.left_edge[0] < epsneg[0] or entry.left_edge[1] < epsneg[1]))
                    or (symmetry is SymmetryTypes.rotating180 and entry.left_edge[0] < epsneg[0])):
                    log_handler.warning(warning_string)
                else:
                    self.dataset_entries.extend(entry.create_symmetry_copy(symmetry))
            self.num_grids = len(self.dataset_entries)

        # Remaining reflection symmetries
        reflections = [sym for sym in self.symmetries if sym not in (SymmetryTypes.rotating90, SymmetryTypes.rotating180)]
        for symmetry in reflections:
            for index in range(self.num_grids):
                entry = self.dataset_entries[index]
                epsneg = -entry.delta/self.refine_by
                if (   (symmetry is SymmetryTypes.reflectx and entry.left_edge[0] < epsneg[0])
                    or (symmetry is SymmetryTypes.reflecty and entry.left_edge[1] < epsneg[1])
                    or (symmetry is SymmetryTypes.reflectz and entry.left_edge[2] < epsneg[2])):
                    log_handler.warning(warning_string)
                else:
                    self.dataset_entries.extend(entry.create_symmetry_copy(symmetry))
            self.num_grids = len(self.dataset_entries)

    def remove_duplicate_grids(self, entry_set):
        target_set = set()

        while len(entry_set) > 0:
            entry = entry_set.pop()
            potential_subgrid_list = [other_entry for other_entry in entry_set.union(target_set) if entry.can_be_subgrid(other_entry)]

            for grid in potential_subgrid_list:
                entry_set.discard(grid)
                target_set.discard(grid)

            if len(potential_subgrid_list) == 0:
                target_set.add(entry)
            else:
                potential_subgrid_list.append(entry)
                target_set.add(max(potential_subgrid_list, key=lambda x: x.npoints.prod()))

        return target_set

    def determine_domain_geometry(self):
        self.dx_max, self.dy_max, self.dz_max = min(self.grid_components, key=lambda x: x.refinement_level).delta
        x_min = min(self.grid_components, key=lambda x: x.x_min).x_min
        x_max = max(self.grid_components, key=lambda x: x.x_max).x_max
        y_min = min(self.grid_components, key=lambda x: x.y_min).y_min
        y_max = max(self.grid_components, key=lambda x: x.y_max).y_max
        z_min = min(self.grid_components, key=lambda x: x.z_min).z_min
        z_max = max(self.grid_components, key=lambda x: x.z_max).z_max

        self.domain_left_edge = np.array([x_min, y_min, z_min])
        self.domain_right_edge = np.array([x_max, y_max, z_max])
        self.domain_delta = np.array([self.dx_max, self.dy_max, self.dz_max])
        self.domain_width = self.domain_right_edge - self.domain_left_edge
        self.domain_dimensions = np.rint(self.domain_width/self.domain_delta).astype(int)

    def initialize_grid_arrays(self):
        self.grid_left_edge = np.empty(
            shape=(self.num_grids, 3), dtype=np.float64)
        self.grid_right_edge = np.empty(
            shape=(self.num_grids, 3), dtype=np.float64)
        self.grid_dimensions = np.ones(
            shape=(self.num_grids, 3), dtype=np.int32)
        self.grid_delta = np.empty(
            shape=(self.num_grids, 3), dtype=np.float64)
        self.grid_particle_count = np.zeros(
            shape=(self.num_grids, 1), dtype=np.int32)
        self.grid_levels = np.empty(
            shape=(self.num_grids, 1), dtype=np.int32)
        self.grid_components = np.empty(self.num_grids, dtype=object)

    def initialize_grid_components(self):
        self.initialize_grid_arrays()

        for index, component in enumerate(self.dataset_entries):
            self.grid_left_edge[index] = component.left_edge
            self.grid_right_edge[index] = component.right_edge
            self.grid_dimensions[index] = component.ncells_ghost
            self.grid_levels[index] = component.refinement_level
            self.grid_delta[index] = component.delta
            self.grid_components[index] = component
        
    def get_refinement_factor(self):
        match = re.search(
            'Carpet::refinement_factor[\s*]=[\s*](\d+)', self.parameter_string, re.IGNORECASE)

        if match is None:
            log_handler.warning('Cannot determine refinement factor from parameters. Assuming value of 2.')
            return 2

        return int(match.groups()[0])

    def get_all_iterations(self):
        iterations = set()

        for dataset_name in self.file_handle.keys():
            match = re.search('.*it=(\d+)', dataset_name)
            if match is None:
                continue

            iterations.add(int(match.groups()[0]))

        iterations = list(iterations)
        iterations.sort()
        return iterations

    def get_file_metadata(self):
        log_handler.info('Compiling metadata for input file %s. This may take a while...' % (
            self.parameter_filename))
        total_datasets = len(self.file_handle) - 1
        all_datasets = [key for key in self.file_handle if key != global_parameters]

        first_time = self.file_handle[min(
            all_datasets, key=lambda x: self.file_handle[x].attrs['time'])].attrs['time']
        final_time = self.file_handle[max(
            all_datasets, key=lambda x: self.file_handle[x].attrs['time'])].attrs['time']

        all_iterations = self.get_all_iterations()

        log_handler.info('\tFile path:                   %s'%self.parameter_filename)
        log_handler.info('\tSize on disk:                %s'%(format_filesize_output(os.path.getsize(self.parameter_filename))))
        log_handler.info('\tTotal number of datasets:    %s' % (total_datasets))
        log_handler.info('\tDataset basenames found:     %s' %
                        (self.dataset_basenames))
        log_handler.info('\tData dimensionality:         %d' % (self.dimensionality))
        if self.dimensionality == 2:
            log_handler.info('\tSlicing plane:               %s' % (self.slice_plane))
        log_handler.info('\tTime range:                  [%f, %f]' % (first_time, final_time))
        log_handler.info('\tIteration range:             [%d, %d]' % (
            all_iterations[0], all_iterations[-1]))
        if len(all_iterations) > 1:
            log_handler.info('\tOutput frequency:            %d' %
                            (all_iterations[1] - all_iterations[0]))
        log_handler.info('\tNumber of output iterations: %d' %
                        (len(all_iterations)))
        log_handler.info('\tRefinement levels:           %d' % (self.index.max_level))

    def dimensionality_fill(self, source, fill_value):
        if len(source) == 3:
            return source

        if not isinstance(fill_value, (list, tuple, np.ndarray)):
            fill_value = [fill_value] * 3

        if self.slice_plane == 'xy':
            return np.array([source[0], source[1], fill_value[2]], dtype=source.dtype)
        elif self.slice_plane == 'xz':
            return np.array([source[0], fill_value[1], source[1]], dtype=source.dtype)
        elif self.slice_plane == 'yz':
            return np.array([fill_value[0], source[0], source[1]], dtype=source.dtype)

    def dimensionality_shape(self, source):
        if len(source) == 3:
            return np.array([source[3 - 1 - index] for index in range(3)])

        if self.slice_plane == 'xy':
            return np.array([source[1], source[0], 1])
        elif self.slice_plane == 'xz':
            return np.array([source[1], 1, source[0]])
        elif self.slice_plane == 'yz':
            return np.array([1, source[1], source[0]])

    def convert_read_shape(self, data):
        if self.dimensionality == 3:
            return np.transpose(data, (2, 1, 0))
        elif self.slice_plane == 'xy':
            target = np.empty(
                (data.shape[1], data.shape[0], 1), dtype=data.dtype)
            target[:, :, 0] = np.transpose(data)
            return target
        elif self.slice_plane == 'xz':
            target = np.empty(
                (data.shape[1], 1, data.shape[0]), dtype=data.dtype)
            target[:, 0, :] = np.transpose(data)
            return target
        elif self.slice_plane == 'yz':
            target = np.empty(
                (1, data.shape[1], data.shape[0]), dtype=data.dtype)
            target[0, :, :] = np.transpose(data)
            return target

    @classmethod
    def _is_valid(self, *args, **kwargs):
        path = args[0]
        try:
            file_handle = h5py.File(path, 'r')
            if global_parameters in file_handle.keys():
                file_handle.close()
                return True

            file_handle.close()
        except:
            pass

        return False

