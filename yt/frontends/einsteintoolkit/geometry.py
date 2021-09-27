"""
EinsteinToolkit geometry code
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

import numpy as np
import weakref
from enum import Enum
from math import sin, cos, pi

variable_placeholder = 'VARNAMEPLACEHOLDER'

class SymmetryTypes(Enum):
    rotating90 = 1
    rotating180 = 2
    reflectx = 3
    reflecty = 4
    reflectz = 5

# A couple of helper functions
def transform_vector(transformation_matrix, vector):
    return np.append(np.dot(transformation_matrix, vector[:2]), vector[2])

def transform_vector_absolute(transformation_matrix, vector):
    return np.append(np.absolute(np.dot(transformation_matrix, vector[:2])), vector[2])

class DatasetEntry:
    def __init__(self, ds, dataset_name, generic_dataset_name, iteration, level):
        if type(ds) is weakref.ProxyType:
            self.ds = ds
        else:
            self.ds = weakref.proxy(ds)

        self.file_handle = ds.file_handle
        self.dataset_name = dataset_name
        self.generic_dataset_name = generic_dataset_name
        self.refinement_level = level
        self.iteration = iteration
        dataset = self.file_handle[dataset_name]
        attributes = dataset.attrs

        self.symmetry_copy = False
        self.symmetry_type = None
        self.symmetry_original = None
        self.symmetry_quadrant = None

        self.ndims = len(dataset.shape)
        self.npoints = self.ds.dimensionality_shape(dataset.shape)
        self.delta = self.ds.dimensionality_fill(attributes['delta'], 1.0)
        self.iorigin = self.ds.dimensionality_fill(attributes['iorigin'], 0)
        self.origin = self.ds.dimensionality_fill(attributes['origin'], 0.)
        self.ghost_zones = self.ds.dimensionality_fill(
            attributes['cctk_nghostzones'], 0)

        self.ghosts = np.zeros(6, dtype=np.int)
        for index in range(3):
            self.ghosts[2 * index] = self.ghost_zones[index]
            self.ghosts[2 * index + 1] = self.ghost_zones[index]

        self.npoints_ghost = np.empty(3, dtype=np.int)
        self.ncells_ghost = np.empty(3, dtype=np.int)
        self.iorigin_ghost = np.empty(3, dtype=np.int)
        self.origin_ghost = np.empty(3, dtype=np.float)

        for index in range(3):
            if self.ghosts[2 * index] > 0 or self.ghosts[2 * index + 1] > 0:
                self.npoints_ghost[index] = self.npoints[index] - \
                    self.ghosts[2 * index] - self.ghosts[2 * index + 1] + 1
            else:
                self.npoints_ghost[index] = self.npoints[index]

            self.iorigin_ghost[index] = self.iorigin[index] + \
                self.ghosts[2 * index]
            self.origin_ghost[index] = self.origin[index] + \
                self.ghosts[2 * index] * self.delta[index]

            if self.ghosts[2 * index] > 0 or self.ghosts[2 * index + 1] > 0:
                if self.refinement_level > 0 and self.iorigin_ghost[index] % 2 > 0:
                    self.npoints_ghost[index] += 1
                    self.iorigin_ghost[index] -= 1
                    self.origin_ghost[index] -= self.delta[index]
                    self.ghosts[2 * index] -= 1
                if self.refinement_level > 0 and self.npoints_ghost[index] % 2 == 0:
                    self.npoints_ghost[index] -= 1

        #Subtract one from npoints to convert to cell-centered (except when only one cell)
        self.ncells_ghost = np.maximum(1, self.npoints_ghost - 1)
        self.disk_slice = (slice(self.ghosts[0], self.ghosts[0] + self.ncells_ghost[0]), slice(
            self.ghosts[2], self.ghosts[2] + self.ncells_ghost[1]), slice(self.ghosts[4], self.ghosts[4] + self.ncells_ghost[2]))

        self.x_min = 0. if self.npoints_ghost[0] == 1 else self.origin_ghost[0]
        self.y_min = 0. if self.npoints_ghost[1] == 1 else self.origin_ghost[1]
        self.z_min = 0. if self.npoints_ghost[2] == 1 else self.origin_ghost[2]
        self.left_edge = np.array([self.x_min, self.y_min, self.z_min])

        self.x_max = 1. if self.npoints_ghost[0] == 1 else self.origin_ghost[0] + (
            self.npoints_ghost[0] - 1) * self.delta[0]
        self.y_max = 1. if self.npoints_ghost[1] == 1 else self.origin_ghost[1] + (
            self.npoints_ghost[1] - 1) * self.delta[1]
        self.z_max = 1. if self.npoints_ghost[2] == 1 else self.origin_ghost[2] + (
            self.npoints_ghost[2] - 1) * self.delta[2]
        self.right_edge = np.array([self.x_max, self.y_max, self.z_max])

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.refinement_level ==  other.refinement_level and np.array_equal(self.npoints, other.npoints) and np.array_equal(self.iorigin, other.iorigin)
        else:
            raise TypeError('Cannot compare types %s and %s.'%(str(type(self)), str(type(other))))

    def __hash__(self):
        return hash((self.refinement_level, self.npoints[0], self.npoints[1], self.npoints[2], self.iorigin[0], self.iorigin[1], self.iorigin[2]))

    def can_be_subgrid(self, other):
        if isinstance(other, self.__class__):
            return self.refinement_level == other.refinement_level and (np.array_equal(self.iorigin, other.iorigin) or np.array_equal(self.iorigin + self.npoints, other.iorigin + other.npoints))
        else:
            raise TypeError('Cannot compare types %s and %s.'%(str(type(self)), str(type(other))))

    def read_data(self, field):
        field_name = field[1]
        dataset_name = self.generic_dataset_name.replace(variable_placeholder, field_name)

        if not self.symmetry_copy:
            return self.ds.convert_read_shape(self.file_handle[dataset_name])[self.disk_slice]
        elif self.symmetry_type is SymmetryTypes.rotating90:
            return np.rot90(self.symmetry_original.read_data(field), self.symmetry_quadrant)
        elif self.symmetry_type is SymmetryTypes.rotating180:
            return np.rot90(self.symmetry_original.read_data(field), 2)
        elif self.symmetry_type is SymmetryTypes.reflectx:
            return np.flip(self.symmetry_original.read_data(field), 0)
        elif self.symmetry_type is SymmetryTypes.reflecty:
            return np.flip(self.symmetry_original.read_data(field), 1)
        elif self.symmetry_type is SymmetryTypes.reflectz:
            return np.flip(self.symmetry_original.read_data(field), 2)
        else:
            raise TypeError('Invalid symmetry type encountered in read_data().')
            return None

    def create_symmetry_copy(self, symmetry_type):
        # Rotating90 symmetry is unique in that it returns three copies.
        if symmetry_type is SymmetryTypes.rotating90:
            copies = []
            for quadrant in range(1, 4):
                # First calculate the rotation matrix targeting this specific quadrant.
                # With theta in (pi/2, pi, 3pi/2), it is given by
                # cos(theta)     -sin(theta)
                # sin(theta)      cos(theta) 
                theta = quadrant*pi/2.
                transformation_matrix = np.array([[cos(theta), -sin(theta)], 
                                                  [sin(theta),  cos(theta)]]).astype(int)
                
                copy = DatasetEntry(self.ds, self.dataset_name, self.generic_dataset_name, self.iteration, self.refinement_level)
                copy.symmetry_copy = True
                copy.symmetry_original = weakref.proxy(self)
                copy.symmetry_type = symmetry_type
                copy.symmetry_quadrant = quadrant

                # Transform the x-y components of the relevant quantities using the two helper functions
                # transform_vector() and transform_vector_absolute()
                copy.npoints      = transform_vector_absolute(transformation_matrix, self.npoints)
                copy.delta        = transform_vector_absolute(transformation_matrix, self.delta)
                copy.ncells_ghost = transform_vector_absolute(transformation_matrix, self.ncells_ghost)
                copy.iorigin      = transform_vector(transformation_matrix, self.iorigin)
                copy.origin       = transform_vector(transformation_matrix, self.origin)

                # Cycle through the corners of the original entry to determine what the corners of the copy
                # should be.
                corners = [self.left_edge[:2], 
                           np.array([self.left_edge[0], self.right_edge[1]]), 
                           self.right_edge[:2], 
                           np.array([self.right_edge[0], self.left_edge[1]])]

                # Finally determine the left and right edges by transforming the corners.
                copy.left_edge =  np.append(np.dot(transformation_matrix, corners[quadrant]), np.array([copy.z_min]))
                copy.right_edge = np.append(np.dot(transformation_matrix, corners[(2+quadrant)%len(corners)]), np.array([copy.z_max]))

                copy.x_min = copy.left_edge[0]
                copy.x_max = copy.right_edge[0]

                copy.y_min = copy.left_edge[1]
                copy.y_max = copy.right_edge[1]

                copies.append(copy)
            return copies

        # Remaining symmetries only create one copy.
        copy = DatasetEntry(self.ds, self.dataset_name, self.generic_dataset_name, self.iteration, self.refinement_level)
        copy.symmetry_copy = True
        copy.symmetry_original = weakref.proxy(self)
        copy.symmetry_type = symmetry_type

        if symmetry_type is SymmetryTypes.rotating180:
            copy.x_min = -self.x_max
            copy.x_max = -self.x_min
            copy.y_min = -self.y_max
            copy.y_max = -self.y_min

            copy.left_edge[0] = copy.x_min
            copy.left_edge[1] = copy.y_min
            copy.right_edge[0] = copy.x_max
            copy.right_edge[1] = copy.y_max

            return [copy]

        # If this is a copy of a copy, we need to make sure that it has all
        # of the relevant properties of the first copy _not_ the original dataset.
        if self.symmetry_copy:
            np.copyto(copy.left_edge, self.left_edge)
            np.copyto(copy.right_edge, self.right_edge)
            copy.x_min = self.x_min
            copy.y_min = self.y_min
            copy.z_min = self.z_min
            copy.x_max = self.x_max
            copy.y_max = self.y_max
            copy.z_max = self.z_max

        if symmetry_type is SymmetryTypes.reflectx:
            copy.x_min = -self.x_max
            copy.x_max = -self.x_min
            copy.left_edge[0] = copy.x_min
            copy.right_edge[0] = copy.x_max
        elif symmetry_type is SymmetryTypes.reflecty:
            copy.y_min = -self.y_max
            copy.y_max = -self.y_min
            copy.left_edge[1] = copy.y_min
            copy.right_edge[1] = copy.y_max
        elif symmetry_type is SymmetryTypes.reflectz:
            copy.z_min = -self.z_max
            copy.z_max = -self.z_min
            copy.left_edge[2] = copy.z_min
            copy.right_edge[2] = copy.z_max
        else:
            raise TypeError('Unrecognized symmetry type %s passed to create_symmetry_copy().'%(str(symmetry_type)))
            return None

        return [copy]
        

class RefinementLevel:
    class Chunk:
        def __init__(self, index, entry):
            self.index = index
            self.entry = entry
            self.center = (entry.right_edge + entry.left_edge)/2.
            self.volume = entry.delta.prod()
            self.overlap = 0.
            self.distance = None

        def determine_overlap(self, others):
            self.overlap = 0.
            for other in others:
                if self == other:
                    continue
                common = np.minimum(self.entry.right_edge, other.entry.right_edge) - np.maximum(self.entry.left_edge, other.entry.left_edge)
                if np.any(common < 0):
                    continue
                overlap_volume = common.prod()
                self.overlap += round(overlap_volume/self.volume)
        
        def __eq__(self, other):
            return self.index == other.index

    def __init__(self, level, entry_set):
        self.level = level
        entries = [entry for entry in entry_set if entry.refinement_level == level]
        self.chunks = [self.Chunk(index, entry) for index,entry in enumerate(entries)]
        self.base_set = []

    def require_removal(self):
        self.update_overlaps(self.chunks)
        for chunk in self.chunks:
            if chunk.overlap > 0.:
                return True

        return False

    def determine_base_set(self):
        self.update_overlaps(self.chunks)
        self.chunks.sort(key=lambda x:x.overlap)

        seed = self.chunks[0]
        self.base_set.append(seed)
        self.chunks.remove(seed)

        self.update_distances(seed)
        self.update_overlaps(self.base_set)

        while True:
            addition = next((chunk for chunk in self.chunks if chunk.overlap == 0.), None)
            if addition is None:
                break
            self.base_set.append(addition)
            self.chunks.remove(addition)
            self.update_overlaps(self.base_set)
        
        return [chunk.entry for chunk in self.base_set]

    def update_distances(self, seed):
        for chunk in self.chunks:
            chunk.distance = np.linalg.norm(seed.center - chunk.center)

        self.chunks.sort(key=lambda x:x.distance)

    def update_overlaps(self, others):
        for chunk in self.chunks:
            chunk.determine_overlap(others)
