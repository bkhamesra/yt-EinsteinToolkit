"""
EinsteinToolkit-specific IO functions
"""

# -----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

import numpy as np

from yt.utilities.io_handler import \
    BaseIOHandler

from yt.funcs import mylog as logHandler

# Basic function to format filesize output.
# NOTE: filesize should be in bytes.
def format_filesize_output(filesize):
    for unit in ['B', 'kB', 'MB', 'GB', 'TB', 'PB']:
        if filesize < 1000.:
            return '%3.1f %s'%(filesize, unit)
        filesize /= 1000.

class IOHandlerEinsteinToolkit(BaseIOHandler):
    _particle_reader = False
    _dataset_type = "EinsteinToolkit"

    def _read_particle_coords(self, chunks, ptf):
        raise NotImplementedError(
            "EinsteinToolkit does not currently support particles.")

    def _read_particle_fields(self, chunks, ptf, selector):
        raise NotImplementedError(
            "EinsteinToolkit does not currently support particles.")

    def _read_fluid_selection(self, chunks, selector, fields, size):
        for ftype, fname in fields:
            if ftype != self._dataset_type:
                raise NotImplementedError('Data is not EinsteinToolkit data.')

        if size is None:
            size = self.getSize(chunks, selector)

        rv = {}
        for field in fields:
            rv[field] = np.empty(size, dtype=np.float64)

        logHandler.info("Reading %s cells for fields %s",
                        size, [f2 for f1, f2 in fields])

        offset = 0
        for chunk in chunks:
            data = self._read_chunk_data(chunk, fields)
            for g in chunk.objs:
                for field in fields:
                    ds = data[g.id].pop(field)
                    if selector.__class__.__name__ == "GridSelector":
                        length = ds.size
                        rv[field][offset:offset +
                                  length] = np.reshape(ds, length)
                    else:
                        length = g.select(selector, ds, rv[field], offset)
                offset += length
                data.pop(g.id)
        if offset != size:
            raise ValueError("Offset is wrong")
        return rv

    def _read_chunk_data(self, chunk, fields):
        data = {}
        if len(chunk.objs) == 0:
            return data

        for grid in chunk.objs:
            data[grid.id] = {}
            for field in fields:
                data[grid.id][field] = grid.component.read_data(field)

        return data

    def getSize(self, chunks, selector):
        return sum(grid.count(selector) for chunk in chunks for grid in chunk.objs)
