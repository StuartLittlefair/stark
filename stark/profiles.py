from __future__ import (print_function, absolute_import)
import pkg_resources
import os

import numpy as np

from .utils import read_n_items


class ProfileMetaData:

    def __init__(self, nlower, nupper, log_alpha_min, log_ne_min, log_t_min,
                 log_alpha_increment, log_ne_increment, log_t_increment,
                 num_alpha, num_ne, num_temp):
        self.nlower = int(nlower)
        self.nupper = int(nupper)
        self.log_alpha_min = float(log_alpha_min)
        self.log_ne_min = float(log_ne_min)
        self.log_t_min = float(log_t_min)
        self.log_alpha_increment = float(log_alpha_increment)
        self.log_ne_increment = float(log_ne_increment)
        self.log_t_increment = float(log_t_increment)
        self.num_alpha = int(num_alpha)
        self.num_ne = int(num_ne)
        self.num_temp = int(num_temp)

    @classmethod
    def from_file_handle(cls, fh):
        items = read_n_items(11, fh)
        return cls(*items)

    def __repr__(self):
        return "ProfileMetaData(nlower={}, nupper={})".format(self.nlower, self.nupper)


def get_profiles(nlower=2, nupper=3, with_doppler=False):
    files = ['lyman', 'balmer', 'pasch', 'brack']
    file = files[nlower-1]
    if not with_doppler:
        file += 'nd'
    data_dir = pkg_resources.resource_filename('stark', 'data')
    path = os.path.join(data_dir, file)

    with open(path) as fh:
        nlines = int(read_n_items(1, fh))
        metadata = []
        for i in range(nlines):
            metadata.append(ProfileMetaData.from_file_handle(fh))
        for i in range(nlines):
            meta = metadata[i]
            data_size = (meta.num_alpha+1) * meta.num_ne * meta.num_temp
            line_info = read_n_items(2, fh)
            # make sure we're still in the right place in the file
            check_line_info(line_info, meta)
            data = read_n_items(data_size, fh)

            # now if this is the right transition, reshape and return
            if meta.nlower == nlower and meta.nupper == nupper:
                break
        data = np.array(data).astype('float')
        data = data.reshape((meta.num_ne, meta.num_temp, meta.num_alpha+1))
        flags = data[:, :, 0]
        data = data[:, :, 1:]
        return meta, flags, data


def check_line_info(line_info, meta):
    lower_info, upper_info = line_info
    if lower_info != 'nl={:02d}'.format(meta.nlower):
        raise Exception('File reading has gone wrong, unexpected line info')
    if upper_info != 'nu={:02d}'.format(meta.nupper):
        raise Exception('File reading has gone wrong, unexpected line info')
