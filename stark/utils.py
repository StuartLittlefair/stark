from __future__ import (print_function, absolute_import)

import pkg_resources
import os

import numpy as np


def get_central_wavelength(nlower, nupper):
    """
    Get the central wavelength of a Hydrogen line from the lower and upper transitions

    Parameters
    ----------
    nlower : int
        lower energy level involved
    nupper : int
        upper energy level involved

    Returns
    -------
    wave : float
        central wavelength in angstroms
    """
    data_dir = pkg_resources.resource_filename('stark', 'data')
    path = os.path.join(data_dir, 'H_linelist.dat')
    wave, nlo, nup, *rest = np.loadtxt(path).T
    mask = (nlo == nlower) & (nup == nupper)
    if mask.sum() != 1:
        raise ValueError("Could not find a unique match to these energy levels")
    return wave[mask][0]


def read_by_whitespace(fileobj):
    """
    Create a generator that reads in file items split by any whitespace

    Parameters
    ----------
    fileobj : open file handler

    Returns
    -------
    generator : generator object
    """
    for line in fileobj:
        for token in line.split():
            yield token


def read_n_items(n, fileobj):
    tokenized = read_by_whitespace(fileobj)
    if n == 1:
        return next(tokenized)
    else:
        return [next(tokenized) for i in range(n)]

