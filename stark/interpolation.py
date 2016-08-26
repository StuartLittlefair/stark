from __future__ import (print_function, absolute_import)

import numpy as np

from .profiles import get_profiles


def interpolate_profile(nlower, nupper, nelec, temp, with_doppler=False):
    """
    Interpolate profile tables of Lemke 1997 to get a Stark broadened line profile

    Parameters
    ----------
    nlower : int
        lower level of transition
    nupper : int
        upper level of transition
    nelec : float
        number density of electrons in cm**-3
    temp : float
        temperature in K

    Returns
    -------
    log_alpha : `np.ndarray`
        alpha values
    profile : `np.ndarray`
        stark profile
    fo : float
        conversion between delta-alpha and delta-lambda
    """
    meta, flags, data = get_profiles(nlower, nupper, with_doppler)
    f0 = 1.25e-9*nelec**(2./3.)
    log_ne = np.log10(nelec)
    log_t = np.log10(temp)
    log_ne_index = (log_ne - meta.log_ne_min) / meta.log_ne_increment
    log_t_index = (log_t - meta.log_t_min) / meta.log_t_increment

    low_ne_index = int(np.floor(log_ne_index))
    high_ne_index = int(np.ceil(log_ne_index))
    low_t_index = int(np.floor(log_t_index))
    high_t_index = int(np.ceil(log_t_index))

    # check we are within bounds
    if low_ne_index < 0 or high_ne_index > meta.num_ne:
        raise ValueError("electron density outside allowed range 10**10 to 10**18 cm**-3")
    if low_t_index < 0 or high_t_index > meta.num_temp:
        raise ValueError("temperature outside allowed range 2500 to 160000 K")

    # points bracketing requested values
    ne1 = meta.log_ne_min + low_ne_index*meta.log_ne_increment
    ne2 = meta.log_ne_min + high_ne_index*meta.log_ne_increment
    t1 = meta.log_t_min + low_t_index*meta.log_t_increment
    t2 = meta.log_t_min + high_t_index*meta.log_t_increment

    # profiles at these points
    p1 = data[low_ne_index, low_t_index]
    p2 = data[high_ne_index, low_t_index]
    p3 = data[low_ne_index, high_t_index]
    p4 = data[high_ne_index, high_t_index]

    if ne1 == ne2 and t1 == t2:
        # no interpolation needed
        profile = p1
    elif ne1 == ne2:
        # interpolate in temp
        profile = p1 + (p3 - p1) * (log_t - t1) / (t2 - t1)
    elif t1 == t2:
        # interpolate in nelec
        profile = p1 + (p2 - p1) * (log_ne - ne1) / (ne2 - ne1)
    else:
        # otherwise do the full bilinear interpolation
        # interpolate in temp at low_ne
        r1 = p1 + (p3 - p1) * (log_t - t1) / (t2 - t1)
        # interpolate in temp at high_ne
        r3 = p2 + (p4 - p2) * (log_t - t1) / (t2 - t1)
        # interpolate in ne
        profile = r1 + (r3-r1) * (log_ne - ne1) / (ne2 - ne1)

    # OK - now find matching alpha values
    alpha = meta.log_alpha_min + np.arange(meta.num_alpha) * meta.log_alpha_increment
    return alpha, profile, f0

