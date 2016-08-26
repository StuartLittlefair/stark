#!/usr/bin/env python
from __future__ import (print_function, absolute_import, division)

import numpy as np

from astropy.convolution import convolve
from astropy.convolution import Gaussian1DKernel
from astropy.stats import gaussian_fwhm_to_sigma

from .interpolation import interpolate_profile
from .utils import get_central_wavelength


def make_hydrogen_line_profile(nlower, nupper, temp, nelec, v_inst, v_macro,
                               thermal_broadening=False):
    """
    Make a synthetic, Stark broadened line profile for Hydrogen

    Parameters
    ----------
    nlower : int
        lower level of transition
    nupper : int
        upper level of transition
    temp : float
        temperature in kelvin
    nelec : float
        electron number density in cm**-3
    v_inst : float
        instrumental line broadening (km/s)
    v_macro : float
        local broadening of line profile due to macroturbulence (km/s)
    thermal_broadening : bool, default False
        include thermal broadening at specified temperature, as implemented
        by Lemke 1997. There seems to be something up with this. For T=2500K
        I would expect around 5 km/s of broadening, whereas it looks closer
        to 20 km/s. Use with caution.

    Returns
    -------
    velocity : `np.ndarray`
        velocity shift from line centre
    line_profile : `np.ndarray`
        line profile of hydrogen
    """

    central_wavelength = get_central_wavelength(nlower, nupper)
    wave = central_wavelength + np.linspace(-100, 100, 10000)

    log_alpha, stark_profile, f0 = interpolate_profile(nlower, nupper, nelec, temp,
                                                       with_doppler=thermal_broadening)
    actual_log_alpha_values = np.log10(np.fabs((wave-central_wavelength)/f0))

    # Stark broadened line profile
    line_profile = np.interp(actual_log_alpha_values, log_alpha, stark_profile)

    # now convolve with local and instrumental Gaussian profile
    v_total = np.sqrt(v_inst**2 + v_macro**2)
    velocity = (wave - central_wavelength) * 3.e5 / central_wavelength
    kernel_width_fwhm = v_total / np.diff(velocity).mean()
    kernel_width_sigma = kernel_width_fwhm * gaussian_fwhm_to_sigma
    kernel = Gaussian1DKernel(kernel_width_sigma)

    line_profile_convolved = convolve(line_profile, kernel, boundary='extend')

    # normalise
    line_profile_convolved /= np.fabs(line_profile_convolved.sum())
    return velocity, line_profile_convolved
