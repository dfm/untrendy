#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import os
import pyfits
import numpy as np
from scipy.interpolate import interp1d

from .untrendy import untrend, fit_trend


def _load_kepler_lc(name):
    fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_data",
                      name)
    with pyfits.open(fn) as f:
        time, flux, ferr = (f[1].data["TIME"], f[1].data["SAP_FLUX"],
                            f[1].data["SAP_FLUX_ERR"])

    return time, flux / np.median(flux), ferr / np.median(flux)


def test_SW_condition():
    np.seterr(all="raise")
    time, flux, ferr = _load_kepler_lc("kplr010874614-2009131105131_llc.fits")
    flux, ferr = untrend(time, flux, ferr, fill_times=10 ** -1.25)


def test_crazy_artifact():
    time, flux, ferr = _load_kepler_lc("kplr009002278-2012179063303_llc.fits")
    flux, ferr = untrend(time, flux, ferr, fill_times=10 ** -1.25)


def test_fake_data():
    np.random.seed(458)

    # Generate some fake data.
    tmax = 64.
    t0 = np.arange(0., tmax + 5, 4.)
    f0 = 0.1 * np.random.randn(len(t0)) - 0.001 * t0 * t0
    p0 = interp1d(t0, f0, kind="cubic")
    t = np.arange(0, tmax, 30. / 60. / 24.)
    f = p0(t)

    # Add some breakpoints.
    breakpoints = np.linspace(0, tmax, 5)[1:-1]
    for b in breakpoints:
        f[t > b] += 2 * np.random.randn()

    # Normalize to some reasonable light curve like values.
    f += 100
    f /= np.median(f)

    # Copy the truth.
    truth = np.array(f)

    # Add some noise.
    sigma = 0.0001 + 0.001 * np.random.rand(len(t))
    f += np.random.randn(len(sigma)) * sigma

    # Do the fit.
    trend = fit_trend(t, f, sigma)

    assert np.mean(np.abs((truth - trend(t)) / truth)) < 1e-3
