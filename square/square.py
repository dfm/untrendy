#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Various tools and constants for working with Kepler data.

"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["detrend"]

import numpy as np

try:
    from scipy.interpolate import LSQUnivariateSpline
    LSQUnivariateSpline = LSQUnivariateSpline
except ImportError:
    LSQUnivariateSpline = None

from . import _square


def detrend(x, y, yerr=None, Q=4, dt=4., tol=1.25e-3, maxiter=15,
            fill_times=True, maxditer=4, nfill=4):
    """
    Use iteratively re-weighted least squares to fit a spline to the base
    trend in a time series. This is especially useful (and specifically
    tuned) for de-trending Kepler light curves.

    :param x:
        The sampled times.

    :param y:
        The fluxes corresponding to the times in ``x``.

    :param yerr: (optional)
        The 1-sigma error bars on ``y``.

    :param Q: (optional)
        The parameter controlling the severity of the re-weighting.

    :param dt: (optional)
        The initial spacing between time control points.

    :param tol: (optional)
        The convergence criterion.

    :param maxiter: (optional)
        The maximum number of re-weighting iterations to run.

    :param maxditer: (optional)
        The maximum number of discontinuity search iterations to run.

    :param nfill: (optional)
        The number of knots to use to fill in the gaps.

    """
    if LSQUnivariateSpline is None:
        raise ImportError("scipy is required for spline de-trending.")

    if yerr is None:
        yerr = np.ones_like(y)

    inds = np.argsort(x)
    x, y, yerr = x[inds], y[inds], yerr[inds]
    ivar = 1. / yerr / yerr
    w = np.array(ivar)

    # Build the list of knot locations.
    N = (x[-1] - x[0]) / dt + 2
    t = np.linspace(x[0], x[-1], N)[1:-1]

    # Refine knot locations around break points.
    if fill_times:
        inds = x[1:] - x[:-1] > 10 ** (-1.25)
        for i in np.arange(len(x))[inds]:
            t = _add_knots(t, x[i], x[i + 1], N=nfill)

    for j in range(maxditer):
        s0 = None
        for i in range(maxiter):
            # Fit the spline.
            extra_t = np.append(t, [x[0], x[-1]])
            x0 = np.append(x, extra_t)
            inds = np.argsort(x0)
            y0 = np.append(y, np.ones_like(extra_t))[inds]
            w0 = np.append(w, np.ones_like(extra_t))[inds]
            p = LSQUnivariateSpline(x0[inds], y0, t, k=3, w=w0)

            # Compute chi_i ^2.
            chi = (y - p(x)) / yerr
            chi2 = chi * chi

            # Check for convergence.
            sigma = np.median(chi2)
            if s0 is not None and np.abs(s0 - sigma) < tol:
                break
            s0 = sigma

            # Re compute weights.
            w = ivar * Q / (chi2 + Q)

        # Find any discontinuities.
        i = _square.discontinuities(x, chi, 0.5 * dt, Q, 1.0)
        if i < 0:
            return p

        t = _add_knots(t, x[i], x[i + 1], N=np.max([nfill, 4]))

    return p


def _add_knots(t, t1, t2, N=3):
    """
    A hack for adding ``N`` samples in a given region and removing any other
    samples in that region.

    """
    return np.sort(np.append(t[(t < t1) + (t > t2)], np.linspace(t1, t2, N)))
