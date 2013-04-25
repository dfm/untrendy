#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tools for fitting and removing trends from light curves using the ￭ algorithm.
The default settings are tuned to work well for Kepler light curves but the
same algorithm might be useful for other datasets.

"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = ["untrend", "fit_trend", "discontinuity_scalar", "median"]

import logging
import numpy as np

try:
    from scipy.interpolate import LSQUnivariateSpline
    LSQUnivariateSpline = LSQUnivariateSpline
except ImportError:
    LSQUnivariateSpline = None

from . import _untrendy


def untrend(x, y, yerr=None, **kwargs):
    """
    Use iteratively re-weighted least squares to remove the out-of-transit
    trends in a light curve. Unlike ``fit_trend``, this function masks bad
    data (``NaN``) and normalizes the data before fitting.

    **Parameters**

    :``x``:          The sampled times.
    :``y``:          The fluxes corresponding to the times in ``x``.
    :``yerr``:       (optional) The 1-sigma error bars on ``y``.
    :``**kwargs``:   (optional) Other arguments passed to the ``fit_trend``
                     function.

    **Returns**

    :``flux``:       The de-trended relative fluxes.
    :``ferr``:       The de-trended uncertainties on ``flux``.

    """
    if yerr is None:
        yerr = np.ones_like(y)

    x, y, yerr = np.atleast_1d(x), np.atleast_1d(y), np.atleast_1d(yerr)

    # Mask bad data.
    inds = ~(np.isnan(x) + np.isnan(y) + np.isnan(yerr))
    x0, y0 = x[inds], y[inds]

    # Normalize the data.
    mu = np.median(y0)

    # Fit the trend.
    trend = fit_trend(x0, y0 / mu, yerr=yerr[inds] / mu, **kwargs)

    # De-trend the fluxes.
    factor = mu * trend(x0)
    y[inds] /= factor
    yerr[inds] /= factor
    return y, yerr


def fit_trend(x, y, yerr=None, Q=12, dt=3., tol=1.25e-3, maxiter=15,
              fill_times=None, maxditer=4, nfill=4):
    """
    Use iteratively re-weighted least squares to fit a spline to the
    out-of-transit trends in a time series. The input data should be "clean".
    In other words, bad data should be masked and it often helps to normalize
    the fluxes (by the median or something).

    **Parameters**

    :``x``:          The sampled times.
    :``y``:          The fluxes corresponding to the times in ``x``.
    :``yerr``:       (optional) The 1-sigma error bars on ``y``.
    :``Q``:          (optional) The parameter controlling the severity of the
                     re-weighting.
    :``dt``:         (optional) The initial spacing between time control
                     points.
    :``tol``:        (optional) The convergence criterion.
    :``maxiter``:    (optional) The maximum number of re-weighting iterations
                     to run.
    :``fill_times``: (optional) If provided, this number sets the minimum time
                     spacing between adjacent samples that is acceptable. If
                     the spacing is larger, knots will be added to fill in
                     the gap.
    :``maxditer``:   (optional) The maximum number of discontinuity search
                     iterations to run.
    :``nfill``:      (optional) The number of knots to use to fill in the
                     gaps.

    **Returns**

    :``trend``:      A callable representation of the trend.

    """
    if LSQUnivariateSpline is None:
        raise ImportError("scipy is required for spline de-trending.")

    assert maxditer > 0

    if yerr is None:
        yerr = np.ones_like(y)

    x, y, yerr = np.atleast_1d(x), np.atleast_1d(y), np.atleast_1d(yerr)

    # Mask bad data.
    mask = ~(np.isnan(x) + np.isnan(y) + np.isnan(yerr))
    x_masked, y_masked, yerr_masked = x[mask], y[mask], yerr[mask]

    # Remove any points that aren't sequential.
    delta = np.append(x_masked[1:] - x_masked[:-1], True)
    x_masked = x_masked[delta > 0]
    y_masked = y_masked[delta > 0]
    yerr_masked = yerr_masked[delta > 0]

    # Get the median flux.
    mu = np.median(y_masked)

    # Compute the initial weights.
    ivar = 1. / yerr_masked / yerr_masked
    w = np.array(ivar)

    # Build the list of knot locations.
    N = (x[-1] - x[0]) / dt + 2
    t = np.linspace(x[0], x[-1], N)[1:-1]

    # Refine knot locations around break points.
    if fill_times is not None:
        inds = x[1:] - x[:-1] > fill_times
        logging.info("Filling in {0} time gaps.".format(np.sum(inds)))
        for i in np.arange(len(x))[inds]:
            t = _add_knots(t, x[i], x[i + 1], N=2)

    discontinuity = -1
    for j in range(maxditer):
        s0 = None
        for i in range(maxiter):
            # Remove any knots that are at exactly the same point.
            delta = t[1:] - t[:-1]
            if np.any(delta <= 0):
                t = t[delta > 0]

            # Add "data" at the positions of all the knots at the median of
            # the fluxes. This should help keep the values reasonable even
            # in time gaps.
            x0 = np.append(x_masked, t)
            inds = np.argsort(x0)
            x0 = x0[inds]
            y0 = np.append(y_masked, mu * np.ones_like(t))[inds]
            w0 = np.append(w, np.ones_like(t))[inds]

            # Fit the spline.
            p = LSQUnivariateSpline(x0, y0, t, k=3, w=w0)

            # Compute chi_i ^2.
            chi = (y_masked - p(x_masked)) / yerr_masked
            chi2 = chi * chi

            # Check for convergence.
            sigma = np.median(chi2)
            if s0 is not None and np.abs(s0 - sigma) < tol:
                logging.info("Converged after {0} re-weighting iterations."
                             .format(i))
                break
            s0 = sigma

            # Re compute weights.
            w = ivar * Q / (chi2 + Q)

        # Find any discontinuities.
        i = _untrendy.find_discontinuities(x_masked[1:-2], chi[1:-2],
                                           0.5 * dt, Q, 2)
        if i < 0 or discontinuity == x_masked[i + 1]:
            return p

        discontinuity = x_masked[i + 1]

        logging.info("Discontinuity found at t={0}".format(x_masked[i + 1]))
        t = _add_knots(t, x_masked[i + 1], x_masked[i + 2],
                       N=np.max([nfill, 4]))

    return p


def _add_knots(t, t1, t2, N=3):
    """
    A hack for adding ``N`` samples in a given region and removing any other
    samples in that region.

    """
    return np.sort(np.append(t[(t < t1) + (t > t2)], np.linspace(t1, t2, N)))


def discontinuity_scalar(x, y, yerr=None, **kwargs):
    kwargs["maxditer"] = kwargs.get("maxditer", 15) - 1
    trend = fit_trend(x, y, yerr=yerr, **kwargs)

    Q = kwargs.get("Q", 12)
    dt = kwargs.get("dt", 4.)

    chi = (y - trend(x)) / yerr
    softr = np.sqrt(Q / (Q + chi * chi)) * chi

    tmid = 0.5 * (x[1:] + x[:-1])
    k = (x[:, None] - tmid[None, :]) / dt
    k = (k - 1) ** 2 * (0 <= k) * (k <= 1) - (k + 1) ** 2 * (-1 <= k) * (k < 0)
    val = np.sum(k * softr[:, None], axis=0) / np.sum(k * k, axis=0)

    return tmid, val * val


def median(x, y, dt=4.):
    """
    De-trend a light curve using a windowed median.

    """
    x, y = np.atleast_1d(x), np.atleast_1d(y)
    assert len(x) == len(y)
    r = np.empty(len(y))
    for i, t in enumerate(x):
        inds = (x >= t - 0.5 * dt) * (x <= t + 0.5 * dt)
        r[i] = np.median(y[inds])
    return r
