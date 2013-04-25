#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Untrendy
========

In an age where studying exoplanets is just the hippest thing ever, sometimes
it's good to step out of line and be a little untrendy! This library is a set
of hacks that can robustly remove the out-of-transit trends in light curve
data.

Installation
------------

**Untrendy** depends on ``numpy`` and ``scipy`` so make sure that you install
those first. Then, you can install using ``pip``:

::

    pip install untrendy

Usage
-----

**Untrendy** is really complicated. It has approximately *one* function and
about *200 lines of code (including documentation)*. It mostly runs on love
and magic (more complete details are given below if you want).

Let's say that you have a light curve with time samples ``t``, flux
measurements ``f`` and uncertainties ``sigma``. You can simply run:

.. code-block:: python

    import untrendy
    f_detrend, sigma_detrend = untrendy.untrend(t, f, sigma)

to find a robust estimate of the global trends of the time series and remove
it. The default settings are tuned to work well for finding the
"out-of-transit" trends in Kepler data but a detailed description of the
options is listed below. You can also just fit for the trends and get a
callable representation of the trend:

.. code-block:: python

    trend = untrendy.fit_trend(t, f, ferr)

In this case, you can find the background level at some time ``t0`` by calling
the function:

.. code-block:: python

    bkg = trend(t0)

Notes
-----

1. The spline sometimes goes to hell in regions where you don't have any
   samples so be careful with that.
2. This whole procedure introduces correlated errors. You've been warned.

"""

from __future__ import absolute_import, unicode_literals

__all__ = ["fit_trend", "untrend", "discontinuity_scalar", "median"]
__version__ = "0.0.2"
__author__ = "Dan Foreman-Mackey (danfm@nyu.edu)"
__copyright__ = "Copyright 2013 Daniel Foreman-Mackey"
__contributors__ = []

from .untrendy import fit_trend, untrend, discontinuity_scalar, median


def test():
    from inspect import getmembers, isfunction
    from . import tests

    for o in getmembers(tests):
        if isfunction(o[1]) and o[0].startswith("test"):
            print("{0} ...".format(o[0]))
            try:
                o[1]()
            except Exception as e:
                print("Failed with:\n    {0.__class__.__name__}: {0}"
                      .format(e))
            else:
                print("    Passed.")
