#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import, unicode_literals

__all__ = ["fit_trend", "untrend", "discontinuity_scalar", "median"]
__version__ = "0.0.1"
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
