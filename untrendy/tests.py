#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__all__ = []


import os
import pyfits
from .untrendy import detrend


def _load_kepler_lc(name):
    fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_data",
                      name)
    with pyfits.open(fn) as f:
        time, flux, ferr = (f[1].data["TIME"], f[1].data["SAP_FLUX"],
                            f[1].data["SAP_FLUX_ERR"])

    return time, flux, ferr


def test_SW_condition():
    time, flux, ferr = _load_kepler_lc("kplr010874614-2009131105131_llc.fits")
    detrend(time, flux, ferr)
