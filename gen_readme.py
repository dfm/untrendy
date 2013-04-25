#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Automatically construct the README for untrendy by scraping the docstrings
from the module.

"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import inspect
import textwrap

import untrendy


def get_function_decription(f):
    doc = textwrap.dedent(f.__doc__)
    name = "*untrendy.*\ **{0}** (".format(f.__name__)
    argspec = inspect.getargspec(f)

    args = argspec.args[::-1]
    defaults = argspec.defaults[::-1]
    fmt_args = []

    if argspec.keywords:
        fmt_args += ["**{0}".format(argspec.keywords)]

    for i, v in enumerate(defaults):
        fmt_args.append("{0}={1}".format(args[i], v))
    fmt_args += args[i + 1:]

    name += ", ".join(["``{0}``".format(a) for a in fmt_args[::-1]]) + ")"
    name = textwrap.fill(name, 80)
    return (name + "\n" + doc).strip()


ns = {}
exec open("bin/untrend").read() in ns

parts = [untrendy.__doc__.strip(), ns["__doc__"].strip(), "API\n---",
         "Fit the trend\n+++++++++++++",
         get_function_decription(untrendy.fit_trend),
         "Remove the trend\n++++++++++++++++",
         get_function_decription(untrendy.untrend)]

with open("README.rst", "w") as f:
    f.write("\n\n".join(parts))
