#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import os
import sys

from numpy.distutils.core import setup, Extension


if sys.argv[-1] == "publish":
    os.system("python setup.py sdist upload")
    sys.exit()

# First, make sure that the f2py interfaces exist.
interface_exists = os.path.exists("untrendy/untrendy.pyf")
if "interface" in sys.argv or not interface_exists:
    # Generate the Fortran signature/interface.
    cmd = ("cd untrendy;f2py discontinuities.f90 -m _untrendy -h untrendy.pyf"
           " --overwrite-signature")
    os.system(cmd)
    if "interface" in sys.argv:
        sys.exit(0)

# Define the Fortran extension.
untrendy = Extension("untrendy._untrendy", ["untrendy/untrendy.pyf",
                                            "untrendy/discontinuities.f90"])

# Get the version number.
vre = re.compile("__version__ = \"(.*?)\"")
m = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "untrendy", "__init__.py")).read()
version = vre.findall(m)[0]

setup(
    name="untrendy",
    url="https://github.com/dfm/untrendy",
    version=version,
    author="Dan Foreman-Mackey",
    author_email="danfm@nyu.edu",
    description="De-trending Kepler light curves in style",
    long_description=open("README.rst").read(),
    packages=["untrendy"],
    scripts=["bin/untren"],
    package_data={"": ["README.rst", "LICENSE.rst"]},
    ext_modules=[untrendy],
    classifiers=[
        # "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ],
)
