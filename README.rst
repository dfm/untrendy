￭ Square²
=========

**That's a seriously untrendy verb™**

    "De-trending Kepler light curves in style."

Installation
------------

**Square** depends on ``numpy`` and ``scipy`` so make sure that you install
those first. You'll also need a proper Fortran compiler (say ``gfortran``) so
have one of those too. Then, you can install using ``pip``:

::

    pip install square

Usage
-----

**Square** is really complicated. It has *one* function and about *200 lines
of code (including documentation)*. It mostly runs on love and magic (more
complete details `are available <http://dan.iel.fm/square>`_ if you want).

Let's say that you have a light curve with time samples ``t``, flux
measurements ``f`` and uncertainties ``sigma``. You can simply run:

::

    import square
    spline = square.detrend(t, f, sigma)

to find a robust estimate of the global trends of the time series. The default
settings are tuned to work well for finding the "out-of-transit" trends in
Kepler data but a detailed description of the options is given in `the
documentation <http://dan.iel.fm/square>`_. As the name suggests, ``spline``
is a callable cubic spline representation of the trends.To de-trend your data,
just do something like:

::

    factor = spline(t)
    f /= factor
    sigma /= factor

Notes
-----

1. The spline sometimes goes to hell in regions where you don't have any
   samples so be careful with that.
2. This whole procedure introduces correlated errors. You've been warned.

License
-------

Copyright 2013 Dan Foreman-Mackey & contributors.

Bart is free software made available under the MIT License. For details see
`the LICENSE file <https://raw.github.com/dfm/square/master/LICENSE.rst>`_.
