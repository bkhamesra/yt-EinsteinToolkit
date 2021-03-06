.. _cosmology-calculator:

Cosmology Calculator
====================

The cosmology calculator can be used to calculate cosmological distances and
times given a set of cosmological parameters.  A cosmological dataset, `ds`,
will automatically have a cosmology calculator configured with the correct
parameters associated with it as `ds.cosmology`.  A standalone
:class:`~yt.utilities.cosmology.Cosmology` calculator object can be created
in the following way:

.. code-block:: python

   from yt.utilities.cosmology import Cosmology

   co = Cosmology(hubble_constant=0.7, omega_matter=0.3,
                  omega_lambda=0.7, omega_curvature=0.0)

Once created, various distance calculations as well as conversions between
redshift and time are available:

.. notebook-cell::

   from yt.utilities.cosmology import Cosmology

   co = Cosmology(hubble_constant=0.7, omega_matter=0.3,
                  omega_lambda=0.7, omega_curvature=0.0)

   # Hubble distance (c / h)
   print("hubble distance", co.hubble_distance())

   # distance from z = 0 to 0.5
   print("comoving radial distance", co.comoving_radial_distance(0, 0.5).in_units("Mpc/h"))

   # transverse distance
   print("transverse distance", co.comoving_transverse_distance(0, 0.5).in_units("Mpc/h"))

   # comoving volume
   print("comoving volume", co.comoving_volume(0, 0.5).in_units("Gpc**3"))

   # angulare diameter distance
   print("angular diameter distance", co.angular_diameter_distance(0, 0.5).in_units("Mpc/h"))

   # angular scale
   print("angular scale", co.angular_scale(0, 0.5).in_units("Mpc/degree"))

   # luminosity distance
   print("luminosity distance", co.luminosity_distance(0, 0.5).in_units("Mpc/h"))

   # time between two redshifts
   print("lookback time", co.lookback_time(0, 0.5).in_units("Gyr"))

   # age of the Universe at a given redshift
   print("hubble time", co.hubble_time(0).in_units("Gyr"))

   # critical density
   print("critical density", co.critical_density(0))

   # Hubble parameter at a given redshift
   print("hubble parameter", co.hubble_parameter(0).in_units("km/s/Mpc"))

   # convert time after Big Bang to redshift
   my_t = co.quan(8, "Gyr")
   print("z from t", co.z_from_t(my_t))

   # convert redshift to time after Big Bang (same as Hubble time)
   print("t from z", co.t_from_z(0.5).in_units("Gyr"))

Note, that all distances returned are comoving distances.  All of the above
functions accept scalar values and arrays.  The helper functions, `co.quan`
and `co.arr` exist to create unitful `YTQuantities` and `YTArray` with the
unit registry of the cosmology calculator.  For more information on the usage
and meaning of each calculation, consult the reference documentation at
:ref:`cosmology-calculator-ref`.
