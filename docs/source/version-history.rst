.. _version_history:

Version History
================

.. _new_latest:

Version 2021.11.20
------------------------------


- **New**

  1. Added option ``pressure_level`` in `gen_forcing_era5`.


- **Improvement**

  None.

- **Changes**

  None.


- **Fix**

  1. Fixed an issue in generating ERA5 forcing due to xarray update in ``merge``.


- **Known issue**

  1. ESTM is not supported yet.
  2. BLUEWS, a CBL modules in SUEWS, is not supported yet.
  3. Simulation in parallel mode is NOT supported on Windows
     due to system limitation.



Version 2021.7.22
------------------------------


- **New**

  None.


- **Improvement**

  None.

- **Changes**

  None.


- **Fix**

  1. Fixed an issue in loading parameter table caused by recent update of ``pandas`` to ``1.3.x``.
  2. Fixed an issue in ERA5 download due to renaming of ``orography`` to ``geopotential``.


- **Known issue**

  1. ESTM is not supported yet.
  2. BLUEWS, a CBL modules in SUEWS, is not supported yet.
  3. Simulation in parallel mode is NOT supported on Windows
     due to system limitation.


Version 2021.5.26
------------------------------


- **New**

  1. Update supy-driver to `2021a` iteration.


- **Improvement**

  1. a new method for calculating roughness length for momentum and displacement height (`roughlenmommethod=4`) based on plan area index as illustrated in figure 1a of GO99.

- **Changes**

  None.


- **Fix**

  None.


- **Known issue**

  1. ESTM is not supported yet.
  2. BLUEWS, a CBL modules in SUEWS, is not supported yet.
  3. Simulation in parallel mode is NOT supported on Windows
     due to system limitation.





Version 2020.11.3
------------------------------


- **New**

  1. Update supy-driver to `2020b` iteration.
  2. Add function for plotting RSL variables `supy.util.plot_rsl`.


- **Improvement**

  1. The RSL related functions are more robust in dealing with broader urban morphology settings.
  2. Internal changes to conform with recent upgrades in `pandas`.

- **Changes**

  None.


- **Fix**

  1. Fix an issue in `supy.util.read_forcing` that improper resampling could be conducted
  if input temporal resolution is the same as the desirable resampling time step `tstep_mod`.


- **Known issue**

  1. ESTM is not supported yet.
  2. BLUEWS, a CBL modules in SUEWS, is not supported yet.
  3. Simulation in parallel mode is NOT supported on Windows
     due to system limitation.




Version 2020.5.29
------------------------------


- **New**

  1. Update supy-driver to `2020a` iteration.
  2. Add function for plotting RSL variables `supy.util.plot_rsl`.


- **Improvement**

  None.

- **Changes**

  None.


- **Fix**

  1. Fix the humidity variable in ERA5-based forcing generation.
  2. Fix the impact study tutorial.


- **Known issue**

  1. ESTM is not supported yet.
  2. BLUEWS, a CBL modules in SUEWS, is not supported yet.
  3. Simulation in parallel mode is NOT supported on Windows
     due to system limitation.






Version 2020.2.2
------------------------------


- **New**

  1. A checker to validate input `DataFrame`s.
     See option `check_input` in `run_supy`.

  2. Utilities to generate forcing data using ERA-5 data.
     See `download_era5` and `gen_forcing_era5`.

- **Improvement**

  1. Improved performance of the parallel mode.

- **Changes**

  None.


- **Fix**

  None.


- **Known issue**

  1. ESTM is not supported yet.
  2. BLUEWS, a CBL modules in SUEWS, is not supported yet.
  3. Simulation in parallel mode is NOT supported on Windows
     due to system limitation.



.. _new_20190829:


Version 2019.8.29
------------------------------


- **New**

  1. added WRF-SUEWS related functions.
  2. added `diagnostics of canyon profiles <https://suews-docs.readthedocs.io/en/test-dev/parameterisations-and-sub-models.html#wind-temperature-and-humidity-profiles-in-the-roughness-sublayer>`_.

- **Improvement**

  None.

- **Changes**

  1. synchronised with v2019a interface: minimum supy_driver v2019a2.


- **Fix**

  None.


- **Known issue**

  1. ESTM is not supported yet.
  2. BLUEWS, a CBL modules in SUEWS, is not supported yet.
  3. Performance in parallel mode can be worse than
     serial mode sometimes due to heavy (de)-serialisation loads.



.. _new_20190717:

Version 2019.7.17
------------------------------


- **New**

  1. added OHM related functions.
  2. added surface conductance related functions.

- **Improvement**

  None.

- **Changes**

  None.


- **Fix**

  1. Fixed a bug in unit conversion for TMY data generation.


- **Known issue**

  ESTM is not supported yet.




.. _new_20190608:

Version 2019.6.8
------------------------------


- **New**

  None.

- **Improvement**

  None.

- **Changes**

  None.


- **Fix**

  1. Fixed a bug in rescaling Kdown when loading forcing data.


- **Known issue**

  ESTM is not supported yet.



.. _new_20190528:

Version 2019.5.28
------------------------------

Spring house cleaning with long-await command line tools (more on the way!).

- **New**

  1. Added version info function: `show_version`.
  2. Added command line tools:

    - `suews-run`: SuPy wrapper to mimic SUEWS-binary-based simulation.

    - `suews-convert`: convert input tables from older versions to newer ones (one-way only).


- **Improvement**

  None.

- **Changes**

  None.


- **Fix**

  1. Fixed a bug in writing out multi-grid output files
  caused by incorrect dropping of temporal information by pandas .

- **Known issue**

  ESTM is not supported yet.



.. _new_20190429:

Version 2019.4.29
------------------------------

Parallel run.

- **New**

  Added support for parallel run on the fly.

- **Improvement**

  None.

- **Changes**

  None.


- **Fix**

  None.

- **Known issue**

  None


.. _new_20190417:

Version 2019.4.17
------------------------------

UMEP compatibility tweaks.

- **New**

  None.

- **Improvement**

  None.

- **Changes**

  `problems.txt` will be written out in addition to the console error message similarly as SUEWS binary.


- **Fix**

  Incorrect caching of input libraries.

- **Known issue**

  None



.. _new_20190415:

Version 2019.4.15
------------------------------

ERA-5 download.

- **New**

  Added experimental support for downloading and processing ERA-5 data to force supy simulations.

- **Improvement**

  Improved compatibility with earlier `pandas` version in resampling output.

- **Changes**

  None.


- **Fix**

  None.

- **Known issue**

  None


.. _new_20190321:

Version 2019.3.21
------------------------------

TMY generation.

- **New**

  Added preliminary support for generating TMY dataset with SuPy output.

- **Improvement**

  None.

- **Changes**

  None.


- **Fix**

  None.

- **Known issue**

  None

.. _new_20190314:

Version 2019.3.14
------------------------------

This release improved memory usage.

- **New**

  None.

- **Improvement**

  Optimised memory consumption for longterm simulations.


- **Changes**

  None.


- **Fix**

  None.

- **Known issue**

  None


.. _new_20190225:

Version 2019.2.25
------------------------------

This release dropped support for Python 3.5 and below.

- **New**

  None.

- **Improvement**

  None.


- **Changes**

  Dropped support for Python 3.5 and below.


- **Fix**

  None.

- **Known issue**

  None


.. _new_20190224:

Version 2019.2.24
------------------------------

This release added the ability to save output files.

- **New**

  1. Added support to save output files. See: :py:func:`supy.save_supy`
  2. Added support to initialise SuPy from saved ``df_state.csv``. See: :py:func:`supy.init_supy`

- **Improvement**

  None.


- **Changes**

  None.


- **Fix**

  None.

- **Known issue**

  None


.. _new_20190219:

Version 2019.2.19
------------------------------

This is a release that improved the exception handling due to
fatal error in `supy_driver`.

- **New**

  Added support to handle python kernel crash caused by
  fatal error in `supy_driver` kernel; so python kernel won't
  crash any more even `supy_driver` is stopped.

- **Improvement**

  None.


- **Changes**

  None


- **Fix**

  None.

- **Known issue**

  None


.. _new_20190208:

Version 2019.2.8
------------------------------

This is a release that fixes recent bugs found in SUEWS that
may lead to abnormal simulation results of storage heat flux,
in particular when `SnowUse` is enabled (i.e., `snowuse=1`).

- **New**

  None.

- **Improvement**

  Improved the performance in loading
  initial model state from a large number of grids (>1k)


- **Changes**

  Updated `SampleRun` dataset by:
  1. setting surface fractions (`sfr`) to a
  more realistic value based on London KCL case;
  2. enabling snow module (`snowuse=1`).


- **Fix**

  1. Fixed a bug in the calculation of storage heat flux.
  2. Fixed a bug in loading ``popdens`` for calculating anthropogenic heat flux.

- **Known issue**

  None

.. _new_20190101:

Version 2019.1.1 (preview release, 01 Jan 2019)
------------------------------


- **New**

  #. Slimmed the output groups by excluding unsupported :term:`ESTM` results

  #. SuPy documentation

    * Key IO data structures documented:

     - `df_output_var` (:issue:`9`)
     - `df_state_var` (:issue:`8`)
     - `df_forcing_var` (:issue:`7`)


    * Tutorial of parallel SuPy simulations for impact studies



- **Improvement**

  #. Improved calculation of OHM-related radiation terms


- **Changes**

  None.


- **Fix**

  None

- **Known issue**

  None

.. _new_20181215:

Version 2018.12.15 (internal test release in December 2018)
------------------------------=====


- **New**

  #. Preview release of SuPy based on the computation kernel of SUEWS 2018b


- **Improvement**

  #. Improved calculation of OHM-related radiation terms


- **Changes**

  None.


- **Fix**

  None

- **Known issue**

  #. The heat storage modules AnOHM and ESTM are not supported yet.





