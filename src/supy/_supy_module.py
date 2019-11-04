# ###########################################################################
# SuPy: SUEWS for Python
#
# Authors:
# Ting Sun, ting.sun@reading.ac.uk
#
# History:
# 20 Jan 2018: first alpha release
# 01 Feb 2018: performance improvement
# 03 Feb 2018: improvement in output processing
# 08 Mar 2018: pypi packaging
# 04 Oct 2018: overhaul of structure
# 05 Oct 2018: added sample run data
# 28 Apr 2019: added support for parallel run
###########################################################################


import logging
import multiprocessing
import time
# from multiprocessing import Pool, cpu_count, freeze_support

import os
import sys
# import functools
from pathlib import Path
from typing import Tuple

import pandas
import pathlib


import numpy as np
import pandas as pd

from ._env import path_supy_module
from ._load import (load_InitialCond_grid_df,
                    load_SUEWS_dict_ModConfig,
                    # load_SUEWS_Forcing_ESTM_df_raw,
                    load_SUEWS_Forcing_met_df_raw,
                    load_df_state,
                    resample_forcing_met,
                    resample_linear)
from ._post import pack_df_output, pack_df_output_array, pack_df_state
from ._run import run_supy_ser, run_supy_par
from ._save import get_save_info, save_df_output, save_df_state, save_initcond_nml
from ._check import check_forcing, check_state
from ._env import logger_supy

# set up logging module
logger_supy.setLevel(logging.INFO)


##############################################################################
# 1. compact wrapper for loading SUEWS settings
# @functools.lru_cache(maxsize=16)
def init_supy(path_init: str, force_reload=True) -> pd.DataFrame:
    '''Initialise supy by loading initial model states.

    Parameters
    ----------
    path_init : str
        Path to a file that can initialise SuPy, which can be either of the follows:
            * SUEWS :ref:`RunControl.nml<suews:RunControl.nml>`: a namelist file for SUEWS configurations
            * SuPy `df_state.csv`: a CSV file including model states produced by a SuPy run via :py:func:`supy.save_supy`

    force_reload: boolean, optional
        Flag to force reload all initialisation files by clearing all cached states, with default value `True` (i.e., force reload all files).
        Note: If the number of simulation grids is large (e.g., > 100), `force_reload=False` is strongly recommended for better performance.


    Returns
    -------
    df_state_init: pandas.DataFrame
        Initial model states.
        See `df_state_var` for details.

    Examples
    --------
    1. Use :ref:`RunControl.nml<suews:RunControl.nml>` to initialise SuPy

    >>> path_init = "~/SUEWS_sims/RunControl.nml"
    >>> df_state_init = supy.init_supy(path_init)

    2. Use ``df_state.csv`` to initialise SuPy

    >>> path_init = "~/SuPy_res/df_state_test.csv"
    >>> df_state_init = supy.init_supy(path_init)

    '''

    try:
        path_init_x = Path(path_init).expanduser().resolve()
    except FileNotFoundError:
        logger_supy.exception(f'{path_init_x} does not exists!')
    else:
        if path_init_x.suffix == '.nml':
            # SUEWS `RunControl.nml`:
            df_state_init = load_InitialCond_grid_df(
                path_init_x, force_reload=force_reload)
        elif path_init_x.suffix == '.csv':
            # SuPy `df_state.csv`:
            df_state_init = load_df_state(path_init_x)
        else:
            logger_supy.critical(
                f'{path_init_x} is NOT a valid file to initialise SuPy!')
            sys.exit()
        return df_state_init


# # TODO:
# def load_forcing(path_pattern: str, grid: int = 0) -> pd.DataFrame:
#     pass


# TODO:
# to be superseded by a more generic wrapper: load_forcing
def load_forcing_grid(path_runcontrol: str, grid: int) -> pd.DataFrame:
    '''Load forcing data for a specific grid included in the index of `df_state_init </data-structure/supy-io.ipynb#df_state_init:-model-initial-states>`.

    Parameters
    ----------
    path_runcontrol : str
        Path to SUEWS :ref:`RunControl.nml <suews:RunControl.nml>`
    grid : int
        Grid number

    Returns
    -------
    df_forcing: pandas.DataFrame
        Forcing data. See `df_forcing_var` for details.

    Examples
    --------
    >>> path_runcontrol = "~/SUEWS_sims/RunControl.nml"  # a valid path to `RunControl.nml`
    >>> df_state_init = supy.init_supy(path_runcontrol) # get `df_state_init`
    >>> grid = df_state_init.index[0] # first grid number included in `df_state_init`
    >>> df_forcing = supy.load_forcing_grid(path_runcontrol, grid) # get df_forcing


    '''

    try:
        path_runcontrol = Path(path_runcontrol).expanduser().resolve()
    except FileNotFoundError:
        logger_supy.exception(f'{path_runcontrol} does not exists!')
    else:
        dict_mod_cfg = load_SUEWS_dict_ModConfig(path_runcontrol)
        df_state_init = init_supy(path_runcontrol)

        # load setting variables from dict_mod_cfg
        (
            filecode,
            kdownzen,
            tstep_met_in,
            tstep_ESTM_in,
            multiplemetfiles,
            multipleestmfiles,
            dir_input_cfg
        ) = (dict_mod_cfg[x] for x in
             [
            'filecode',
            'kdownzen',
            'resolutionfilesin',
            'resolutionfilesinestm',
            'multiplemetfiles',
            'multipleestmfiles',
            'fileinputpath'
        ]
        )
        tstep_mod, lat, lon, alt, timezone = df_state_init.loc[
            grid,
            [(x, '0') for x in ['tstep', 'lat', 'lng', 'alt', 'timezone']]
        ].values

        path_site = path_runcontrol.parent
        path_input = path_site / dict_mod_cfg['fileinputpath']

        # load raw data
        # met forcing
        df_forcing_met = load_SUEWS_Forcing_met_df_raw(
            path_input, filecode, grid, tstep_met_in, multiplemetfiles)

        # resample raw data from tstep_in to tstep_mod
        df_forcing_met_tstep = resample_forcing_met(
            df_forcing_met, tstep_met_in, tstep_mod,
            lat, lon, alt, timezone, kdownzen)

        # coerced precision here to prevent numerical errors inside Fortran
        df_forcing = df_forcing_met_tstep.round(10)

        # new columns for later use in main calculation
        df_forcing[['iy', 'id', 'it', 'imin']] = df_forcing[[
            'iy', 'id', 'it', 'imin']].astype(np.int64)

    return df_forcing


# load sample data for quickly starting a demo run
# TODO: to deprecate this by renaming for case consistency: load_SampleData-->load_sample_data
def load_SampleData() -> Tuple[pandas.DataFrame, pandas.DataFrame]:
    '''Load sample data for quickly starting a demo run.

    Returns
    -------
    df_state_init, df_forcing: Tuple[pandas.DataFrame, pandas.DataFrame]
        - df_state_init: `initial model states <df_state_var>`
        - df_forcing: `forcing data <df_forcing_var>`

    Examples
    --------

    >>> df_state_init, df_forcing = supy.load_SampleData()

    '''

    path_SampleData = Path(path_supy_module) / 'sample_run'
    path_runcontrol = path_SampleData / 'RunControl.nml'
    df_state_init = init_supy(path_runcontrol, force_reload=False)
    # path_input = path_runcontrol.parent / ser_mod_cfg['fileinputpath']
    df_forcing = load_forcing_grid(
        path_runcontrol,
        df_state_init.index[0]
    )
    return df_state_init, df_forcing

# input processing code end here
##############################################################################


##############################################################################
# 2. compact wrapper for running a whole simulation
# # main calculation
# input as DataFrame
def run_supy(
        df_forcing: pandas.DataFrame,
        df_state_init: pandas.DataFrame,
        save_state=False,
        n_yr=10,
        logging_level=logging.INFO,) -> Tuple[pandas.DataFrame, pandas.DataFrame]:
    '''Perform supy simulation.

    Parameters
    ----------
    df_forcing : pandas.DataFrame
        forcing data for all grids in `df_state_init`.
    df_state_init : pandas.DataFrame
        initial model states;
        or a collection of model states with multiple timestamps, whose last temporal record will be used as the initial model states.
    save_state : bool, optional
        flag for saving model states at each time step, which can be useful in diagnosing model runtime performance or performing a restart run.
        (the default is False, which instructs supy not to save runtime model states).
    n_yr : int, optional
        chunk size (`n_yr` years) to split simulation periods so memory usage can be reduced.
        (the default is 10, which implies 10-year forcing chunks used in simulations).
    logging_level: logging level
        one of these values [50 (CRITICAL), 40 (ERROR), 30 (WARNING), 20 (INFO), 10 (DEBUG)].
        A lower value informs SuPy for more verbose logging info.

    Returns
    -------
    df_output, df_state_final : Tuple[pandas.DataFrame, pandas.DataFrame]
        - df_output: `output results <df_output_var>`
        - df_state_final: `final model states <df_state_var>`

    Examples
    --------

    >>> df_output, df_state_final = supy.run_supy(df_forcing, df_state_init)


    '''

    # set up a timer for simulation time
    start = time.time()

    # adjust logging level
    logger_supy.setLevel(logging_level)

    # save df_init without changing its original data
    # df.copy() in pandas works as a standard python deepcopy
    # df_init = df_state_init.copy()

    # print some diagnostic info
    logger_supy.info(f'====================')
    logger_supy.info(f'Simulation period:')
    logger_supy.info(f'  Start: {df_forcing.index[0]}')
    logger_supy.info(f'  End: {df_forcing.index[-1]}')
    logger_supy.info('')
    list_grid = df_state_init.index.get_level_values('grid').unique()
    n_grid = list_grid.size
    logger_supy.info(f'No. of grids: {n_grid}')

    if n_grid > 1 and os.name != 'nt':
        logger_supy.info(f'SuPy is running in parallel mode')
        df_output, df_state_final = run_supy_par(
            df_forcing, df_state_init, save_state, n_yr)
    else:
        logger_supy.info(f'SuPy is running in serial mode')
        df_output, df_state_final = run_supy_ser(
            df_forcing, df_state_init, save_state, n_yr)

    # show simulation time
    end = time.time()
    logger_supy.info(f'Execution time: {(end - start):.1f} s')
    logger_supy.info(f'====================\n')

    return df_output, df_state_final


##############################################################################
# 3. save results of a supy run
def save_supy(
        df_output: pandas.DataFrame,
        df_state_final: pandas.DataFrame,
        freq_s: int = 3600,
        site: str = '',
        path_dir_save: str = Path('.'),
        path_runcontrol: str = None,
        logging_level=50) -> list:
    '''Save SuPy run results to files

    Parameters
    ----------
    df_output : pandas.DataFrame
        DataFrame of output
    df_state_final : pandas.DataFrame
        DataFrame of final model states
    freq_s : int, optional
        Output frequency in seconds (the default is 3600, which indicates hourly output)
    site : str, optional
        Site identifier (the default is '', which indicates site identifier will be left empty)
    path_dir_save : str, optional
        Path to directory to saving the files (the default is Path('.'), which indicates the current working directory)
    path_runcontrol : str, optional
        Path to SUEWS :ref:`RunControl.nml <suews:RunControl.nml>`, which, if set, will be preferably used to derive `freq_s`, `site` and `path_dir_save`.
        (the default is None, which is unset)
    logging_level: logging level
        one of these values [50 (CRITICAL), 40 (ERROR), 30 (WARNING), 20 (INFO), 10 (DEBUG)].
        A lower value informs SuPy for more verbose logging info.

    Returns
    -------
    list
        a list of paths of saved files

    Examples
    --------
    1. save results of a supy run to the current working directory with default settings

    >>> list_path_save = supy.save_supy(df_output, df_state_final)


    2. save results according to settings in :ref:`RunControl.nml <suews:RunControl.nml>`

    >>> list_path_save = supy.save_supy(df_output, df_state_final, path_runcontrol='path/to/RunControl.nml')


    3. save results of a supy run at resampling frequency of 1800 s (i.e., half-hourly results) under the site code ``Test`` to a customised location 'path/to/some/dir'

    >>> list_path_save = supy.save_supy(df_output, df_state_final, freq_s=1800, site='Test', path_dir_save='path/to/some/dir')
    '''
    # adjust logging level
    logger_supy.setLevel(logging_level)

    # get necessary information for saving procedure
    if path_runcontrol is not None:
        freq_s, path_dir_save, site = get_save_info(path_runcontrol)

    # save df_output to several files
    list_path_save = save_df_output(df_output, freq_s, site, path_dir_save)

    # save df_state
    if path_runcontrol is not None:
        # save as nml as SUEWS binary
        list_path_nml = save_initcond_nml(df_state_final, site, path_dir_save)
        list_path_save = list_path_save+list_path_nml
    else:
        # save as supy csv for later use
        path_state_save = save_df_state(df_state_final, site, path_dir_save)
        # update list_path_save
        list_path_save.append(path_state_save)

    return list_path_save
