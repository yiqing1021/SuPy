import numpy as np

from .._env import logger_supy
from ._atm import cal_cp




def cal_neutral(df_val, z_meas, h_sfc):
    """ Calculates the rows associated with neutral condition (threshold=0.01)


    Parameters
    ----------
    df_val: pd.DataFrame
        Index should be time with columns: 'H', 'USTAR', 'TA', 'RH', 'PA', 'WS'
    z_meas
        measurement height in m
    h_sfc
        vegetation height in m

    Returns
    -------
    ser_ws: pd.series
        observation time series of WS (Neutral conditions)
    ser_ustar: pd.series
        observation time series of u* (Neutral conditions)
    """

    # calculate Obukhov length
    ser_Lob = df_val.apply(
        lambda ser: cal_Lob(ser.H, ser.USTAR, ser.TA, ser.RH, ser.PA * 10), axis=1
    )

    # zero-plane displacement: estimated using rule f thumb `d=0.7*h_sfc`

    z_d = 0.7 * h_sfc

    if z_d >= z_meas:
        logger_supy.exception(
            "vegetation height is greater than measuring height. Please fix this before continuing . . ."
        )

    # calculate stability scale
    ser_zL = (z_meas - z_d) / ser_Lob

    # determine periods under quasi-neutral conditions
    limit_neutral = 0.01
    ind_neutral = ser_zL.between(-limit_neutral, limit_neutral)

    ind_neutral = ind_neutral[ind_neutral]

    df_sel = df_val.loc[ind_neutral.index, ["WS", "USTAR"]].dropna()
    ser_ustar = df_sel.USTAR
    ser_ws = df_sel.WS

    return ser_ws, ser_ustar


# Optimization for calculating z0 and d
def optimize_MO(df_val, z_meas, h_sfc):
    """Calculates surface roughness and zero plane displacement height.
    Refer to https://suews-parameters-docs.readthedocs.io/en/latest/steps/roughness-SuPy.html for example

    Parameters
    ----------
    df_val: pd.DataFrame
        Index should be time with columns: 'H', 'USTAR', 'TA', 'RH', 'PA', 'WS'
    z_meas
        measurement height in m
    h_sfc
        vegetation height in m

    Returns
    -------
    z0
        surface roughness
    d
        zero displacement height
    ser_ws: pd.series
        observation time series of WS (Neutral conditions)
    ser_ustar: pd.series
        observation time series of u* (Neutral conditions)
    """

    from platypus.core import Problem
    from platypus.types import Real, random
    from platypus.algorithms import NSGAIII

    # Calculates rows related to neutral conditions
    ser_ws, ser_ustar = cal_neutral(df_val, z_meas, h_sfc)

    # function to optimize
    def func_uz(params):
        z0 = params[0]
        d = params[1]
        z = z_meas
        k = 0.4
        uz = (ser_ustar / k) * np.log((z - d) / z0)  # logarithmic law

        o1 = abs(1 - np.std(uz) / np.std(ser_ws))  # objective 1: normalized STD
        # objective 2: normalized MAE
        o2 = np.mean(abs(uz - ser_ws)) / (np.mean(ser_ws))

        return [o1, o2], [uz.min(), d - z0]

    problem = Problem(2, 2, 2)
    problem.types[0] = Real(0, 10)  # bounds for first parameter (z0)
    problem.types[1] = Real(0, h_sfc)  # bounds for second parameter (zd)

    problem.constraints[0] = ">=0"  # constrain for first parameter
    problem.constraints[1] = ">=0"  # constrain for second parameter

    problem.function = func_uz
    random.seed(12345)
    algorithm = NSGAIII(problem, divisions_outer=50)
    algorithm.run(30000)

    z0s = []
    ds = []
    os1 = []
    os2 = []
    # getting the solution vaiables
    for s in algorithm.result:
        z0s.append(s.variables[0])
        ds.append(s.variables[1])
        os1.append(s.objectives[0])
        os2.append(s.objectives[1])
    # getting the solution associated with minimum obj2 (can be changed)
    idx = os2.index(min(os2, key=lambda x: abs(x - np.mean(os2))))
    z0 = z0s[idx]
    d = ds[idx]

    return z0, d, ser_ws, ser_ustar
