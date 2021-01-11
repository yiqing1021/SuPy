# from scipy.optimize import least_squares
import numpy as np
import pandas as pd


# atmospheric related utilities
def cal_des_dta(ta, pa, dta=1.0):
    """Calculate slope of es(Ta), i.e., saturation evaporation pressure `es` as function of air temperature `ta [K]`

    Parameters
    ----------
    ta : numeric
        Air temperature [K]
    pa : numeric
        Air pressure [Pa]
    dta : float, optional
        change in ta for calculating that in es, by default 1.0 K
    """
    from atmosp import calculate as ac

    des = ac("es", p=pa, T=ta + dta / 2) - ac("es", p=pa, T=ta - dta / 2)
    des_dta = des / dta
    try:
        # try to pack as Series
        des_dta = pd.Series(des_dta, index=ta.index)
    except AttributeError as ex:
        print(ex, "cannot pack into pd.Series")
        pass
    return des_dta


def cal_rs_obs(qh, qe, ta, rh, pa, ra):
    """Calculate surface resistance based on observations, notably turbulent fluxes.

    Parameters
    ----------
    qh : numeric
        sensible heat flux [W m-2]
    qe : numeric
        latent heat flux [W m-2]
    ta : numeric
        air temperature [degC]
    rh : numeric
        relative humidity [%]
    pa : numeric
        air pressure [Pa]
    ra : numeric
        aerodynamic resistance [m s-1]

    Returns
    -------
    numeric
        Surface resistance based on observations [s m-1]
    """
    from atmosp import calculate as ac

    # psychrometric constant [Pa K-1] as a function of air pressure
    ser_gamma = 0.665e-3 * pa

    # air density [kg m-3]
    val_rho = 1.27

    # heat capacity of air [J kg-1 K-1]
    val_cp = 1005

    # convert temp from C to K
    ta_K = ta + 273.15

    # slope of es(Ta) curve at Ta
    ser_des_dTa = cal_des_dta(ta_K, pa, dta=1.0)
    #
    arr_e = ac("e", p=pa, T=ta_K, RH=rh)
    arr_es = ac("es", p=pa, T=ta_K)
    arr_vpd = arr_es - arr_e
    #
    ser_rs_1 = (ser_des_dTa / ser_gamma) * (qh / qe - 1) * ra
    ser_rs_2 = val_rho * val_cp * arr_vpd / (ser_gamma * qe)
    ser_rs = ser_rs_1 + ser_rs_2

    try:
        # try to pack as Series
        ser_rs = pd.Series(ser_rs, index=ta_K.index)
    except AttributeError as ex:
        print(ex, "cannot pack into pd.Series")
        pass

    return ser_rs


def cal_gs_obs(qh, qe, ta, rh, pa, ra):
    """Calculate surface conductance based on observations, notably turbulent fluxes.

    Parameters
    ----------
    qh : numeric
        Sensible heat flux [W m-2]
    qe : numeric
        Latent heat flux [W m-2]
    ta : numeric
        Air temperature [degC]
    rh : numeric
        Relative humidity [%]
    pa : numeric
        Air pressure [Pa]

    Returns
    -------
    numeric
        Surface conductance based on observations [mm s-1]
    """
    rs_obs = cal_rs_obs(qh, qe, ta, rh, pa, ra)
    gs_obs = 1e3 / rs_obs
    return gs_obs


def cal_g_lai(lai, g_lai, lai_max):
    """Calculate LAI-related correction coefficient for surface conductance.

    Parameters
    ----------
    lai : numeric
        Leaf area index [m2 m-2]
    g_lai : numeric
        LAI-related correction parameter [-]
    lai_max : numeric
        Maximum LAI [m2 m-2]

    Returns
    -------
    numeric
        LAI-related correction coefficient [-]
    """
    g_lai = lai / lai_max * g_lai
    return g_lai


def cal_g_kd(kd, g2, kd_max=1200.0):
    """Calculate solar radiation-related correction coefficient for surface conductance.

    Parameters
    ----------
    kd : numeric
        Incoming solar radiation [W m-2]
    g2 : numeric
        Solar radiation-related correction parameter [-]
    kd_max : numeric, optional
        Maximum incoming solar radiation [W m-2], by default 1200.

    Returns
    -------
    numeric
        Solar radiation-related correction coefficient [-]
    """
    g_kd_nom = kd / (g2 + kd)
    g_kd_denom = kd_max / (g2 + kd_max)
    g_kd = g_kd_nom / g_kd_denom
    return g_kd


def cal_g_dq(dq, g_q1, g_q2):
    """Calculate air humidity-related correction coefficient for surface conductance.

    Parameters
    ----------
    dq : numeric
        Specific humidity deficit [g kg-1]
    g_q1 : numeric
        Specific humidity-related correction parameter [-]
    g_q2 : numeric
        Specific humidity-related correction parameter [-]

    Returns
    -------
    numeric
        Air humidity-related correction coefficient
    """
    g_dq = g_q1 + (1 - g_q1) * g_q2 ** dq
    return g_dq


def cal_g_ta(ta_c, g_ta, tl=-10.0, th=55.0):
    """Calculate air temperature-related correction coefficient for surface conductance.

    Parameters
    ----------
    ta_c : numeric
        Air temperature [degC]
    g_ta : numeric
        Air temperature-related correction parameter
    tl : numeric, optional
        Low temperature limit [degC], by default -10.
    th : numeric, optional
        High temperature limit [degC], by default 55.

    Returns
    -------
    numeric
        Air temperature-related correction coefficient
    """

    tc = (th - g_ta) / (g_ta - tl)
    # set a threshold for avoiding numeric difficulty
    tc = np.min([tc, 20])
    # g_ta = ((ta_c-tl)*(th-ta_c)**tc)/((g_ta-tl)*(th-g_ta)**tc)
    g_ta_nom = (ta_c - tl) * np.power((th - ta_c), tc)
    g_ta_denom = (g_ta - tl) * np.power((th - g_ta), tc)
    g_ta = g_ta_nom / g_ta_denom

    return g_ta


def cal_g_smd(smd, g_smd, wp_smd):
    """Calculate soil moisture-related correction coefficient for surface conductance.

    Parameters
    ----------
    smd : numeric
        Soil moisture deficit [mm].
    g_smd : numeric
        Soil moisture-related correction parameter.
    wp_smd : numeric, optional
        Wilting point indicated by soil moisture deficit [mm]

    Returns
    -------
    numeric
        Soil moisture-related correction coefficient
    """
    # Wilting point calculated following SUEWS
    # wp = wp_smd / g_smd

    g_smd_nom = 1 - np.exp(g_smd * (smd - wp_smd))
    g_smd_denom = 1 - np.exp(g_smd * (0 - wp_smd))
    g_smd = g_smd_nom / g_smd_denom
    return g_smd


def cal_gs_mod(kd, ta_c, rh, pa, smd, lai, g_cst, g_max, lai_max, wp_smd, debug=False):
    """Model surface conductance/resistance using phenology and atmospheric forcing conditions.

    Parameters
    ----------
    kd : numeric
        Incoming solar radiation [W m-2]
    ta_c : numeric
        Air temperature [degC]
    rh : numeric
        Relative humidity [%]
    pa : numeric
        Air pressure [Pa]
    smd : numeric
        Soil moisture deficit [mm]
    lai : numeric
        Leaf area index [m2 m-2]
    g_cst : size-6 array
        Parameters to determine surface conductance/resistance:
        g_lai (LAI related), g2 (solar radiation related),
        g_q1 (humidity related), g_q2 (humidity related),
        g_ta (air temperature related),
        g_smd (soil moisture related)
    g_max : numeric
        Maximum surface conductance [mm s-1]
    lai_max : numeric
        Maximum LAI [m2 m-2]
    wp_smd : numeric
        Wilting point indicated as soil moisture deficit [mm]

    Returns
    -------
    numeric
        Modelled surface conductance [mm s-1]
    """
    from atmosp import calculate as ac

    # broadcast g_lai – g_smd
    # print('g_cst', g_cst)
    g_lai, g2, g_q1, g_q2, g_ta, g_smd = g_cst
    # print(g_lai, g2, g_q1, g_q2, g_ta, g_smd)
    # lai related
    g_lai = cal_g_lai(lai, g_lai, lai_max)
    # print('g_lai', g_lai)

    # kdown related
    g_kd = cal_g_kd(kd, g2)
    # print('g_kd', g_kd)
    # dq related
    ta_k = ta_c + 273.15
    dq = cal_dq(rh, ta_c, pa / 100) * 1000 + rh * 0
    g_dq = cal_g_dq(dq, g_q1, g_q2)
    # print('g_dq', g_dq)
    # ta related
    ta_c = ta_k - 273.15
    g_ta = cal_g_ta(ta_c, g_ta)
    # print('g_ta', g_ta)
    # smd related
    g_smd = cal_g_smd(smd, g_smd, wp_smd)
    # print('g_smd', g_smd)
    # combine all corrections
    gs_c = g_lai * g_kd * g_dq * g_ta * g_smd
    gs = g_max * gs_c

    if debug:

        # pack results into a dataframe
        df_gs = pd.concat(
            {
                "gs_mod": pd.Series(gs),
                "gs_c": pd.Series(gs_c),
                "g_lai": pd.Series(g_lai),
                "g_kd": pd.Series(g_kd),
                "g_dq": pd.Series(g_dq),
                "g_ta": pd.Series(g_ta),
                "g_smd": pd.Series(g_smd),
            },
            axis=1,
        )

        return df_gs
    else:
        return gs


def calib_g(
    df_fc_suews,
    ser_ra,
    g_max,
    lai_max,
    wp_smd,
    method="cobyla",
    prms_init=None,
    debug=False,
):
    """Calibrate parameters for modelling surface conductance over vegetated surfaces using `LMFIT <https://lmfit.github.io/lmfit-py/model.html>`.

    Parameters
    ----------
    df_fc_suews : pandas.DataFrame
        DataFrame in `SuPy forcing <https://supy.readthedocs.io/en/latest/data-structure/df_forcing.html>`_ format
    ser_ra: pandas.Series
        Series with RA, aerodynamic resistance, [s m-1]
    g_max : numeric
        Maximum surface conductance [mm s-1]
    lai_max : numeric
        Maximum LAI [m2 m-2]
    wp_smd : numeric
        Wilting point indicated as soil moisture deficit [mm]
    method: str, optional
        Method used in minimisation by `lmfit.minimize`: details refer to its `method<lmfit:minimize>`.
    prms_init: lmfit.Parameters, optional
        Initial parameters for calibration
    debug : bool, optional
        Option to output final calibrated `ModelResult <lmfit:ModelResult>`, by default False

    Returns
    -------
    dict, or `ModelResult <lmfit:ModelResult>` if `debug==True`
        1. dict: {parameter_name -> best_fit_value}
        2. `ModelResult`

        Note:
            Parameters for surface conductance:
            g_lai (LAI related), g2 (solar radiation related),
            g_q1 (humidity related), g_q2 (humidity related),
            g_ta (air temperature related),
            g_smd (soil moisture related)

    Note
    ----
    For calibration validity, turbulent fluxes, QH and QE, in `df_fc_suews` should ONLY be observations, i.e., interpolated values should be avoided.
    To do so, please place `np.nan` as missing values for QH and QE.

    """
    from lmfit import Model, Parameters, Parameter

    list_var_sel = ["qh", "qe", "Tair", "RH", "pres", "kdown", "xsmd", "lai"]
    df_obs = df_fc_suews[list_var_sel].copy().dropna()

    # convert to Pa
    df_obs.pres *= 100

    gs_obs = cal_gs_obs(
        df_obs.qh, df_obs.qe, df_obs.Tair, df_obs.RH, df_obs.pres, ser_ra
    )

    def func_fit_g(kd, ta, rh, pa, smd, lai, g_lai, g_kd, g_q1, g_q2, g_ta, g_smd):
        gs = cal_gs_mod(
            kd,
            ta,
            rh,
            pa,
            smd,
            lai,
            [g_lai, g_kd, g_q1, g_q2, g_ta, g_smd],
            g_max,
            lai_max,
            wp_smd,
        )
        return gs

    gmodel = Model(
        func_fit_g,
        independent_vars=["lai", "kd", "ta", "rh", "pa", "smd"],
        param_names=["g_lai", "g_kd", "g_q1", "g_q2", "g_ta", "g_smd"],
    )
    if prms_init is None:
        print("Preset parameters will be loaded!")
        print("Please use with caution.")
        prms = Parameters()
        dict_prms_init = {
            "lai": 3.5,
            "kd": 50,
            "q1": 0.1,
            "q2": 0.7,
            "ta": 25,
            "smd": 0.05,
        }
        list_g = (
            Parameter(f"g_{var}", val, True, 0, None, None, None)
            for var, val in dict_prms_init.items()
        )
        prms.add_many(*list_g)
        # set specific bounds:
        # g_lai: LAI related
        prms["g_lai"].set(min=0, max=10)

        prms["g_kd"].set(min=0, max=300)
        # g_q1, g_q2: specific humidity related
        prms["g_q1"].set(min=0, max=1)
        prms["g_q2"].set(min=0, max=1)
        # g_ta: within reasonable temperature ranges
        prms["g_ta"].set(min=-10, max=55)
        # g_smd: within sensitive ranges of SMD
        prms["g_smd"].set(min=0.02, max=0.1)
    else:
        print("User provided parameters are loaded!")
        prms = prms_init

    # pack into a DataFrame for filtering out nan
    df_fit = pd.concat([gs_obs.rename("gs_obs"), df_obs], axis=1).dropna()

    res_fit = gmodel.fit(
        df_fit.gs_obs,
        kd=df_fit.kdown,
        ta=df_fit.Tair,
        rh=df_fit.RH,
        pa=df_fit.pres,
        smd=df_fit.xsmd,
        lai=df_fit.lai,
        params=prms,
        # useful ones: ['nelder', 'powell', 'cg', 'cobyla', 'bfgs', 'trust-tnc']
        method=method,
        #     nan_policy='omit',
        verbose=True,
    )

    # provide full fitted model if debug == True otherwise only a dict with best fit parameters
    res = res_fit if debug else res_fit.best_values

    return res


def fit_g_ta(ser_ta, ser_gs):
    from lmfit import Model, Parameters, Parameter

    model_g_ta = Model(
        cal_g_ta,
        independent_vars=["ta_c"],
        param_names=["g_ta", "tl", "th"],
    )
    prms = Parameters()
    prm_g_ta = Parameter(
        "g_ta", ser_ta.median(), vary=True, min=ser_ta.min(), max=ser_ta.max()
    )
    prm_tl = Parameter("tl", min=ser_ta.min(), vary=False)
    prm_th = Parameter("th", ser_ta.max(), vary=False)
    prms.add_many(prm_g_ta, prm_tl, prm_th)

    res_fit = model_g_ta.fit(ser_gs / ser_gs.max(), ta_c=ser_ta, params=prms)
    return res_fit


def fit_g_smd(ser_smd, ser_gs, wp_smd):
    from lmfit import Model, Parameters, Parameter

    model_g_smd = Model(
        cal_g_smd,
        independent_vars=["smd"],
        param_names=["g_smd", "wp_smd"],
    )
    prms = Parameters()
    prm_g_smd = Parameter("g_smd", 0.02, vary=True, min=1e-4, max=0.5)
    prm_wp_smd = Parameter("wp_smd", wp_smd, vary=False)
    prms.add_many(prm_g_smd, prm_wp_smd)

    res_fit = model_g_smd.fit(ser_gs / ser_gs.max(), smd=ser_smd, params=prms)
    return res_fit


def fit_g_kd(ser_kd, ser_gs):
    from lmfit import Model, Parameters, Parameter

    model_g_kd = Model(
        cal_g_kd,
        independent_vars=["kd"],
        param_names=["g_kd"],
    )
    prms = Parameters()
    prm_g_kd = Parameter("g_kd", 100, vary=True, min=10, max=300)
    prms.add(prm_g_kd)
    res_fit = model_g_kd.fit(ser_gs / ser_gs.max(), kd=ser_kd, params=prms)
    return res_fit


def fit_g_dq(ser_dq, ser_gs):
    from lmfit import Model, Parameters, Parameter


    model_g_dq = Model(
        cal_g_dq,
        independent_vars=["dq"],
        param_names=["g_q1", "g_q2"],
    )
    prms = Parameters()
    prm_g_q1 = Parameter(
        "g_q1",
        ser_gs.min() / ser_gs.max(),
        vary=False,
    )
    prm_g_q2 = Parameter("g_q2", 0.1, vary=True, min=0.01, max=0.95)
    prms.add_many(prm_g_q1, prm_g_q2)
    res_fit = model_g_dq.fit(ser_gs / ser_gs.max(), dq=ser_dq, params=prms)
    return res_fit


# calculate specific humidity using relative humidity
def cal_qa(rh_pct, theta_K, pres_hPa):
    from atmosp import calculate as ac

    qa = ac("qv", RH=rh_pct, p=pres_hPa * 100, theta=theta_K)
    return qa


# calculate specific humidity deficit [kg kg-1] using relative humidity
def cal_dq(rh_pct, ta_c, pres_hPa):
    from atmosp import calculate as ac

    ta_k = ta_c + 273.16
    pa = pres_hPa * 100
    dq = ac("qvs", T=ta_k, p=pa) - ac("qv", T=ta_k, p=pa, RH=rh_pct)

    return dq


# calculate relative humidity using specific humidity
def cal_rh(qa_kgkg, theta_K, pres_hPa):
    from atmosp import calculate as ac

    RH = ac("RH", av=qa_kgkg, p=pres_hPa * 100, theta=theta_K)
    return RH


# calculate latent heat of vaporisation
def cal_lat_vap(qa_kgkg, theta_K, pres_hPa):
    from atmosp import calculate as ac

    # wel-bulb temperature
    tw = ac(
        "Tw", qv=qa_kgkg, p=pres_hPa, theta=theta_K, remove_assumptions=("constant Lv")
    )
    # latent heat [J kg-1]
    Lv = 2.501e6 - 2370.0 * (tw - 273.15)
    return Lv


# calculate specific heat capacity of air [J kg-1 K-1]
def cal_cp(qa_kgkg, ta_K, pres_hPa):
    from atmosp import calculate as ac

    temp_C = ta_K - 273.15

    rh_pct = ac("RH", qv=qa_kgkg, T=ta_K, p=pres_hPa * 100)

    # # Garratt equation a20(1992)
    # cpd = 1005.0 + ((temp_C + 23.16) ** 2) / 3364.0

    # # Beer(1990) for water vapour
    # cpm = (
    #     1859
    #     + 0.13 * rh_pct
    #     + (19.3 + 0.569 * rh_pct) * (temp_C / 100.0)
    #     + (10.0 + 0.5 * rh_pct) * (temp_C / 100.0) ** 2
    # )

    # # air density
    # rho = ac("rho", qv=qa_kgkg, T=ta_K, p=pres_hPa * 100)

    # # water vapour mixing ratio
    # rv = ac("rv", qv=qa_kgkg, T=ta_K, p=pres_hPa * 100)

    # # dry air density
    # rho_d = rv / (1 + rv) * rho

    # # water vapour density
    # rho_v = rho - rho_d

    # # heat capacity of air
    # cp = cpd * (rho_d / (rho_d + rho_v)) + cpm * (rho_v / (rho_d + rho_v))

    cpa = cal_cp_with_rh(temp_C, rh_pct, pres_hPa)

    return cpa


# calculate specific heat capacity of air [J kg-1 K-1]
def cal_cp_with_rh(temp_C, rh_pct, pres_hPa):
    from atmosp import calculate as ac

    # Garratt equation a20(1992)
    cpd = 1005.0 + ((temp_C + 23.16) ** 2) / 3364.0

    # Beer(1990) for water vapour
    cpm = (
        1859
        + 0.13 * rh_pct
        + (19.3 + 0.569 * rh_pct) * (temp_C / 100.0)
        + (10.0 + 0.5 * rh_pct) * (temp_C / 100.0) ** 2
    )

    # density of dry air
    rho_d = cal_dens_dry(rh_pct, temp_C, pres_hPa)

    # density of vapour
    rho_v = cal_dens_vap(rh_pct, temp_C, pres_hPa)

    # heat capacity of air
    cpa = cpd * (rho_d / (rho_d + rho_v)) + cpm * (rho_v / (rho_d + rho_v))

    return cpa


# stability correction for momentum
def cal_psi_mom(zoL):
    # limit for neutral condition
    lim_neutral = 1e-5

    zoL = np.where(np.abs(zoL) > 5, 5 * np.sign(zoL), zoL)

    # stable, zoL>0
    zoL_stab = np.where(zoL > lim_neutral, zoL, 0)
    psim_stab = (-6) * np.log(1 + zoL_stab)

    # unstable, zoL<0
    zoL_unstab = np.where(zoL < -lim_neutral, zoL, 0)
    psim_unstab = 0.6 * (2) * np.log((1 + (1 - 16 * zoL_unstab) ** 0.5) / 2)

    # populate values with respect to stability
    psim = np.where(zoL > lim_neutral, psim_stab, psim_unstab)
    psim = np.where(np.abs(zoL) <= lim_neutral, 0, psim)

    return psim


# stability correction for heat
def cal_psi_heat(zoL):
    # limit for neutral condition
    lim_neutral = 1e-5

    zoL = np.where(np.abs(zoL) > 5, 5 * np.sign(zoL), zoL)

    # stable, zoL>0
    zoL_stab = np.where(zoL > lim_neutral, zoL, 0)
    psih_stab = -4.5 * zoL_stab

    # unstable, zoL<0
    zoL_unstab = np.where(zoL < -lim_neutral, zoL, 0)
    psih_unstab = (2) * np.log((1 + (1 - 16 * zoL_unstab) ** 0.5) / 2)

    # populate values with respect to stability
    psih = np.where(zoL > lim_neutral, psih_stab, psih_unstab)
    psih = np.where(np.abs(zoL) <= lim_neutral, 0, psih)

    return psih


# saturation vapour pressure [hPa]
def cal_vap_sat(Temp_C, Press_hPa):
    ser_pres_kPa = pd.Series(Press_hPa / 10)
    ser_TaC = pd.Series(Temp_C)
    # if ser_TaC is close to zero degC, set to 0.001 degC
    ser_TaC = ser_TaC.where(ser_TaC.abs() > 0.001, 0.001)

    # ser_es_hPa = cal_vap_sat(0.001, Press_hPa)

    # 0.001000 <= ser_TaC < 50:
    ser_emb_pos = 6.1121 * np.exp(
        ((18.678 - ser_TaC / 234.5) * ser_TaC) / (ser_TaC + 257.14)
    )
    ser_f_pos = 1.00072 + ser_pres_kPa * (3.2e-6 + 5.9e-10 * ser_TaC ** 2)

    # -40 < ser_TaC <= -0.001000:
    ser_emb_neg = 6.1115 * np.exp(
        ((23.036 - ser_TaC / 333.7) * ser_TaC) / (ser_TaC + 279.82)
    )
    ser_f_neg = 1.00022 + ser_pres_kPa * (3.83e-6 + 6.4e-10 * ser_TaC ** 2)

    # combine both conditions
    ser_emb = ser_emb_pos.where(ser_TaC >= 0.001, ser_emb_neg)
    ser_f = ser_f_pos.where(ser_TaC >= 0.001, ser_f_neg)
    ser_es_hPa = ser_emb * ser_f
    return ser_es_hPa


# density of dry air [kg m-3]
def cal_dens_dry(RH_pct, Temp_C, Press_hPa):
    gas_ct_dry = 8.31451 / 0.028965  # dry_gas/molar
    es_hPa = cal_vap_sat(Temp_C, Press_hPa)
    Ea_hPa = RH_pct / 100 * es_hPa
    dens_dry = ((Press_hPa - Ea_hPa) * 100) / (gas_ct_dry * (273.16 + Temp_C))
    return dens_dry


# density of vapour [kg m-3]
def cal_dens_vap(RH_pct, Temp_C, Press_hPa):
    gas_ct_wv = 8.31451 / 0.0180153  # dry_gas/molar_wat_vap
    es_hPa = cal_vap_sat(Temp_C, Press_hPa)
    Ea_hPa = RH_pct / 100 * es_hPa
    vap_dens = Ea_hPa * 100 / ((Temp_C + 273.16) * gas_ct_wv)
    return vap_dens


#
# # specific heat capacity of air mass [J kg-1 K-1]
# def cal_cpa(Temp_C, RH_pct, Press_hPa):
#     # heat capacity of dry air depending on air temperature
#     cpd = 1005.0 + ((Temp_C + 23.16) ** 2) / 3364.0
#     # heat capacity of vapour
#     cpm = (
#         1859
#         + 0.13 * RH_pct
#         + (19.3 + 0.569 * RH_pct) * (Temp_C / 100.0)
#         + (10.0 + 0.5 * RH_pct) * (Temp_C / 100.0) ** 2
#     )
#
#     # density of dry air
#     rho_d = cal_dens_dry(RH_pct, Temp_C, Press_hPa)
#
#     # density of vapour
#     rho_v = cal_dens_vap(RH_pct, Temp_C, Press_hPa)
#
#     # specific heat
#     cpa = cpd * (rho_d / (rho_d + rho_v)) + cpm * (rho_v / (rho_d + rho_v))
#     return cpa


# air density [kg m-3]
def cal_dens_air(Press_hPa, Temp_C):
    # dry_gas/molar
    gas_ct_dry = 8.31451 / 0.028965

    # air density [kg m-3]
    dens_air = (Press_hPa * 100) / (gas_ct_dry * (Temp_C + 273.16))
    return dens_air


# Obukhov length
def cal_Lob(QH, UStar, Temp_C, RH_pct, Press_hPa, g=9.8, k=0.4):
    # gravity constant/(Temperature*Von Karman Constant)
    G_T_K = (g / (Temp_C + 273.16)) * k

    # air density [kg m-3]
    rho = cal_dens_air(Press_hPa, Temp_C)

    # specific heat capacity of air mass [J kg-1 K-1]
    cpa = cal_cp_with_rh(Temp_C, RH_pct, Press_hPa)

    # Kinematic sensible heat flux [K m s-1]
    H = QH / (rho * cpa)

    # temperature scale
    # uStar = np.max([0.01, UStar])
    TStar = -H / UStar

    # Obukhov length
    Lob = (UStar ** 2) / (G_T_K * TStar)

    return Lob


# aerodynamic resistance
def cal_ra_obs(zm, zd, z0m, z0v, ws, Lob):
    """
    calculate aerodynamic resistance

    Parameters
    ----------
    zm : numeric
        measurement height
    zd : numeric
        displacement height
    z0m : numeric
        roughness length for momumtum
    z0v : numeric
        roughness length for vapour
    ws : numeric
        wind speed
    Lob : numeric
        Obukhov length

    Returns
    -------
    numeric
        aerodynamic resistance
    """
    # von Karman constant
    k = 0.4

    # effective measurement height
    zmd = zm - zd

    # z/L, stability scale
    zoL = zmd / Lob

    ra = (
        1
        / (k ** 2 * ws)
        * (np.log(zmd / z0m) - (cal_psi_mom(zoL) - cal_psi_mom(z0m / Lob)))
        * (np.log(zmd / z0v) - (cal_psi_heat(zoL) - cal_psi_heat(z0v / Lob)))
    )
    return ra