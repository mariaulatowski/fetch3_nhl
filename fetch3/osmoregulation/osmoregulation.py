"""
Osmotic regulation module
"""
from fetch3.model_config import ConfigParams
from fetch3.nhl_transpiration.NHL_functions import *


R_GAS = 8.31446261815324 # Universal gas constant [J mol^-1 K^-1].


#this equation should be multipled by the feddes_root_stress equation in the roots.py script
def calc_osmotic_stress_roots(salinity, iv, R, T, h1,h2):
    """
    Calculates osmotic root water uptake stress function. Fit from data in greenhouse experiment, using the salt limiting function from Perri et al 2017
    Parameters
    ----------
    salinity: salinity in the soil [ppt]
    iv: van't hoff coefficient
    R: Universal Gas constant
    T: Temperature
    h: np.ndarray
        osmotic potential of soil [Pa]
    h1:osmotic potential which reaches maximum [Pa]
    h2: osmotic potential where decrease from max starts [Pa]

    Returns
    -------
    osmotic_stress_roots: np.ndarray
        Output of osmotic root water stress reduction function [unitless]
    """
    # Calculate osmotic potential from salinity
    C = salinity * 1000 / 58.44
    osmotic_potential = -C * iv * R * T

    # Define piecewise linear equations for different osmotic potential ranges
    if osmotic_potential >= h1:  # Use parameters for Segment 1
        slope, intercept = -5.4854372e-08,9.3140000e-01
    elif osmotic_potential<=h1 and osmotic_potential>=h2:  # Use parameters for Segment 2
        slope, intercept = 7.18522819e-07, 1.89857314e+00
    else:  # Use parameters for Segment 3
        slope, intercept = 1.01223963e-07, 6.17014944e-01

    # Calculate the linear function value
    osmotic_stress_roots = slope * osmotic_potential + intercept
    osmotic_stress_roots=max(0,osmotic_stress_roots)
    # Ensure result is non-negative
    return  osmotic_stress_roots


#needs to be added to the matric potential in the soil
def calc_osmotic_potential(salinity,E, R, iv, T):
    """
    Calculates osmotic potential at the root zone based on solute concentration and the
    filtration efficiency of NaCl from water at the roots from Perri et al 2018.

    Parameters
    ----------
    alpha_o : float
        parameter for E, filtration efficiency.
    beta_o: float
        parameter for E, filtration efficiency.
    Salinity : float
        salinity [ppt]
    R : float, optional
        Gas constant [J mol^-1 K^-1].
    iv: float, van't hoff coefficient
        default is 2 for NaCl
    T : float, optional
        Temperature [K].

    Returns
    -------
    osmotic_potential : float
        Osmotic potential at the root zone  [Pa]
    """
    C = salinity * 1000 / 58.44  # Convert salinity from g/L to mol/m^3 nacl
    osmotic_potential = -E*C * R *iv* T #calculate stem osmostic potential in pascals [Pa]
    return osmotic_potential


#needs to be added to the water potential in the stem
def calc_osmotic_potential_stem(salinity,E, R, iv, T):
    """
    Calculates osmotic potential in the xylem based on solute concentration and the
    filtration efficiency of NaCl from water at the roots.

    Parameters
    ----------
    alpha_o : float
        parameter for E, filtration efficiency.
    beta_o: float
        parameter for E, filtration efficiency.
    Salinity : float
        salinity [ppt]
    R : float, optional
        Gas constant [J mol^-1 K^-1]. Default is 8.314.
    iv: float, van't hoff coefficient
        default is 2 for NaCl
    T : float, optional
        Temperature [K].

    Returns
    -------
    osmotic_potential_stem : float
        Osmotic potential in the xylem  [Pa]
    """
    C = salinity * 1000 / 58.44  # Convert salinity from g/L to mol/m^3 nacl
    osmotic_potential_stem = -(1-E)*C * R *iv* T #calculate stem osmotic potential in pascals [Pa]
    return osmotic_potential_stem


def calc_Vcmax25(Vcmax25,a_sal,kvc,salinity):
    """
    Calculates Vcmax25 based off salinity in the soil. Equation fit from data from Suarez et al 2006

    Parameters
    ----------
    Vcmax25 : float
        maximium carboxyliation rate- Vcmax at 25 C at 0 ppt salinity [umol/umol]
    a_sal: float
        shape parameter for vcmax25 reduction function based off salinity (1/ppt)
    kvc: float
        shape parameter for vcmax25 reduction function based off salinity
    Salinity : float
        salinity [ppt]
    Returns
    -------
    Vcmax25_psu : float
        maximum carboxylation rate at 25 degrees at current salinity  [umol/umol]
    """
    Vcmax25_psu = Vcmax25 * (1 - ((a_sal * salinity) / (1 + a_sal * salinity)) ** kvc)
    return Vcmax25_psu


def calc_wp50_params(osmotic_potential, L_op, k_op, x0_op):
    """
    Calculates wp_s50 and c3 shape parameter for the stomatal conductance vs stem water potential curve according to
    the salinity (ppt) in the soil (using soil osmotic potential). Calculates the turgor loss point based off salinity in the soil and
    assumes wp_s50 is 85% of the tlp in accordance with christofferson et al 2016

    Parameters
    ----------
    L_op : float
        Empirical shape parameter of turgor loss point equation
    k_op: float
        Empirical shape parameter of turgor loss point equation
    x0_op: float
        Empirical shape parameter of turgor loss point equation [MPa]
    c3_b:
        b parameter of the c3 power function
    osmotic_potential : float
        soil osmotic potential [Pa]
    Returns
    -------
    wp_s50 : float
        stem water potential at 50% stomatal closure at current salinity [Pa]
    C3 : float
        shape parameter for the stem water potential vs stomatal response curve function at current salinity
    """

    turgor_loss_point= -L_op/ (1 + np.exp(k_op * (osmotic_potential*1e-6 - x0_op)))# turgor loss point--christofferson et al 2016, Bartlett et al 2012
    wp_s50 = 0.85 * turgor_loss_point #water potential at 50% stomatal closure -christofferson et al 2016
    c3= -2.406*wp_s50*(-wp_s50)**(-1.25) #xu et al 2023, christofferson et al 2016
    wp_s50=wp_s50*1e6
    return wp_s50, c3

def calc_NHL_osmo(cfg: ConfigParams, met_data, LADnorm_df, timestep, salinity_data):
    """
    Calculate NHL transpiration
    #TODO make docstring
    """
    # unpack config parameters
    dz = cfg.model_options.dz
    h = cfg.parameters.Hspec
    Cd = cfg.parameters.Cd
    Vcmax25 = calc_Vcmax25(cfg.parameters.Vcmax25, cfg.parameters.a_sal, cfg.parameters.kvc, salinity_data.loc[timestep, 'Salinity'])
    m = cfg.parameters.m
    alpha_p = cfg.parameters.alpha_p
    LAIp_sp = cfg.parameters.LAI
    LAIc_sp = cfg.parameters.LAIc_sp
    latitude = cfg.model_options.latitude
    longitude = cfg.model_options.longitude
    zenith_method = cfg.model_options.zenith_method
    x = cfg.parameters.x
    Cf = cfg.parameters.Cf
    time_offset = cfg.model_options.time_offset

    # unpack met data
    U_top = met_data.WS_F[timestep]
    ustar = met_data.USTAR[timestep]
    PAR = met_data.PPFD_IN[timestep]
    Ca = met_data.CO2_F[timestep]
    Tair = met_data.TA_F[timestep]
    Press = met_data.PA_F[timestep]
    doy = met_data.Timestamp.iloc[timestep].dayofyear
    time_of_day = met_data.Timestamp[timestep].hour + met_data.Timestamp[timestep].minute / 60

    # Look up the correct LADnorm for the model tree
    if cfg.model_options.LAD_column_labels is not None:
        LADnorm = LADnorm_df[cfg.model_options.LAD_column_labels[cfg.species]]
    else:
        LADnorm = LADnorm_df[cfg.species]

    z_h_LADnorm = LADnorm_df.z_h

    VPD = met_data.VPD_F_kPa[timestep]

    # Set up vertical grid
    zmin = 0
    z = np.arange(zmin, h, dz)  # [m]

    # Distrubute leaves vertically, and assign leaf area to stem
    LAD = calc_LAI_vertical(LADnorm, z_h_LADnorm, LAIc_sp, dz, h)  # [m2leaf m-2crown m-1stem]

    # Calculate wind speed at each layer
    U, Km = solve_Uz(z, dz, Cd, LAD, U_top, h=h)

    # Adjust the diffusivity and velocity by Ustar
    # Eqn A.5 from Mirfenderesgi et al 2016

    U = U * ustar
    Km = Km * ustar

    # Calculate radiation at each layer
    P0, Qp, zenith_angle = calc_rad_attenuation(
        PAR,
        LAD,
        dz,
        Cf=Cf,
        x=x,
        lat=latitude,
        long=longitude,
        doy=doy,
        time_of_day=time_of_day,
        time_offset=time_offset,
        zenith_method=zenith_method,
    )

    # Solve conductances
    A, gs, Ci, Cs, gb, geff = solve_leaf_physiology(Tair, Qp, Ca, Vcmax25, alpha_p, VPD=VPD, m=m, uz=U)

    # Calculate the transpiration per m-1 [ kg H2O s-1 m-1_stem]
    NHL_trans_leaf = calc_transpiration_leaf(VPD, Tair, geff, Press)  # [kg H2O m-2leaf s-1]
    NHL_trans_sp_stem = NHL_trans_leaf * LAD

    # Add data to dataset
    ds = xr.Dataset(
        data_vars=dict(
            U=(["z"], U),
            Km=(["z"], Km),
            P0=(["z"], P0),
            Qp=(["z"], Qp),
            A=(["z"], A),
            gs=(["z"], gs),
            Ci=(["z"], Ci),
            Cs=(["z"], Cs),
            gb=(["z"], gb),
            geff=(["z"], geff),
            NHL_trans_leaf=(["z"], NHL_trans_leaf),
            NHL_trans_sp_stem=(["z"], NHL_trans_sp_stem),
        ),
        coords=dict(z=(["z"], z)),
        attrs=dict(description="Model output"),
    )

    # Add metadata to dataset
    ds.NHL_trans_leaf.attrs = dict(
        units="kg H2O m-2leaf s-1", description="NHL transpiration per unit leaf area"
    )
    ds.NHL_trans_sp_stem.attrs = dict(
        units="kg H2O s-1 m-1_stem", description="NHL transpiration per unit height of stem"
    )

    return ds, LAD, zenith_angle

def calc_NHL_timesteps_osmo(cfg: ConfigParams, met_data, LADnorm_df, salinity_data,  **kwargs):
    """
    Calls NHL for each timestep in the met data
    #TODO docstring

    Parameters
    ----------
    dz : _type_
        _description_
    h : _type_
        _description_
    Cd : _type_
        _description_
    met_data : _type_
        _description_
    Vcmax25 : _type_
        _description_
    alpha_p : _type_
        _description_
    total_LAI_spn : _type_
        _description_
    plot_area : _type_
        _description_
    total_crown_area_spn : _type_
        _description_
    mean_crown_area_spn : _type_
        _description_
    LAD_norm : _type_
        _description_
    z_h_LADnorm : _type_
        _description_
    lat : _type_
        _description_
    long : _type_
        _description_
    time_offset : int, optional
        _description_, by default -5

    Returns
    -------
    _type_
        _description_
    """

    h = cfg.parameters.Hspec
    dz = cfg.model_options.dz
    zmin = 0
    z = np.arange(zmin, h, dz)  # [m]

    NHL_tot_trans_sp_tree_all = np.empty((len(met_data)))
    zenith_angle_all = np.empty((len(met_data)))

    datasets = []
    for i in range(0, len(met_data)):
        ds, LAD, zenith_angle = calc_NHL_osmo(cfg, met_data, LADnorm_df, i, salinity_data)

        # ds, LAD, zenith_angle = calc_NHL(
        #     dz, h, Cd, met_data.WS_F.iloc[i], met_data.USTAR.iloc[i], met_data.PPFD_IN.iloc[i], met_data.CO2_F.iloc[i], Vcmax25, alpha_gs, alpha_p,
        #     total_LAI_spn, plot_area, total_crown_area_spn, mean_crown_area_spn, LAD_norm, z_h_LADnorm,
        #     met_data.RH.iloc[i], met_data.TA_F.iloc[i], met_data.PA_F.iloc[i], doy = met_data.Timestamp.iloc[i].dayofyear, lat = lat,
        #     long= long, time_offset = time_offset, time_of_day = met_data.Timestamp[i].hour + met_data.Timestamp[i].minute/60, **kwargs)

        zenith_angle_all[i] = zenith_angle
        datasets.append(ds)
    d2 = xr.concat(datasets, pd.Index(met_data.Timestamp, name="time"))
    return d2, LAD, zenith_angle_all

#####################################################################################################################
def calc_A_gs_wplimit(cfg: ConfigParams, met_data, LADnorm_df, timestep, salinity_data, H):
    """
    Calculate NHL transpiration
    #TODO make docstring
    """
    # unpack config parameters
    dz = cfg.model_options.dz
    h = cfg.parameters.Hspec
    Cd = cfg.parameters.Cd
    Vcmax25 = calc_Vcmax25(cfg.parameters.Vcmax25, cfg.parameters.a_sal, cfg.parameters.kvc, salinity_data.loc[timestep, 'Salinity'])
    m = cfg.parameters.m
    alpha_p = cfg.parameters.alpha_p
    LAIp_sp = cfg.parameters.LAI
    LAIc_sp = cfg.parameters.LAIc_sp
    latitude = cfg.model_options.latitude
    longitude = cfg.model_options.longitude
    zenith_method = cfg.model_options.zenith_method
    x = cfg.parameters.x
    Cf = cfg.parameters.Cf
    time_offset = cfg.model_options.time_offset

    #unpack osmoregulation parameters
    L_op = cfg.parameters.L_op 
    k_op = cfg.parameters.k_op
    x0_op =cfg.parameters.x0_op
    E_f = cfg.parameters.filt_eff
    R = R_GAS 
    iv = cfg.parameters.iv
    Salinity = salinity_data.loc[timestep, 'Salinity'] 
    H = H[timestep]
    # unpack met data
    U_top = met_data.WS_F[timestep]
    ustar = met_data.USTAR[timestep]
    PAR = met_data.PPFD_IN[timestep]
    Ca = met_data.CO2_F[timestep]
    Tair = met_data.TA_F[timestep]
    Press = met_data.PA_F[timestep]
    doy = met_data.Timestamp.iloc[timestep].dayofyear
    time_of_day = met_data.Timestamp[timestep].hour + met_data.Timestamp[timestep].minute / 60

    # Look up the correct LADnorm for the model tree
    if cfg.model_options.LAD_column_labels is not None:
        LADnorm = LADnorm_df[cfg.model_options.LAD_column_labels[cfg.species]]
    else:
        LADnorm = LADnorm_df[cfg.species]

    z_h_LADnorm = LADnorm_df.z_h

    VPD = met_data.VPD_F_kPa[timestep]

    # Set up vertical grid
    zmin = 0
    z = np.arange(zmin, h, dz)  # [m]

    # Distrubute leaves vertically, and assign leaf area to stem
    LAD = calc_LAI_vertical(LADnorm, z_h_LADnorm, LAIc_sp, dz, h)  # [m2leaf m-2crown m-1stem]

    # Calculate wind speed at each layer
    U, Km = solve_Uz(z, dz, Cd, LAD, U_top, h=h)

    # Adjust the diffusivity and velocity by Ustar
    # Eqn A.5 from Mirfenderesgi et al 2016

    U = U * ustar
    Km = Km * ustar

    # Calculate radiation at each layer
    P0, Qp, zenith_angle = calc_rad_attenuation(
        PAR,
        LAD,
        dz,
        Cf=Cf,
        x=x,
        lat=latitude,
        long=longitude,
        doy=doy,
        time_of_day=time_of_day,
        time_offset=time_offset,
        zenith_method=zenith_method,
    )

    # Solve conductances
    A, gs, Ci, Cs, gb, geff = solve_leaf_physiology_wplimit(Tair, Qp, Ca, Vcmax25, alpha_p, VPD, m, uz=U, H=H, L_op=L_op, k_op=k_op, x0_op=x0_op, E_f=E_f, Salinity=Salinity, iv=iv, R=R)


    # Add data to dataset
    ds = xr.Dataset(
        data_vars=dict(
            U=(["z"], U),
            Km=(["z"], Km),
            P0=(["z"], P0),
            Qp=(["z"], Qp),
            A=(["z"], A),
            gs=(["z"], gs),
            Ci=(["z"], Ci),
            Cs=(["z"], Cs),
            gb=(["z"], gb),
            geff=(["z"], geff),
        ),
        coords=dict(z=(["z"], z)),
        attrs=dict(description="Model output"),
    )

    return ds, LAD, zenith_angle


def solve_leaf_physiology_wplimit(Tair, Qp, Ca, Vcmax25, alpha_p, VPD, m, H, L_op, k_op, x0_op, E_f, Salinity, iv, R, **kwargs):
    """
    Calculates photosynthesis and stomatal conductance
    Uses Leuning model for stomatal conductance

    Parameters
    ----------
    Tair : float
        Air temperature [deg C]
    Qp : float
        absorbed photosynthetically active radiation at each level within the canopy [µmol m-2 s-1]
    Ca : float
        CO2 concentration [µmol/mol]
    Vcmax25 : float
        Farquhar model parameter
    alpha_p : float
        Farquhar model parameter
    VPD : float
        Vapor pressure deficit [kPa]
    **kwargs for calc_gb

    Returns
    -------
    A : float
        photosynthesis [µmol m-2 s-1]
    gs : float
        stomatal conductance [mol m-2 s-1]
    Ci : float
        intracellular CO2 concentration [µmol mol-1]
    Cs : float
        CO2 concentration at leaf surface [µmol mol-1]
    gb : float
        boundary layer conductance [mol m-2 s-1]
    geff : float
        effective leaf conductance [mol m-2 s-1]
    """
    # Parameters
    Kc25 = 300  # [µmol mol-1] Michaelis-Menten constant for CO2, at 25 deg C
    Ko25 = 300  # [mmol mol-1] Michaelis-Menten constant for O2, at 25 deg C
    e_m = 0.08  # [mol mol-1]
    o = 210  # [mmol mol-1]
    g0 = 0.01  # [mol m-2 s-1]

    # Adjust the Farquhar model parameters for temperature
    Vcmax = Vcmax25 * np.exp(0.088 * (Tair - 25)) / (1 + np.exp(0.29 * (Tair - 41)))
    Kc = Kc25 * np.exp(0.074 * (Tair - 25))
    Ko = Ko25 * np.exp(0.018 * (Tair - 25))

    # Calculate gamma_star and Rd
    Rd = 0.015 * Vcmax  # Dark respiration [µmol m-2 s-1]
    gamma_star = (3.69 + 0.188 * (Tair - 25) + 0.0036 * (Tair - 25) ** 2) * 10

    def calc_Ac(Vcmax, Ci, gamma_star, Kc, o, Ko, Rd):
        return Vcmax * (Ci - gamma_star) / (Ci + Kc * (1 + o / Ko)) - Rd

    def calc_Aj(alpha_p, e_m, Qp, Ci, gamma_star, Rd):
        return alpha_p * e_m * Qp * (Ci - gamma_star) / (Ci + 2 * gamma_star) - Rd

    # Solve for An, gs, and Ci
    Ci = 0.99 * Ca
    Cs = Ca  # CO2 concentration at the surface
    err = 10000
    count = 0
    while (err > 0.01) & (count < 200):
        Aj = calc_Aj(alpha_p, e_m, Qp, Ci, gamma_star, Rd)
        Ac = np.full(len(Aj), calc_Ac(Vcmax, Ci, gamma_star, Kc, o, Ko, Rd))

        A = np.minimum(Ac, Aj)

        # Calculate stomatal conductance
        gs = calc_gs_wplimit(g0, m, A, Cs, gamma_star, VPD, H, L_op, k_op, x0_op, E_f, iv, Salinity, R, D0=3)

        # Calculate leaf boundary layer resistance
        gb, rb = calc_gb(**kwargs)
        Cs = np.maximum(Ca - A * rb, np.full(len(A), 0.1 * Ca))

        # Update Ci while enforcing the limits
        Ci2 = Cs - A / gs
        Ci2 = np.clip(Ci2, 150, 600)  # Enforce Ci limits
        
        err = max(np.abs(Ci - Ci2))
        Ci = Ci2
        count += 1

    geff = calc_geff(gb, gs)

    A[0] = A[1]
    Ci[0] = Ci[1]
    Cs[0] = Cs[1]
    gs[0] = gs[1]
    gb[0] = gb[1]
    geff[0] = geff[1]

    return A, gs, Ci, Cs, gb, geff


def calc_gs_wplimit(g0, m, A, c_s, gamma_star, VPD, H, L_op, k_op, x0_op, E_f, iv, Salinity, R, D0=3):
    """
    Calculates gs according to Anderegg 2017, includes stem wp limit

    Parameters
    ----------
    g0 : float
        cuticular conductance [mol m-2 s-1], residual stomatal conductance at the
        light compensation point (empirically fitted parameter)
    m : float
        empirically fitted parameter [unitless]
    A : float
        net CO2 assimilation rate [µmol CO2 m-2 s-1]
    c_s : float
        atmospheric CO2 concentration [µmol mol-1]
    gamma_star : float
        CO2 compensation point [µmol mol-1]
    VPD : float
        VPD [kPa]
    D0 : float
        reference vapor pressure [kPa], by default assumed to be 3.0 kPa

    Returns
    -------
    gs : float
        stomatal conductance [mol H2O m-2 s-1]
    """
    osmotic_potential = calc_osmotic_potential(Salinity, E_f, R, iv, 293)
    wp_s50, c3 = calc_wp50_params(osmotic_potential, L_op, k_op, x0_op)
    wp_response = np.exp(-(((H*1e6) / wp_s50) ** c3))
    gs = g0 + (m * abs(A)*wp_response) / ((c_s - gamma_star) * (1 + VPD / D0))
    return gs

def calc_A_gs_wplimit_timesteps_osmo(cfg: ConfigParams, met_data, LADnorm_df, salinity_data, H,  **kwargs):
    """
    Calls NHL for each timestep in the met data
    #TODO docstring

    Parameters
    ----------
    dz : _type_
        _description_
    h : _type_
        _description_
    Cd : _type_
        _description_
    met_data : _type_
        _description_
    Vcmax25 : _type_
        _description_
    alpha_p : _type_
        _description_
    total_LAI_spn : _type_
        _description_
    plot_area : _type_
        _description_
    total_crown_area_spn : _type_
        _description_
    mean_crown_area_spn : _type_
        _description_
    LAD_norm : _type_
        _description_
    z_h_LADnorm : _type_
        _description_
    lat : _type_
        _description_
    long : _type_
        _description_
    time_offset : int, optional
        _description_, by default -5

    Returns
    -------
    _type_
        _description_
    """

    h = cfg.parameters.Hspec
    dz = cfg.model_options.dz
    zmin = 0
    z = np.arange(zmin, h, dz)  # [m]

    NHL_tot_trans_sp_tree_all = np.empty((len(met_data)))
    zenith_angle_all = np.empty((len(met_data)))

    datasets = []
    for i in range(0, len(met_data)):
        ds, LAD, zenith_angle = calc_A_gs_wplimit(cfg, met_data, LADnorm_df, i, salinity_data, H)
        # ds, LAD, zenith_angle = calc_NHL(
        #     dz, h, Cd, met_data.WS_F.iloc[i], met_data.USTAR.iloc[i], met_data.PPFD_IN.iloc[i], met_data.CO2_F.iloc[i], Vcmax25, alpha_gs, alpha_p,
        #     total_LAI_spn, plot_area, total_crown_area_spn, mean_crown_area_spn, LAD_norm, z_h_LADnorm,
        #     met_data.RH.iloc[i], met_data.TA_F.iloc[i], met_data.PA_F.iloc[i], doy = met_data.Timestamp.iloc[i].dayofyear, lat = lat,
        #     long= long, time_offset = time_offset, time_of_day = met_data.Timestamp[i].hour + met_data.Timestamp[i].minute/60, **kwargs)

        zenith_angle_all[i] = zenith_angle
        datasets.append(ds)
    d2 = xr.concat(datasets, pd.Index(met_data.Timestamp, name="time"))

    # Extract A and gs data from the concatenated dataset
    A_data = d2['A'].to_dataframe().reset_index()
    gs_data = d2['gs'].to_dataframe().reset_index()
    Ci_data =d2['Ci'].to_dataframe().reset_index()
    # Merge A and gs data with z coordinates
    A_data.rename(columns={'A': 'A'}, inplace=True)
    gs_data.rename(columns={'gs': 'gs'}, inplace=True)
    Ci_data.rename(columns={'Ci': 'Ci'}, inplace=True)
    # Combine A and gs data
    combined_data = A_data.merge(gs_data, on=['time', 'z'])
    combined_data = combined_data.merge(Ci_data, on = ['time', 'z'])

    return combined_data, d2, LAD, zenith_angle_all