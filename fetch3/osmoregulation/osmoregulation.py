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


def calc_Vcmax25(Vcmax25_m,Vcmax25_b,salinity):
    """
    Calculates Vcmax25 based off salinity in the soil. Equation fit from data from Suarez et al 2006

    Parameters
    ----------
    Vcmax25_m : float
        slope of the Vcmax equation
    Vcmax25_b: float
        intercept of the Vcmax equation
    Salinity : float
        salinity [ppt]
    Returns
    -------
    Vcmax25 : float
        maximum carboxylation rate at 25 degrees at current salinity  [umol/umol]
    """
    Vcmax25 = Vcmax25_m * salinity + Vcmax25_b
    return Vcmax25


def calc_wp50_params(E_m, osmotic_potential):
    """
    Calculates wp_s50 and c3 shape parameter for the stomatal conductance vs stem water potential curve according to
    the salinity (ppt) in the soil (using soil osmotic potential)

    Parameters
    ----------
    wp_s50m : float
        slope of the wp_s50 equation
    wp_s50b: float
        intercept of the wp_s50 equation
    c3_m: float
        m paramater of the c3 power function
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

    op_full_turgor=2.82*osmotic_potential/10e6-3.04 #osmotic potential at full turgor--fit from data from Nguyen et al 2017
    turgor_loss_point=(op_full_turgor*E_m)/(op_full_turgor + E_m) # turgor loss point--christofferson et al 2016, Bartlett et al 2012
    wp_s50 = 0.85 * turgor_loss_point #water potential at 50% stomatal closure -christofferson et al 2016
    c3= -2.406*wp_s50*(-wp_s50)**(-1.25) #xu et al 2023, christofferson et al 2016
    wp_s50=wp_s50*10e6
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
    Vcmax25 = calc_Vcmax25(cfg.parameters.Vcmax25_m, cfg.parameters.Vcmax25_b, salinity_data.loc[timestep, 'Salinity'])
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



