"""
Osmotic regulation module
"""

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
def calc_osmotic_potential(salinity,alpha_o,beta_o, R, iv, T):
    """
    Calculates osmotic potential at the root zone based on solute concentration and the
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
    E= alpha_o+beta_o*C #filtration efficiency of salt at the root zone
    osmotic_potential = -E*C * R *iv* T #calculate stem osmostic potential in pascals [Pa]
    return osmotic_potential


#needs to be added to the water potential in the stem
def calc_osmotic_potential_stem(salinity,alpha_o,beta_o, R, iv, T):
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
    E= alpha_o+beta_o*C #filtration efficiency of salt at the root zone
    osmotic_potential_stem = -(1-E)*C * R *iv* T #calculate stem osmostic potential in pascals [Pa]
    return osmotic_potential_stem


def calc_Vcmax25(Vc_m,Vc_b,salinity):
    """
    Calculates Vcmax25 based off salinity in the soil. Equation fit from data from Suarez et al 2006

    Parameters
    ----------
    Vc_m : float
        slope of the Vcmax equation
    Vc_b: float
        intercept of the Vcmax equation
    Salinity : float
        salinity [ppt]
    Returns
    -------
    Vcmax25 : float
        maximum carboxylation rate at 25 degrees at current salinity  [umol/umol]
    """
    Vcmax25 = Vc_m * salinity + Vc_b
    return Vcmax25


def calc_wp50_params(wp_s50m, wp_s50b, c3_m, c3_b, salinity):
    """
    Calculates wp_s50 and c3 shape parameter for the stomatal conductance vs stem water potential curve according to
    the salinity (ppt) in the soil

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
    Salinity : float
        salinity [ppt]
    Returns
    -------
    wp_s50 : float
        stem water potential at 50% stomatal closure at current salinity [Pa]
    C3 : float
        shape parameter for the stem water potential vs stomatal response curve function at current salinity
    """
    wp_s50 = wp_s50m * salinity + wp_s50b
    c3= c3_m*salinity**c3_b
    return wp_s50, c3


#add this to the water potential in the xylem when passed to the capacitance equation in porous_media_xylem in the model_functions script
def calc_osmoregulation(arg,m_o,b_o):
    """
    Calculates osmotic potential in the storage based off the relationship between sap salinity and stem water potential from
    LÓPEZ-PORTILLO et al 2005

    Parameters
    ----------
    arg : float
        stem water potential [Pa]
    m_o: float
        shape parameter for osmoregulation function.
    b_o : float
        shape parameter for osmoregulation function.

    Returns
    -------
    os : float
        omsoregulation in the xylem. Osmotic potential in xylem storage based off sap salinity concentration [Pa]
    """
    os = m_o*arg-b_o
    return os
