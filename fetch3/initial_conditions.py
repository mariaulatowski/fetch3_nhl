import numpy as np
from fetch3.model_config import ConfigParams

def calc_potential_vangenuchten(theta, theta_r, theta_s, alpha, m, n, rho, g, osmotic_potential=0):
    """
    Calculates water potential from soil moisture, including osmotic potential, using van Genuchten equation.

    Parameters
    ----------
    theta : float or np.ndarray
        soil water content [m3 m-3]
    theta_r : float
        residual water content [m3 m-3]
    theta_s : float
        saturated water content [m3 m-3]
    alpha : float
        empirical van Genuchten parameter [m-1]
    m : float
        empirical van Genuchten parameter [unitless]
    n : float
        empirical van Genuchten parameter [unitless]
    rho : float
        density of water [kg m-3]
    g : float
        gravitational constant [m s-2]
    osmotic_potential : float
        osmotic potential [Pa], usually a negative value that reduces the total water potential

    Returns
    -------
    water_potential_Pa : float or np.ndarray
        total water potential [Pa], combining matric and osmotic potential
    """
    # Calculate matric potential using van Genuchten equation
    effective_saturation = (theta - theta_r) / (theta_s - theta_r)
    water_potential_m = -((((1 / effective_saturation) ** (1 / m) - 1) ** (1 / n)) / alpha)
    water_potential_Pa = water_potential_m * rho * g

    # Add osmotic potential (if any) to matric potential
    total_water_potential_Pa = water_potential_Pa + osmotic_potential

    return total_water_potential_Pa  # [Pa]

def initial_conditions(cfg: ConfigParams, q_rain, zind):
    """
    Calculate initial water potential conditions, including osmotic potential.

    Parameters
    ----------
    cfg : ConfigParams
        model configuration
    q_rain : np.ndarray
        array of rain data
    zind : dataclass
        z index dataclass

    Returns
    -------
    H_initial: np.ndarray
        initial values for water potential [Pa] over the concatenated z domain (soil, roots, xylem)
    Head_bottom_H: np.ndarray
        water potential [Pa] for the bottom boundary. size is len(number of timesteps)
    """
    # Example osmotic potential values [Pa]
    osmotic_potential_clay = -2
    osmotic_potential_sand = -2
    osmotic_potential_roots = -1.5
    osmotic_potential_xylem = -1.2

    # Soil initial conditions with osmotic potential
    H_initial_soil = np.piecewise(
        zind.z_soil,
        [zind.z_soil <= cfg.parameters.clay_d, zind.z_soil > cfg.parameters.clay_d],
        [
            calc_potential_vangenuchten(
                cfg.parameters.initial_swc_clay,
                cfg.parameters.theta_R1,
                cfg.parameters.theta_S1,
                cfg.parameters.alpha_1,
                cfg.parameters.m_1,
                cfg.parameters.n_1,
                cfg.Rho,
                cfg.g,
                osmotic_potential=osmotic_potential_clay
            ),
            calc_potential_vangenuchten(
                cfg.parameters.initial_swc_sand,
                cfg.parameters.theta_R2,
                cfg.parameters.theta_S2,
                cfg.parameters.alpha_2,
                cfg.parameters.m_2,
                cfg.parameters.n_2,
                cfg.Rho,
                cfg.g,
                osmotic_potential=osmotic_potential_sand
            ),
        ],
    )

    # Root layer initial conditions with osmotic potential
    z_root_start = np.round(cfg.parameters.Soil_depth - cfg.parameters.Root_depth, decimals=5)
    H_initial_root_bottom = H_initial_soil[zind.z_soil == z_root_start] + osmotic_potential_roots
    H_initial_root = H_initial_root_bottom - (zind.z_root - z_root_start) * cfg.Rho * cfg.g

    # Xylem initial conditions with osmotic potential
    H_initial_xylem_bottom = H_initial_root_bottom + osmotic_potential_xylem
    H_initial_xylem = H_initial_xylem_bottom - (zind.z_upper - z_root_start) * cfg.Rho * cfg.g

    # Concatenate initial conditions for full domain
    H_initial = np.concatenate((H_initial_soil, H_initial_root, H_initial_xylem))

    # Bottom boundary condition with osmotic potential for soil
    Head_bottom_H = np.full(
        len(q_rain),
        calc_potential_vangenuchten(
            cfg.parameters.soil_moisture_bottom_boundary,
            cfg.parameters.theta_R1,
            cfg.parameters.theta_S1,
            cfg.parameters.alpha_1,
            cfg.parameters.m_1,
            cfg.parameters.n_1,
            cfg.Rho,
            cfg.g,
            osmotic_potential=osmotic_potential_clay
        ),
    )

    # Set bottom boundary for initial condition
    if cfg.model_options.BottomBC == 0:
        H_initial[0] = Head_bottom_H[0]
    
    return H_initial, Head_bottom_H
