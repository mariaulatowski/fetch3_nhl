"""
Root functions
"""
import numpy as np

from fetch3.model_config import ConfigParams


def feddes_root_stress(theta, theta1, theta2):
    """
    Calculates Feddes root water uptake stress function, from equations S.73, S.74, and S.75 of Silva et al. 2022

    Parameters
    ----------
    theta: np.ndarray
        Volumetric soil water content [m3 m-3]
    theta1: np.ndarray
        Soil water content below which root water uptake ceases [m3 m-3]
    theta2: np.ndarray
        Soil water content below which root water uptake starts decreasing [m3 m-3]

    Returns
    -------
    stress_roots: np.ndarray
        Output of Feddes root water stress reduction function [unitless]
    """

    conditions = [theta <= theta1, (theta > theta1) & (theta <= theta2), theta > theta2]
    outputs = [0, (theta - theta1) / (theta2 - theta1), 1]
    stress_roots = np.select(conditions, outputs)

    return stress_roots


def verma_root_mass_dist(cfg: ConfigParams):
    """
    Root mass distribution following Verma et al. 2014

    Parameters
    ----------
    cfg : dataclass
        Model configuration dataclass

    Returns
    -------
    r_dist: np.ndarray
        Root mass distribution

    """

    z_dist = np.arange(0, cfg.parameters.Root_depth + cfg.model_options.dz, cfg.model_options.dz)
    z_dist = np.flipud(z_dist)

    r_dist = (
        np.exp(cfg.parameters.qz - ((cfg.parameters.qz * z_dist) / cfg.parameters.Root_depth))
        * cfg.parameters.qz**2
        * (cfg.parameters.Root_depth - z_dist)
    ) / (cfg.parameters.Root_depth**2 * (1 + np.exp(cfg.parameters.qz) * (-1 + cfg.parameters.qz)))
    return r_dist


def calc_root_K(r_dist, stress_roots, cfg: ConfigParams):
    """

    Parameters
    ----------
    r_dist :
    stress_roots :
    cfg :

    Returns
    -------
    Effective root radial conductivity
    """

    # specific radial conductivity under saturated soil conditions
    Ksrad = stress_roots * cfg.parameters.Kr  # stress function is unitless

    # effective root radial conductivity
    Kerad = Ksrad * r_dist  # [1/sPa] #Kr is already divided by Rho*g

    return Kerad
