"""
Extract salinity data from the met input file
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from fetch3.model_config import ConfigParams
from fetch3.utils import interpolate_2d

def import_salinity_data(filein):
    """
    Imports salinity data to a DataFrame.
    Parameters
    ----------
    filein : str
        Filepath to data

    Returns
    -------
    df: pd.DataFrame
        salinity data
    """

    df = pd.read_csv(filein, parse_dates=["TIMESTAMP_START", "TIMESTAMP_END"], na_values=["-9999"])

    return df


def prepare_salinity_data(filein, cfg: ConfigParams):

    df = import_salinity_data(filein)
    df = df.rename(columns={"TIMESTAMP_START": "Timestamp"})

    # Rename
    if cfg.model_options.met_column_labels is not None:
        df = df.rename(columns=cfg.model_options.met_column_labels)


    varlist = [
        "Timestamp",
        "Salinity"
        ]
 

    # Keep only variables needed
    df = df[varlist]

    # fill any gaps in the data
    df = df.set_index('Timestamp')
    df = df.interpolate(method='linear', limit_direction="both", limit_area=None)
    df = df.reset_index()

    # Select data for length of run
    df = df[(df.Timestamp >= cfg.model_options.start_time) & (df.Timestamp <= cfg.model_options.end_time)].reset_index(
        drop=True
    )

    return df

# Helper functions
def calc_model_time_grid(df, cfg: ConfigParams):
    tmax = len(df) * cfg.model_options.dt
    t_data = np.arange(cfg.model_options.tmin, tmax, cfg.model_options.dt)  # data time grids for input data
    t_data = list(t_data)
    nt_data = len(t_data)  # length of input data
    return tmax, t_data, nt_data


def interp_to_model_res(var, tmax, t_data, dt0):
    return np.interp(np.arange(0, tmax + dt0, dt0), t_data, var)


def prepare_osmo_data(cfg: ConfigParams, data_dir, z_upper):
    ###########################################################
    # Load and format input data
    ###########################################################

    # Input file
    data_path = data_dir / cfg.model_options.input_fname

    start_time = pd.to_datetime(cfg.model_options.start_time)
    end_time = pd.to_datetime(cfg.model_options.end_time)

    # read input data
    df = prepare_salinity_data(data_path, cfg)

    df = df.set_index("Timestamp")

    tmax, t_data, nt_data = calc_model_time_grid(df, cfg)

    Salinity=df['Salinity']
    Salinity=Salinity.interpolate(method='linear')

    Salinity = interp_to_model_res(Salinity, tmax, t_data, cfg.model_options.dt0)

    ####2d interpolation of met data
    Salinity_2d = interpolate_2d(import_salinity_data, len(z_upper))

    return Salinity_2d, tmax, start_time, end_time


@dataclass
class Salinity:
    """
    Dataclass to hold salinity data
    """

    Salinity_2d: np.ndarray
