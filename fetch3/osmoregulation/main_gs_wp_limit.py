"""
This file runs the NHL_gs_wp_limit module which reruns the nhl functions for gs and A with the stem water potential limitor
to get actual A and gs modeled values.

It writes the output for A and gs as a csv file.

"""

import logging
import time
from pathlib import Path

import numpy as np
import pandas as pd

from fetch3.met_data import prepare_ameriflux_data
from fetch3.nhl_transpiration.NHL_functions import *
from fetch3.model_config import ConfigParams
from fetch3.osmoregulation.osmoregulation import calc_A_gs_wplimit_timesteps_osmo
from fetch3.osmoregulation.salinity_data import prepare_salinity_data



from fetch3.scaling import trans2d_to_tree

def main_gs_wp_limit(cfg: ConfigParams, output_dir, data_dir, H, to_model_res=True, write_output=False):
    """
    Calculate NHL transpiration with the dynamic vcmax that changes with salinity.

    Parameters
    ----------
    cfg : FETCH3 config object
        Model configuration
    output_dir : pathlib.Path
        Directory to write output files to. If output_dir is None, NHL output is not written.
    data_dir : pathlib.Path
        Directory containing input data
    to_model_res : bool, optional
        Whether or not to calculate and return NHL in the resolution needed to run in FETCH, by default True.

    Returns
    -------
    NHL transpiration : xarray.DataArray
        NHL transpiration. If to_model_res is True, NHL is returned in the resolution
        needed to run in FETCH, with units of [m3 H2O m-2crown m-1stem s-1]. If to_model_res is False,
        NHL is returned as tree-level transpiration, with units of [kg H2O s-1].
    LAD : array
        LAD profile
    DS: xarray.Dataset
    """
    logger = logging.getLogger(__name__)

    start = time.time()

    # Read in LAD and met data
    met_data = prepare_ameriflux_data(data_dir / cfg.model_options.input_fname, cfg)
    salinity_data=prepare_salinity_data(data_dir / cfg.model_options.input_fname, cfg)
    import pandas as pd

    LADnorm_df = pd.read_csv(data_dir / cfg.model_options.LAD_norm)

    logger.info("Calculating NHL...")

    combined_data, ds, LAD, zen = calc_A_gs_wplimit_timesteps_osmo(cfg, met_data, LADnorm_df, salinity_data, H)
    combined_data.to_csv(output_dir /'A_gs_wplimit_timesteps.csv', index=False)
    return ds
