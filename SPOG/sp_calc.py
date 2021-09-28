# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import sys
import time
# import ray

__author__ = "Stephan Stock @ ZAH, Landessternwarte Heidelberg"
__version__ = "0.92"
__license__ = "MIT"


def calc_prob(model_list, params):
    df_return = pd.DataFrame({})
    for m in range(len(model_list)):
        for n in range(len(model_list[m])):
            df = get_mean_track(model_list[m][n])
            df['IMF_weight'] = IMF_prior(df)
            df['posterior_weight'] = np.exp(((params['ABL_star']-df['ABL'])**2/params['ABL_star_err']**2)-((params['color_star']-df['color'])**2/params['color_star_err']**2)-(
                (params['met_star']-df['[FE/H]'])**2/params['met_star_err']**2))*df['met_weight']*df['IMF_weight']*df['evol_weight']
            frames = [df_return, df]
            df_return = pd.concat(frames)
    return df_return

# get the mean value of section between evolutionary track points in
# order to derive the evolutionary time along the track


def get_mean_track(df):
    """Short summary.

    Parameters
    ----------
    df : type pd.DataFrame
        Dataframe consisting of evolutionary track model points

    Returns
    -------
    type pd.DataFrame
        Dataframe consisting of evolutionary track sections based on mean values between evolutionary track model points

    """
    df = df.reset_index(drop=True)
    df_dummy = df.diff()[1:]
    evol_time = df_dummy['age']  # get evolutionary time of track section
    evol_time_reset = evol_time.reset_index(drop=True)
    df_dummy_reset = df_dummy.reset_index(drop=True)
    df_dummy2 = df[:len(df)-1]
    df_dummy3 = df_dummy_reset.add(df_dummy2, fill_value=0)
    df_final = df_dummy2.add(df_dummy3, fill_value=0)/2.
    df_final['evol_weight'] = evol_time_reset  # required for timescale prior
    return (df_final)


# Salpeter IMF as prior


def IMF_prior(df):
    return df['mass_ZAMS']**-2.35


def opt_bin(df, weight, maxbin):
    print('WARNING: Automated binning not yet implemented, setting bins to 20')
    return 20
