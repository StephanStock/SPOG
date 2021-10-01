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
    """This function creates a new dataframe that is based on the mean values between evolutionary track model points.

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe consisting of evolutionary track model points

    Returns
    -------
    pd.DataFrame
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


def IMF_prior(df):
    """Calculate the weight of a certain model point based on the initial mass function (IMF) and zero-age main sequence mass.

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe consisting of evolutionary track model points

    Returns
    -------
    pd.DataFrame
         Dataframe with additional weight column based on IMF
    """
    return df['mass_ZAMS']**-2.35


def opt_bin(df, weight, maxbin):
    """Calculate the optimal number of bins, for a histogram consisting of weighted discrete measurements or model points
       using the algorithm provided in Hogg (2008) (arkiv:0807.4820v1) .

    Parameters
    ----------
    df : pd.DataFrame
        Includes the parameter for which the optimal binning is derived
    weight : pd.DataFrame
        The posterior weights of each posterior sample.
    maxbin : float
        sets the maximum number of possible bins

    Returns
    -------
    int
        An integer providing the optimal binning according to the algorithm

    """
    print('Calculating optimal binning for '+str(df.name))
    df_comb = pd.DataFrame({'Parameter': df,
                           'weight': weight}).sort_values(by=['Parameter'], ascending=True).reset_index(drop=True)

    par_range = abs(df_comb['Parameter'].max()-df_comb['Parameter'].min())
    # fix smoothing parameter alpha to be on the order of the mean weight, see Hogg(2008)
    alpha = df_comb['weight'].mean()
    df_bins = pd.DataFrame(0., index=np.arange(maxbin-2),
                           columns=['Number_bins', 'Likelihood'])
    for n in range(2, int(maxbin)):
        df_bins['Number_bins'][n-2] = n
        D_i = par_range/n
        df_bins_weight = pd.DataFrame(0., index=np.arange(n), columns=['Bin_number', 'Weight'])
        low_bound = 0
        for k in range(n):
            df_bins_weight['Bin_number'][k] = k
            index = np.max(np.where(df_comb['Parameter'] <= (df_comb['Parameter'].min()+(k+1)*D_i)))
            df_bins_weight['Weight'][k] = df_comb['weight'][low_bound:index+1].sum()
            low_bound = index+1
        for j in range(len(df_comb['weight'])):
            bin = get_bin(df_comb['Parameter'][j], df_comb['Parameter'].min(), D_i)
            if bin < 0:
                bin = 0
            if bin >= df_bins['Number_bins'][n-2]:
                bin = df_bins['Number_bins'][n-2]-1
            df_bins['Likelihood'][n-2] = df_bins['Likelihood'][n-2] + df_comb['weight'][j]*np.log(
                (df_bins_weight['Weight'][bin]+alpha-df_comb['weight'][j]) / (D_i*(((df_bins_weight['Weight']+alpha).sum())-df_comb['weight'][j])))

    print('Bins set to '+str(df_bins['Number_bins'][df_bins['Likelihood'].idxmax()]))
    return int(df_bins['Number_bins'][df_bins['Likelihood'].idxmax()])


def get_bin(x, xmin, D_i):
    """Get the current bin of a certain value.

    Parameters
    ----------
    x : float
        A Parameter value of the posterior.
    xmin : float
        The minum paramter value of the posterior.
    D_i : float
        Range of a bin.

    Returns
    -------
    float
        The bin in which x is placed.

    """
    return np.floor(((x-xmin)/D_i)-1e-10)
