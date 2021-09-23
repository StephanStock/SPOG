# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import os
import sys
import yaml
import h5py
import time
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import ThreadPoolExecutor

__author__ = "Stephan Stock @ ZAH, Landessternwarte Heidelberg"
__version__ = "0.9"
__license__ = "MIT"

import st_calc
import st_utils


print(f' \n'
      f'Welcome to SPOG V{__version__}\n')
# Input the parameter set in the parameter file - the list of possible
# arguments is defined in the default_params-dictionairy
with open(str(sys.argv[1])) as stream:
    try:
        input_params = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)
default_params = {'object_name': '',
                  'save_path': './',
                  'model_path': [],
                  'mag_star': 9,
                  'mag_star_err': 0.1,
                  'color_star': 0.97,
                  'color_star_err': 0.036,
                  'met_star': 0,
                  'met_star_err': 0.1,
                  'par_star': 2-45,
                  'par_err_star': 0.01,
                  'par_err_dom': False,
                  'use_extinction': False,
                  'A_lambda': 0,
                  'E_color': 0,
                  'parameterization': 'default',
                  'model_sampling': 1}
params = {}
for key in default_params:
    try:
        params[key] = input_params[key]
    except KeyError:
        params[key] = default_params[key]
        print(f'Parameter {key} is missing in the input file.\n'
              f'Value set to default: {key} = {default_params[key]}\n')

# calculate ABL from parallax and magnitude (This parametrization is key to be less biased!)
if params['use_extinction'] == False:
    params['ABL_star'] = 10**(params['mag_star']*0.2+1)*(params['par_star']/1000.)
    # Assumes parallax error dominates photometric uncertainties!
    if params['par_err_dom'] == True:
        params['ABL_star_err'] = 10**(params['mag_star']*0.2+1)*(params['par_err_star']/1000.)
    # Assumes magnitude uncertainty is not negligible! Slight bias due to non-linear transformation between magnitude and ABL
    else:
        params['ABL_star_err'] = 10**(params['mag_star']*0.2+1)*(params['par_err_star']/1000.) + \
            params['mag_star_err']*4.460517*(params['par_star']/1000) * \
            np.exp(0.460517*params['mag_star'])
else:
    print('Applying extinction...')
    params['ABL_star'] = 10**((params['mag_star']-params['A_lambda'])
                              * 0.2+1)*(params['par_star']/1000.)
    params['ABL_star_err'] = 10**(params['mag_star']*0.2+1)*(params['par_err_star']/1000.)
    #params['color_star_err'] = np.sqrt(params['color_star']**2+params['E_color']**2)
    params['color_star'] = params['color_star']-params['E_color']

# load hdf5 file for metallicity list
with h5py.File(params['model_path'], 'r') as hdf:
    metallicities = list(hdf.keys())

# check if parameters within correct bounds
if isinstance(params['model_sampling'], int) == False or params['model_sampling'] < 1:
    params['model_sampling'] = 1
    print('Parameter model_samling wrong format or value (must be integer greater than 0), setting to default value of 1 ')

# apply sparser sampling if requested
metallicities = metallicities[::params['model_sampling']]

Z = []

for i in range(np.size(metallicities)):
    dummy = metallicities[i]
    dummy2 = dummy.split('Z')
    dummy3 = dummy2[1].split('Y')
    Z.append(float(dummy3[0]))

Z = np.array(Z)


# Calculates the models metallicities given in atomic fraction Z to [FE/H]
# using the relation Y=0.2485+1.78Z by the Padova group

Zsun = 0.0152
Ysun = 0.2485+1.78*Zsun
Xsun = 1.-Ysun-Zsun
Y = 0.2485+1.78*Z
X = 1.-Y-Z
Fe_mod = np.log10(Xsun/Zsun)+np.log10(Z/X)


# make a pandas data frame and sort ascending by Fe/H
df_all_met = pd.DataFrame({'folder': metallicities, 'Z': Z, 'Fe/H': Fe_mod})
df_all_met_sort = df_all_met.sort_values(by=['Fe/H'], ascending=True)
df_only_fe = pd.DataFrame({'Fe/H': Fe_mod})
df_only_fe_sort = df_only_fe.sort_values(by=['Fe/H'], ascending=True)


# Important: create metallicity weights of available models for the correct
# probability estimation
df_only_fe_diff = df_only_fe_sort.diff()[1:]


met_weights = []

for n in range(np.size(metallicities)-1):
    if n == 0:
        met_weights.extend([df_only_fe_diff.iloc[n]['Fe/H']])
    elif n == (np.size(metallicities)-2):
        met_weights.extend([df_only_fe_diff.iloc[n]['Fe/H']/2. +
                           df_only_fe_diff.iloc[n-1]['Fe/H']/2.])
        met_weights.extend([df_only_fe_diff.iloc[n]['Fe/H']])
    else:
        met_weights.extend([df_only_fe_diff.iloc[n]['Fe/H']/2. +
                           df_only_fe_diff.iloc[n-1]['Fe/H']/2.])

# update DataFrame
df_all_met_sort['weight'] = np.asarray(met_weights)

# Now load models for parameter estimation of the star (only 5 sigma range of measured metallicity!)


maxfeh = params['met_star']+5*params['met_star_err']
minfeh = params['met_star']-5*params['met_star_err']

metallicities_to_load = df_all_met_sort.loc[(
    df_all_met_sort['Fe/H'] >= minfeh) & (df_all_met_sort['Fe/H'] <= maxfeh)]['folder']
weights_of_metallicities = df_all_met_sort.loc[(
    df_all_met_sort['Fe/H'] >= minfeh) & (df_all_met_sort['Fe/H'] <= maxfeh)]['weight']

weights_of_metallicities = weights_of_metallicities/len(weights_of_metallicities)


if __name__ == '__main__':
    metlist_rgb = []
    metlist_hb = []
    t0 = time.time()
    print('loading...')
    for model, weight in zip(metallicities_to_load, weights_of_metallicities):

        with h5py.File(params['model_path'], 'r') as hdf:
            mass_group_lowmass = hdf.get(model+'/lowmass')
            mass_group_lowmass_items = list(mass_group_lowmass.items())

            mass_group_highmass = hdf.get(model+'/highmass')
            mass_group_highmass_items = list(mass_group_highmass.items())

            mass_group_hb = hdf.get(model+'/hb')
            mass_group_hb_items = list(mass_group_hb.items())

            metlist_rgb.append(st_utils.load_models(hdf, mass_group_lowmass_items,
                                                    model, weight, params, 6., 10.999, 'lowmass'))
            metlist_rgb.append(st_utils.load_models(hdf, mass_group_highmass_items,
                                                    model, weight, params, 6., 10.999, 'highmass'))
            metlist_hb.append(st_utils.load_models(hdf, mass_group_hb_items,
                                                   model, weight, params, 1., 5., 'hb'))
            metlist_hb.append(st_utils.load_models(hdf, mass_group_highmass_items,
                                                   model, weight, params, 11., 14.999, 'highmass'))

        print('loaded models of metallicity: '+model)

    t1 = time.time()
    total = t1-t0
    print('Time to load models was '+str(total)+'s')

    print('All models loaded, starting calculations...')

    with ProcessPoolExecutor(2) as executor:
        future_1 = executor.submit(st_calc.calc_prob, metlist_rgb, params)
        future_2 = executor.submit(st_calc.calc_prob, metlist_hb, params)

    rgb_dataframe = future_1.result()
    hb_dataframe = future_2.result()

    print('Calculations finished')

    if params['plot_corner'] == True:
        print('Creating cornerplot....')

        if len(rgb_dataframe.index) > 0:
            st_utils.plot_cornerplot(rgb_dataframe, params, '_RGB')
        if len(hb_dataframe.index) > 0:
            st_utils.plot_cornerplot(hb_dataframe, params, '_HB')

        print('Cornerplot saved under path: '+str(params['save_path']))

    # calculate probability of HB or RGB evolutionary stage
    rgb_prob = rgb_dataframe['posterior_weight'].sum(
    )/(rgb_dataframe['posterior_weight'].sum()+hb_dataframe['posterior_weight'].sum())
    hb_prob = hb_dataframe['posterior_weight'].sum(
    )/(rgb_dataframe['posterior_weight'].sum()+hb_dataframe['posterior_weight'].sum())

    if params('return_ascii') == True:
        print('Saving output file '+params['object_name']+'RGB.out under '+str(params['save_path']))

print('END')
