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
__version__ = "1.0"
__license__ = "MIT"

import sp_calc
import sp_utils
import sp_plots


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
                  'photometric_band_A': 'V_johnson',
                  'photometric_band_B': 'B_johnson',
                  'reverse': False,
                  'mag_star': [9, 0.1],
                  'color_star': [0.97, 0.14],
                  'met_star': [0, 0.1],
                  'par_star': [2.45, 0.01],
                  'par_err_dom': False,
                  'use_extinction': False,
                  'A_lambda': 0,
                  'E_color': 0,
                  'parameterization': 'default2',
                  'mode': 'new',
                  'smooth': 1.0,
                  'model_sampling': 1,
                  'plot_corner': True,
                  'return_ascii': True,
                  'plot_posterior': True,
                  'posterior_bins': 20,
                  'posterior_fig_kwargs': {},
                  'posterior_plot_kwargs': {},
                  'save_posterior': True}
params = {}
for key in default_params:
    try:
        params[key] = input_params[key]
    except KeyError:
        params[key] = default_params[key]
        print(f'Parameter {key} is missing in the input file.\n'
              f'Value set to default: {key} = {default_params[key]}\n')


if params['mode'] == 'classic':  # if classic mode is used... posterior has to be plotted
    params['plot_posterior'] = True
    params['parameterization'] = 'default'
    print('Warning: You have chosen the classic mode. This requires to plot the posterior and to use the default parametrization! Setting "plot_posterior": True and "parametrization": default')

# calculate ABL from parallax and magnitude (This parametrization is key to be less biased!)
if params['use_extinction'] == False:
    params['ABL_star'] = 10**(params['mag_star'][0]*0.2+1)*(params['par_star'][0]/1000.)
    # Assumes parallax error dominates photometric uncertainties!
    if params['par_err_dom'] == True:
        params['ABL_star_err'] = 10**(params['mag_star'][0]*0.2+1)*(params['par_star'][1]/1000.)
    # Assumes magnitude uncertainty is not negligible! Slight bias due to non-linear transformation between magnitude and ABL
    else:
        params['ABL_star_err'] = 10**(params['mag_star'][0]*0.2+1)*(params['par_star'][1]/1000.) + \
            params['mag_star'][1]*4.460517*(params['par_star'][0]/1000) * \
            np.exp(0.460517*params['mag_star'][0])
else:
    print('Applying extinction and reddening...')
    params['ABL_star'] = 10**((params['mag_star'][0]-params['A_lambda'][0])
                              * 0.2+1)*(params['par_star'][0]/1000.)
    if params['par_err_dom'] == True:
        params['ABL_star_err'] = 10**(params['mag_star'][0]*0.2+1)*(params['par_star'][1]/1000.)
    else:
        params['ABL_star_err'] = 10**(params['mag_star'][0]*0.2+1)*(params['par_star'][1]/1000.) + \
            np.sqrt(params['mag_star'][1]**2+params['A_lambda'][1]**2)*4.460517*(params['par_star'][0]/1000) * \
            np.exp(0.460517*params['mag_star'][0])
    params['color_star'][0] = params['color_star'][0]-params['E_color'][0]
    params['color_star'][1] = np.sqrt(params['color_star'][1]**2+params['E_color'][1]**2)

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


maxfeh = params['met_star'][0]+5*params['met_star'][1]
minfeh = params['met_star'][0]-5*params['met_star'][1]

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
            metlist_rgb.append(sp_utils.load_models(hdf, mass_group_lowmass_items,
                                                    model, weight, params, 6., 10.999, 'lowmass'))
            metlist_rgb.append(sp_utils.load_models(hdf, mass_group_highmass_items,
                                                    model, weight, params, 6., 10.999, 'highmass'))
            metlist_hb.append(sp_utils.load_models(hdf, mass_group_hb_items,
                                                   model, weight, params, 1., 5., 'hb'))
            metlist_hb.append(sp_utils.load_models(hdf, mass_group_highmass_items,
                                                   model, weight, params, 11., 14.999, 'highmass'))

        print('loaded models of metallicity: '+model)

    t1 = time.time()
    total = t1-t0
    print('Time to load models was '+str(round(total, 2))+' s \n')

    print('All necessary models loaded, starting calculations...\n')

    with ProcessPoolExecutor(2) as executor:
        future_1 = executor.submit(sp_calc.calc_prob, metlist_rgb, params)
        future_2 = executor.submit(sp_calc.calc_prob, metlist_hb, params)

    rgb_dataframe = future_1.result()
    hb_dataframe = future_2.result()

    print('Calculations finished \n')

    # calculate probability under prior of each evolutionary stage
    rgb_prob = rgb_dataframe['posterior_weight'].sum(
    )/(rgb_dataframe['posterior_weight'].sum()+hb_dataframe['posterior_weight'].sum())
    hb_prob = hb_dataframe['posterior_weight'].sum(
    )/(rgb_dataframe['posterior_weight'].sum()+hb_dataframe['posterior_weight'].sum())

    if hb_prob < 0.005:
        print('Warning: Horizontal branch solution below 0.5% probability. No output will be created for this solution type!')
    if rgb_prob < 0.005:
        print('Warning: Red giant brach solution below 0.5% probability. No output will be created for this solution type!')
    if params['plot_corner'] == True:
        print('Creating cornerplot...\n')

        if len(rgb_dataframe.index) > 0 and rgb_prob > 0.005:
            sp_plots.plot_cornerplot(rgb_dataframe, params, '_RGB')
        if len(hb_dataframe.index) > 0 and hb_prob > 0.005:
            sp_plots.plot_cornerplot(hb_dataframe, params, '_HB')

        print('Cornerplots saved under path: '+str(params['save_path'])+'\n')

    if params['plot_posterior'] == True:
        print('Creating Posterior plot...\n')
        if len(rgb_dataframe.index) > 0 and rgb_prob > 0.005:
            sp_plots.plot_posterior(rgb_dataframe, params, '_RGB')
        if len(hb_dataframe.index) > 0 and hb_prob > 0.005:
            sp_plots.plot_posterior(hb_dataframe, params, '_HB')
        print('Posterior plots saved under path: '+str(params['save_path'])+'\n')

    if params['return_ascii'] == True:
        print('Saving output file '+params['object_name'] +
              'RGB.out under '+str(params['save_path'])+'\n')
        if len(rgb_dataframe.index) > 0 and rgb_prob > 0.005:
            sp_utils.write_outputfile(rgb_dataframe, params, '_RGB', rgb_prob)
        if len(hb_dataframe.index) > 0 and hb_prob > 0.005:
            sp_utils.write_outputfile(hb_dataframe, params, '_HB', hb_prob)

    if params['save_posterior'] == True:
        print('Saving posterior samples in hdf5 file ' +
              params['object_name']+'_posteriors.h5 under path '+str(params['save_path'])+'\n')
        rgb_dataframe.to_hdf(params['save_path']+params['object_name'] +
                             '_posteriors.h5', key='RGB', mode='w')
        hb_dataframe.to_hdf(params['save_path']+params['object_name'] +
                            '_posteriors.h5', key='HB')
print('END')
