# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import os
import sys
import yaml
import h5py
import time
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import ThreadPoolExecutor

__author__ = "Stephan Stock @ ZAH, Landessternwarte Heidelberg"
__version__ = "1.0"
__license__ = "MIT"

import SPOG.sp_calc as sp_calc
import SPOG.sp_utils as sp_utils
import SPOG.sp_plots as sp_plots


def main():

    print(f' \n'
          f'Welcome to SPOG+ V{__version__}\n')
    # Input the parameter set in the parameter file - the list of possible
    print(f'Author: {__author__}\n')
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
                      'evolutionary_stage_prior': ['RGB', 'HB'],
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

    # check if models have been dowloaded, if not download
    if not os.path.isfile(params['model_path']):
        if not os.path.exists(os.path.dirname(params['model_path'])):
            os.makedirs(os.path.dirname(params['model_path']))
        params['model_path'] = sp_utils.download_models(os.path.dirname(params['model_path']))

    # check if output path exists
    if not os.path.exists(params['save_path']):
        while True:
            check = input(f'Output path does not exist yes. Create? \n'
                          f'Type (y)es or (n)o: \n')
            if check.lower() == 'yes' or check.lower() == 'y':
                try:
                    os.makedirs(params['save_path'], exist_ok=True)
                    print("Directory created successfully")
                except OSError as error:
                    print("Directory can not be created")
                break
            elif check.lower() == 'no' or check.lower() == 'n':
                raise Exception('Sorry, without a path to put the output files I cannot continue!')
            else:
                print("Sorry, I didn't understand that.")
                continue

    # if classic mode is used... posterior has to be plotted
    if params['mode'] == 'classic':
        params['plot_posterior'] = True
        params['parameterization'] = 'default'
        print('WARNING: You have chosen the classic mode. This requires to plot the posterior and to use the default parametrization! Setting "plot_posterior": True and "parametrization": default')

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

    # run some checks
    if 'RGB' not in params['evolutionary_stage_prior'] and 'HB' not in params['evolutionary_stage_prior'] and 'MS' not in params['evolutionary_stage_prior'] and 'PMS' not in params['evolutionary_stage_prior']:
        raise Exception(
            'No valid evolutionary stages applied. Please check the param.yaml file. Allowed are: PMS, MS, RGB, HB')
    if 'V_johnson' not in params['photometric_band_A'] and 'B_johnson' not in params['photometric_band_A'] and 'I_johnson' not in params['photometric_band_A'] and 'J_2mass' not in params['photometric_band_A'] and 'H_2mass' not in params['photometric_band_A'] and 'Ks_2mass' not in params['photometric_band_A'] and 'G_gaia' not in params['photometric_band_A'] and 'G_BP_gaia' not in params['photometric_band_A'] and 'G_RP_gaia' not in params['photometric_band_A']:
        raise Exception('No valid photometric_band_A band applied. Please check the param.yaml file. Allowed are: V_johnson, B_johnson, I_johnson, J_2mass, H_2mass, Ks_2mass, G_gaia, G_BP_gaia, G_RP_gaia ')
    if 'V_johnson' not in params['photometric_band_B'] and 'B_johnson' not in params['photometric_band_B'] and 'I_johnson' not in params['photometric_band_B'] and 'J_2mass' not in params['photometric_band_B'] and 'H_2mass' not in params['photometric_band_B'] and 'Ks_2mass' not in params['photometric_band_B'] and 'G_gaia' not in params['photometric_band_B'] and 'G_BP_gaia' not in params['photometric_band_B'] and 'G_RP_gaia' not in params['photometric_band_B']:
        raise Exception('No valid photometric_band_B band applied. Please check the param.yaml file. Allowed are: V_johnson, B_johnson, I_johnson, J_2mass, H_2mass, Ks_2mass, G_gaia, G_BP_gaia, G_RP_gaia ')
    if 'default' not in params['parameterization'] and 'default2' not in params['parameterization'] and 'linear' not in params['parameterization'] and 'log' not in params['parameterization']:
        raise Exception(
            'No valid parametrization applied. Please check the param.yaml file. Allowed are: default, default2, log, linear')
    if 'new' not in params['mode'] and 'classic' not in params['mode']:
        raise Exception(
            'No valid mode applied. Please check the param.yaml file. Allowed are: new, classic')
    if params['met_star'][0] > np.max(Fe_mod) or params['met_star'][0] < np.min(Fe_mod):
        raise Exception(
            'Metallicity of star outside of model grid. No reliable output can be generated. I am sorry.')

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

    if maxfeh > np.max(Fe_mod) or minfeh < np.min(Fe_mod):
        print('WARNING: Range of 5-Sigma uncertainty of stellar metallicity ['+str(np.round(minfeh, 2))+','+str(
            np.round(maxfeh, 2)) + '] exceeds available model grid['+str(np.round(np.min(Fe_mod), 2))+', '+str(np.round(np.max(Fe_mod), 2))+']. Result might not be reliable. ')

    metallicities_to_load = df_all_met_sort.loc[(
        df_all_met_sort['Fe/H'] >= minfeh) & (df_all_met_sort['Fe/H'] <= maxfeh)]['folder']
    weights_of_metallicities = df_all_met_sort.loc[(
        df_all_met_sort['Fe/H'] >= minfeh) & (df_all_met_sort['Fe/H'] <= maxfeh)]['weight']

    weights_of_metallicities = weights_of_metallicities/len(weights_of_metallicities)

    metlist_rgb = []
    metlist_hb = []
    metlist_ms = []
    metlist_pms = []
    t0 = time.time()
    print('loading relevant models and calculating posterior')
    for model, weight in tqdm(zip(metallicities_to_load, weights_of_metallicities), total=len(metallicities_to_load)):

        with h5py.File(params['model_path'], 'r') as hdf:
            mass_group_lowmass = hdf.get(model+'/lowmass')
            mass_group_lowmass_items = list(mass_group_lowmass.items())

            mass_group_highmass = hdf.get(model+'/highmass')
            mass_group_highmass_items = list(mass_group_highmass.items())

            mass_group_hb = hdf.get(model+'/hb')
            mass_group_hb_items = list(mass_group_hb.items())

            mass_group_ms = hdf.get(model+'/ms')
            mass_group_ms_items = list(mass_group_ms.items())

            mass_group_pms = hdf.get(model+'/pms')
            mass_group_pms_items = list(mass_group_pms.items())

            if 'RGB' in params['evolutionary_stage_prior']:
                metlist_rgb.append(sp_utils.load_models(hdf, mass_group_lowmass_items,
                                                        model, weight, params, 6., 10.999, 'lowmass'))
                metlist_rgb.append(sp_utils.load_models(hdf, mass_group_highmass_items,
                                                        model, weight, params, 6., 10.999, 'highmass'))
            if 'HB' in params['evolutionary_stage_prior']:
                metlist_hb.append(sp_utils.load_models(hdf, mass_group_hb_items,
                                                       model, weight, params, 1., 5., 'hb'))
                metlist_hb.append(sp_utils.load_models(hdf, mass_group_highmass_items,
                                                       model, weight, params, 11., 14.999, 'highmass'))
            if 'MS' in params['evolutionary_stage_prior']:
                metlist_ms.append(sp_utils.load_models(hdf, mass_group_ms_items,
                                                       model, weight, params, 5., 6., 'ms'))
            if 'PMS' in params['evolutionary_stage_prior']:
                metlist_pms.append(sp_utils.load_models(hdf, mass_group_pms_items,
                                                        model, weight, params, 1., 5., 'pms'))

    frames = []

    with ProcessPoolExecutor(4) as executor:
        if 'RGB' in params['evolutionary_stage_prior'] and len(metlist_rgb) > 0:
            future_1 = executor.submit(sp_calc.calc_prob, pd.concat(metlist_rgb), params)
        if 'HB' in params['evolutionary_stage_prior'] and len(metlist_hb) > 0:
            future_2 = executor.submit(sp_calc.calc_prob, pd.concat(metlist_hb), params)
        if 'MS' in params['evolutionary_stage_prior'] and len(metlist_ms) > 0:
            future_3 = executor.submit(sp_calc.calc_prob, pd.concat(metlist_ms), params)
        if 'PMS' in params['evolutionary_stage_prior'] and len(metlist_pms) > 0:
            future_4 = executor.submit(sp_calc.calc_prob, pd.concat(metlist_pms), params)

    if 'RGB' in params['evolutionary_stage_prior'] and len(metlist_rgb) > 0:
        rgb_dataframe = future_1.result()
        frames.append(rgb_dataframe)
    else:
        rgb_dataframe = pd.DataFrame(columns=['posterior_weight'], dtype=object)
    if 'HB' in params['evolutionary_stage_prior'] and len(metlist_hb) > 0:
        hb_dataframe = future_2.result()
        frames.append(hb_dataframe)
    else:
        hb_dataframe = pd.DataFrame(columns=['posterior_weight'], dtype=object)
    if 'MS' in params['evolutionary_stage_prior'] and len(metlist_ms) > 0:
        ms_dataframe = future_3.result()
        frames.append(ms_dataframe)
    else:
        ms_dataframe = pd.DataFrame(columns=['posterior_weight'], dtype=object)
    if 'PMS' in params['evolutionary_stage_prior'] and len(metlist_pms) > 0:
        pms_dataframe = future_4.result()
        frames.append(pms_dataframe)
    else:
        pms_dataframe = pd.DataFrame(columns=['posterior_weight'], dtype=object)

    all_dataframe = pd.concat(frames)

    t1 = time.time()
    total = t1-t0
    print('Models loaded and calculations finished in '+str(round(total, 2))+' s \n')

    # calculate probability under prior of each evolutionary stage
    rgb_prob = rgb_dataframe['posterior_weight'].sum()/(rgb_dataframe['posterior_weight'].sum(
    )+hb_dataframe['posterior_weight'].sum()+ms_dataframe['posterior_weight'].sum()+pms_dataframe['posterior_weight'].sum())
    hb_prob = hb_dataframe['posterior_weight'].sum()/(rgb_dataframe['posterior_weight'].sum(
    )+hb_dataframe['posterior_weight'].sum()+ms_dataframe['posterior_weight'].sum()+pms_dataframe['posterior_weight'].sum())
    ms_prob = ms_dataframe['posterior_weight'].sum()/(rgb_dataframe['posterior_weight'].sum(
    )+hb_dataframe['posterior_weight'].sum()+ms_dataframe['posterior_weight'].sum()+pms_dataframe['posterior_weight'].sum())
    pms_prob = pms_dataframe['posterior_weight'].sum()/(rgb_dataframe['posterior_weight'].sum(
    )+hb_dataframe['posterior_weight'].sum()+ms_dataframe['posterior_weight'].sum()+pms_dataframe['posterior_weight'].sum())

    if hb_prob < 0.005:
        print('WARNING: Horizontal branch solution below 0.5% probability. No extra output will be created for this solution type!')
    if rgb_prob < 0.005:
        print('WARNING: Red giant branch solution below 0.5% probability. No extra output will be created for this solution type!')
    if ms_prob < 0.005:
        print('WARNING: Main Sequence solution below 0.5% probability. No extra output will be created for this solution type!')
    if pms_prob < 0.005:
        print('WARNING: Pre-main sequence solution below 0.5% probability. No extra output will be created for this solution type!')

    if params['plot_corner'] == True:
        print('Creating cornerplot...\n')

        if len(rgb_dataframe.index) > 0 and rgb_prob > 0.005:
            sp_plots.plot_cornerplot(rgb_dataframe, params, '_RGB')
        if len(hb_dataframe.index) > 0 and hb_prob > 0.005:
            sp_plots.plot_cornerplot(hb_dataframe, params, '_HB')
        if len(ms_dataframe.index) > 0 and ms_prob > 0.005:
            sp_plots.plot_cornerplot(ms_dataframe, params, '_MS')
        if len(pms_dataframe.index) > 0 and pms_prob > 0.005:
            sp_plots.plot_cornerplot(pms_dataframe, params, '_PMS')
        if len(all_dataframe.index) > 0:
            sp_plots.plot_cornerplot(all_dataframe, params, '_ALL')

        print('Cornerplots saved under path: '+str(params['save_path'])+'\n')

    if params['plot_posterior'] == True:
        print('Creating Posterior plot...\n')
        if len(rgb_dataframe.index) > 0 and rgb_prob > 0.005:
            sp_plots.plot_posterior(rgb_dataframe, params, '_RGB')
        if len(hb_dataframe.index) > 0 and hb_prob > 0.005:
            sp_plots.plot_posterior(hb_dataframe, params, '_HB')
        if len(ms_dataframe.index) > 0 and ms_prob > 0.005:
            sp_plots.plot_posterior(ms_dataframe, params, '_MS')
        if len(pms_dataframe.index) > 0 and pms_prob > 0.005:
            sp_plots.plot_posterior(pms_dataframe, params, '_PMS')
        if len(all_dataframe.index) > 0:
            sp_plots.plot_posterior(all_dataframe, params, '_ALL')

        print('Posterior plots saved under path: '+str(params['save_path'])+'\n')

    if params['return_ascii'] == True:
        if len(rgb_dataframe.index) > 0 and rgb_prob > 0.005:
            print('Saving parameter output file "'+params['object_name'] +
                  '_RGB.out" under path: '+str(params['save_path'])+'\n')
            sp_utils.write_outputfile(rgb_dataframe, params, '_RGB', rgb_prob)
        if len(hb_dataframe.index) > 0 and hb_prob > 0.005:
            print('Saving parameter output file "'+params['object_name'] +
                  '_HB.out" under path: '+str(params['save_path'])+'\n')
            sp_utils.write_outputfile(hb_dataframe, params, '_HB', hb_prob)
        if len(ms_dataframe.index) > 0 and ms_prob > 0.005:
            print('Saving parameter output file "'+params['object_name'] +
                  '_MS.out" under path: '+str(params['save_path'])+'\n')
            sp_utils.write_outputfile(ms_dataframe, params, '_MS', ms_prob)
        if len(pms_dataframe.index) > 0 and pms_prob > 0.005:
            print('Saving parameter output file "'+params['object_name'] +
                  '_PMS.out" under path: '+str(params['save_path'])+'\n')
            sp_utils.write_outputfile(pms_dataframe, params, '_PMS', pms_prob)
        if len(all_dataframe.index) > 0:
            print('Saving parameter output file "'+params['object_name'] +
                  '_ALL.out" under path: '+str(params['save_path'])+'\n')
            sp_utils.write_outputfile(all_dataframe, params, '_ALL', 1.0)

    if params['save_posterior'] == True:
        print('Saving posterior samples in HDF5 file "' +
              params['object_name']+'_posteriors.h5" under path: '+str(params['save_path'])+'\n')
        rgb_dataframe.to_hdf(params['save_path']+params['object_name'] +
                             '_posteriors.h5', key='RGB', mode='w')
        hb_dataframe.to_hdf(params['save_path']+params['object_name'] +
                            '_posteriors.h5', key='HB')
        ms_dataframe.to_hdf(params['save_path']+params['object_name'] +
                            '_posteriors.h5', key='MS', mode='w')
        pms_dataframe.to_hdf(params['save_path']+params['object_name'] +
                             '_posteriors.h5', key='PMS')
        all_dataframe.to_hdf(params['save_path']+params['object_name'] +
                             '_posteriors.h5', key='ALL')

    print('END')


if __name__ == '__main__':
    main()
