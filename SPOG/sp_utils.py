import numpy as np
import pandas as pd
import h5py
import corner
import sp_plots

Rsun = 6.957e10
__author__ = "Stephan Stock @ ZAH, Landessternwarte Heidelberg"
__version__ = "1.0"
__license__ = "MIT"


def load_models(hdf, group, model, weight, params, phase_low, phase_up, string):
    list = []
    # breakpoint()
    for i in range(len(group)):
        dataset = hdf[model+'/'+string+'/'+str(group[i][0])][()]
        attrs = hdf[model+'/'+string+'/' +
                    str(group[i][0])].attrs['header']
        df = pd.DataFrame(dataset, columns=attrs.split())
        df['mag'] = df[params['photometric_band_A']]
        if params['reverse'] == True:
            df['color'] = df[params['photometric_band_A']]-df[params['photometric_band_B']]
        else:
            df['color'] = df[params['photometric_band_B']]-df[params['photometric_band_A']]

        df['ABL'] = 10**(0.2*df['mag'])
        df['met_weight'] = weight
        # use only 5 sigma range of color and ABL around star
        # (faster calculation and less ram use)
        df_load = df.loc[(df['phase'] >= phase_low) & (df['phase'] <= phase_up) & (df['color'] <= params['color_star'][0]+5*params['color_star'][1]) & (df['color'] >=
                                                                                                                                                        params['color_star'][0]-5*params['color_star'][0]) & (df['ABL'] >= params['ABL_star']-5*params['ABL_star_err']) & (df['ABL'] <= params['ABL_star']+5*params['ABL_star_err'])]
        list.append(df_load)
    return pd.concat(list)


def write_outputfile(df, params, sol_type, probability):
    if params['parameterization'] == 'log':
        df['logM'] = np.log10(df['mass_act'])
        df['log_age'] = np.log10(df['age'])
        df['logRR'] = np.log((10**df['logR'])/Rsun)  # transform from logR in cgs to logR in Rsun
        result_list = []
        parameters = ['logM', 'logRR', 'logg_act', 'log_age', 'logL', 'logT']
        for i in parameters:
            result_list.extend(corner.quantile(
                df[i], [0.16, 0.5, 0.84], weights=df['posterior_weight']))
        result_list.extend([probability])
        np.savetxt(params['save_path']+params['object_name']+sol_type+'.dat', np.array(result_list).reshape((1, 19)),
                   header='log_Mass_q16_[log_Msun] log_Mass_q50_[log_Msun] log_Mass_q84_[log_Msun] log_Radius_q16_[log_Rsun] log_Radius_q50_[log_Rsun] log_Radius_q84_[log_Rsun] logg_q16_[cgs] logg_q50_[cgs] logg_q84_[cgs] log_Age_q16_[log_yr] log_Age_q50_[log_yr] log_Age_q84_[log_yr] log_Luminosity_q16_[log_Lsun] log_Luminosity_q50_[log_Lsun] log_Luminosity_q84_[log_Lsun] log_Temperature_q16_[log_K] log_Temperature_q50_[log_K] log_Temperature_q84_[log_K] Probability'+sol_type)

    elif params['parameterization'] == 'linear':
        df['g'] = 10**(df['logg_act'])
        df['L'] = 10**(df['logL'])
        df['T'] = 10**(df['logT'])
        df['Radius'] = (10**df['logR'])/Rsun
        result_list = []
        parameters = ['mass_act', 'Radius', 'g', 'age', 'L', 'T']
        for i in parameters:
            result_list.extend(corner.quantile(
                df[i], [0.16, 0.5, 0.84], weights=df['posterior_weight']))
        result_list.extend([probability])
        np.savetxt(params['save_path']+params['object_name']+sol_type+'.dat', np.array(result_list).reshape((1, 19)),
                   header='Mass_q16_[Msun] Mass_q50_[Msun] Mass_q84_[Msun] Radius_q16_[Rsun] Radius_q50_[Rsun] Radius_q84_[Rsun] SurfaceGravity_q16_[cgs] SurfaceGravity_q50_[cgs] SurfaceGravity_q84_[cgs] Age_q16_[yr] Age_q50_[yr] Age_q84_[yr] Luminosity_q16_[Lsun] Luminosity_q50_[Lsun] Luminosity_q84_[Lsun] Temperature_q16_[K] Temperature_q50_[K] Temperature_q84_[K] Probability'+sol_type)

    elif params['parameterization'] == 'default2':
        df['Radius'] = (10**df['logR'])/Rsun
        df['L'] = 10**(df['logL'])
        df['T'] = 10**(df['logT'])
        df['log_age'] = np.log10(df['age'])
        result_list = []
        parameters = ['mass_act', 'mass_loss', 'Radius', 'logg_act', 'log_age', 'L', 'T', 'phase']
        for i in parameters:
            result_list.extend(corner.quantile(
                df[i], [0.16, 0.5, 0.84], weights=df['posterior_weight']))
        result_list.extend([probability])
        np.savetxt(params['save_path']+params['object_name']+sol_type+'.dat', np.array(result_list).reshape((1, 25)),
                   header='Mass_q16_[Msun] Mass_q50_[Msun] Mass_q84_[Msun] Massloss_q16_[Msun] Massloss_q50_[Msun] Massloss_q84_[Msun] Radius_q16_[Rsun] Radius_q50_[Rsun] Radius_q84_[Rsun] logg_q16_[cgs] logg_q50_[cgs] logg_q84_[cgs] log_Age_q16_[log_yr] log_Age_q50_[log_yr] log_Age_q84_[log_yr] Luminosity_q16_[Lsun] Luminosity_q50_[Lsun] Luminosity_q84_[Lsun] Temperature_q16_[K] Temperature_q50_[K] Temperature_q84_[K] Phase_q16 Phase_q50 Phase_q84 Probability'+sol_type)

    else:
        df['log_age'] = np.log10(df['age'])
        df['mass_loss'] = pd.Series.abs(df['mass_ZAMS']-df['mass_act'])
        df['Radius'] = (10**df['logR'])/Rsun
        df['L'] = 10**(df['logL'])
        df['T'] = 10**(df['logT'])
        result_list = []
        parameters = ['mass_act', 'Radius', 'logg_act', 'log_age', 'L', 'T']
        if params['mode'] == 'classic':
            result_list = sp_plots.get_global_list_uncertainty()
            result_list.extend([probability])
            np.savetxt(params['save_path']+params['object_name']+sol_type+'.dat', np.array(result_list).reshape((1, 19)),
                       header='Mass_neg_[Msun] Mass_mode_[Msun] Mass_pos_[Msun] Radius_neg_[Rsun] Radius_mode_[Rsun] Radius_pos_[Rsun] logg_neg_[cgs] logg_mode_[cgs] logg_pos_[cgs] log_Age_neg_[log_yr] log_Age_mode_[log_yr] log_Age_pos_[log_yr] Luminosity_neg_[Lsun] Luminosity_mode_[Lsun] Luminosity_pos_[Lsun] Temperature_neg_[K] Temperature_mode_[K] Temperature_pos_[K] Probability'+sol_type)

        else:
            for i in parameters:
                result_list.extend(corner.quantile(
                    df[i], [0.16, 0.5, 0.84], weights=df['posterior_weight']))
            result_list.extend([probability])
            np.savetxt(params['save_path']+params['object_name']+sol_type+'.dat', np.array(result_list).reshape((1, 19)),
                       header='Mass_q16_[Msun] Mass_q50_[Msun] Mass_q84_[Msun] Radius_q16_[Rsun] Radius_q50_[Rsun] Radius_q84_[Rsun] logg_q16_[cgs] logg_q50_[cgs] logg_q84_[cgs] log_Age_q16_[log_yr] log_Age_q50_[log_yr] log_Age_q84_[log_yr] Luminosity_q16_[Lsun] Luminosity_q50_[Lsun] Luminosity_q84_[Lsun] Temperature_q16_[K] Temperature_q50_[K] Temperature_q84_[K] Probability'+sol_type)

    pass
