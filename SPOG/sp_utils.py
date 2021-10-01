import numpy as np
import pandas as pd
import h5py
import corner
Rsun = 6.957e10
__author__ = "Stephan Stock @ ZAH, Landessternwarte Heidelberg"
__version__ = "0.92"
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
            df['color'] = df[params['photometric_band_B']]  # for the moment
        else:
            df['color'] = df[params['photometric_band_B']]

        df['ABL'] = 10**(0.2*df['mag'])
        df['met_weight'] = weight
        # use only 5 sigma range of color and ABL around star
        # (faster calculation and less ram use)
        df_load = df.loc[(df['phase'] >= phase_low) & (df['phase'] <= phase_up) & (df['color'] <= params['color_star']+5*params['color_star_err']) & (df['color'] >=
                                                                                                                                                      params['color_star']-5*params['color_star_err']) & (df['ABL'] >= params['ABL_star']-5*params['ABL_star_err']) & (df['ABL'] <= params['ABL_star']+5*params['ABL_star_err'])]
        list.append(df_load)
    return list


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
        np.savetxt(params['save_path']+params['object_name']+sol_type+'.dat', np.array(result_list).reshape((1, 18)),
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
        np.savetxt(params['save_path']+params['object_name']+sol_type+'.dat', np.array(result_list).reshape((1, 18)),
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
        for i in parameters:
            result_list.extend(corner.quantile(
                df[i], [0.16, 0.5, 0.84], weights=df['posterior_weight']))
        result_list.extend([probability])
        np.savetxt(params['save_path']+params['object_name']+sol_type+'.dat', np.array(result_list).reshape((1, 18)),
                   header='Mass_q16_[Msun] Mass_q50_[Msun] Mass_q84_[Msun] Radius_q16_[Rsun] Radius_q50_[Rsun] Radius_q84_[Rsun] logg_q16_[cgs] logg_q50_[cgs] logg_q84_[cgs] log_Age_q16_[log_yr] log_Age_q50_[log_yr] log_Age_q84_[log_yr] Luminosity_q16_[Lsun] Luminosity_q50_[Lsun] Luminosity_q84_[Lsun] Temperature_q16_[K] Temperature_q50_[K] Temperature_q84_[K] Probability'+sol_type)

    pass
