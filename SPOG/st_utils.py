import numpy as np
import pandas as pd
import h5py
import corner
Rsun = 6.957e10
__author__ = "Stephan Stock @ ZAH, Landessternwarte Heidelberg"
__version__ = "0.9"
__license__ = "MIT"


def load_models(hdf, group, model, weight, params, phase_low, phase_up, string):
    list = []
    # breakpoint()
    for i in range(len(group)):
        dataset = hdf[model+'/'+string+'/'+str(group[i][0])][()]
        attrs = hdf[model+'/'+string+'/' +
                    str(group[i][0])].attrs['header']
        df = pd.DataFrame(dataset, columns=attrs.split())
        df['ABL'] = 10**(0.2*df['mag'])
        df['met_weight'] = weight
        # use only 5 sigma range of color and ABL around star
        # (faster calculation and less ram use)
        df_load = df.loc[(df['phase'] >= phase_low) & (df['phase'] <= phase_up) & (df['color'] <= params['color_star']+5*params['color_star_err']) & (df['color'] >=
                                                                                                                                                      params['color_star']-5*params['color_star_err']) & (df['ABL'] >= params['ABL_star']-5*params['ABL_star_err']) & (df['ABL'] <= params['ABL_star']+5*params['ABL_star_err'])]
        list.append(df_load)
    return list


def plot_cornerplot(df, params, sol_type):
    # define solar radius (cgs)

    if params['parameterization'] == 'log':
        df['logM'] = np.log10(df['mass_act'])
        df['logR'] = np.log10(df['Radius'])
        df['log_age'] = np.log10(df['age'])
        figure = corner.corner(df[['logM', 'logR', 'logg_act', 'log_age', 'logL', 'logT']],
                               weights=df['posterior_weight'], quantiles=[0.16, 0.5, 0.84],
                               labels=[r'log($M [M_\odot])$', r'log($R [R_\odot])$', r'log($g[cm/s^2]$)', r'log($\tau$[yr])', r'log($L [L_\odot]$)', 'log(T[K])'], smooth=1.0, smooth1d=1.0, plot_contours=True, fill_contours=True,
                               plot_datapoints=False, show_titles=True)
        figure.savefig(params['save_path']+params['object_name']+sol_type+'.pdf')
    elif params['parameterization'] == 'linear':
        df['g'] = 10**(df['logg_act'])
        df['L'] = 10**(df['logL'])
        df['T'] = 10**(df['logT'])
        figure = corner.corner(df[['mass_act', 'Radius', 'g', 'age', 'L', 'T']],
                               weights=df['posterior_weight'], quantiles=[0.16, 0.5, 0.84],
                               labels=[r'Mass $[M_\odot]$', r'Radius $[R_\odot]$', r'$g[cm/s^2]$', r'$\tau$[yr]', r'$L [L_\odot]$', 'T[K]'], smooth=1.0, smooth1d=1.0, plot_contours=True, fill_contours=True,
                               plot_datapoints=False, show_titles=True)
        figure.savefig(params['save_path']+params['object_name']+sol_type+'.pdf')
    elif params['parameterization'] == 'default2':
        df['Radius'] = (10**df['Radius'])/Rsun
        df['L'] = 10**(df['logL'])
        df['T'] = 10**(df['logT'])
        df['log_age'] = np.log10(df['age'])
        df['mass_loss'] = pd.Series.abs(df['mass_ZAMS']-df['mass_act'])
        figure = corner.corner(df[['mass_act', 'mass_loss', 'Radius', 'logg_act', 'log_age', 'L', 'T', 'phase']],
                               weights=df['posterior_weight'], quantiles=[0.16, 0.5, 0.84],
                               labels=[r'Mass $[M_\odot]$', r'Mass loss $[M_\odot]$', r'Radius $[R_\odot]$', r'log($g[cm/s^2]$)', r'log($\tau$[yr])', r'$L [L_\odot]$', 'T[K]', 'Phase'], smooth=1.0, smooth1d=1.0, plot_contours=True, fill_contours=True,
                               plot_datapoints=False, show_titles=True)
        figure.savefig(params['save_path']+params['object_name']+sol_type+'.pdf')
    else:
        df['log_age'] = np.log10(df['age'])
        df['mass_loss'] = pd.Series.abs(df['mass_ZAMS']-df['mass_act'])
        df['Radius'] = (10**df['Radius'])/Rsun
        figure = corner.corner(df[['mass_act', 'mass_loss', 'Radius', 'logg_act', 'log_age', 'logL', 'logT', 'phase']],
                               weights=df['posterior_weight'], quantiles=[0.16, 0.5, 0.84],
                               labels=[r'Mass $[M_\odot]$', r'Mass loss $[M_\odot]$', r'Radius $[R_\odot]$', r'log($g[cm/s^2]$)', r'log($\tau$[yr])', r'log($L [L_\odot]$)', 'log(T[K])', 'Phase'], smooth=1.0, smooth1d=1.0, plot_contours=True, fill_contours=True,
                               plot_datapoints=False, show_titles=True)
        figure.savefig(params['save_path']+params['object_name']+sol_type+'.pdf')


def write_outputfile(arg):
    pass
