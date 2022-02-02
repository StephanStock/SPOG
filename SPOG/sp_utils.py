import numpy as np
import pandas as pd
import h5py
import corner
import copy
import SPOG.sp_plots as sp_plots
import SPOG.sp_calc as sp_calc


Rsun = 6.957e10
__author__ = "Stephan Stock @ ZAH, Landessternwarte Heidelberg"
__version__ = "1.0"
__license__ = "MIT"


def download(url, filename):
    """Starts a download using the request package.

    Parameters
    ----------
    url : string
        The url from where the file should be downloaded.
    filename : string
        The name the downloaded file should have.

    Returns
    -------
    string
        Returns path to the downloaded file.

    """
    import requests
    import functools
    import pathlib
    import shutil
    from tqdm.auto import tqdm
    r = requests.get(url, stream=True, allow_redirects=True)
    if r.status_code != 200:
        r.raise_for_status()  # Will only raise for 4xx codes, so...
        raise RuntimeError(f"Request to {url} returned status code {r.status_code}")
    file_size = int(r.headers.get('Content-Length', 0))

    path = pathlib.Path(filename).expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)

    desc = "(Unknown total file size)" if file_size == 0 else ""
    r.raw.read = functools.partial(r.raw.read, decode_content=True)  # Decompress if needed
    with tqdm.wrapattr(r.raw, "read", total=file_size, desc=desc) as r_raw:
        with path.open("wb") as f:
            shutil.copyfileobj(r_raw, f)

    return path


def download_models(out):
    """Downloads the required evolutionary model to the user defined path.

    Parameters
    ----------
    out : string
        User defined path were downloaded models will be stored.

    Returns
    -------
        None
    """
    while True:
        check = input(f'Required evolutionary models were not found or do not exist on the disk.\n'
                      f'Do you want to download the models to the following path: {out} ?\n'
                      f'Type (y)es or (n)o: \n')
        if check.lower() == 'yes' or check.lower() == 'y':
            check2 = input(f'Which version of the models do you want to download?\n'
                           f'(1) Minimal, recommended for testing only (830M) \n'
                           f'(2) Student, a good compromise between accuracy, model load time and disk space usage (11.0G) \n'
                           f'(3) Professional, recommended for scientists who would like to publish their results (54.5G) \n'
                           f'Please type either 1, 2, 3, (c)ancel to abort, or d(etails) get more information about which models to download. \n')
            if check2.lower() == '1':
                print('Downloading minimal number of models...')
                output = out+'/Models_minimal.h5'
                download(
                    "https://heibox.uni-heidelberg.de/d/253b8d99e1324fa2b4f5/files/?p=%2FModels_minimal.h5&dl=1", output)
                print('Download completed!')
                break
            elif check2.lower() == '2':
                print('Downloading medium number of models...')
                output = out+'/Models_student.h5'
                download(
                    "https://heibox.uni-heidelberg.de/d/253b8d99e1324fa2b4f5/files/?p=%2FModels_student.h5&dl=1", output)
                print('Download completed!')
                break
            elif check2.lower() == '3':
                print('Downloading all available models...')
                output = out+'/Models_professional.h5'
                download(
                    "https://heibox.uni-heidelberg.de/d/253b8d99e1324fa2b4f5/files/?p=%2FModels_professional.h5&dl=1", output)
                print('Download completed!')
                break
            elif check2.lower() == 'd' or check.lower() == 'details':
                print(f' (1) The minimal version consists of the following metallicities: Z0.0005, Z0.001, Z0.002, Z0.004, Z0.006, Z0.008, Z0.01, Z0.014, Z0.017, Z0.02, Z0.03, Z0.04, Z0.06.')
                print(f'     The mass grid is 0.05 Msun. \n')
                print('\n')
                print(f' (2) The student version consists of 1/5th of the metallicities provided in the professional version. The mass grid is 0.025 Msun.')
                print(
                    f'     If the uncertainty of the metallicity of your star is larger than 0.1 in [Fe/H] than this grid might be enough, even as a professional user. \n')
                print('\n')
                print(f' (3) The professional version consists of metallicities ranging from Z0.005 to Z0.06 in steps of Z=0.0000125. The mass grid is 0.025 Msun. \n')
                continue
            elif check2.lower() == 'c' or check.lower() == 'cancel':
                raise Exception('Sorry, without models this code does not work!')
            else:
                print("Sorry, I didn't understand that.")
                continue
        elif check.lower() == 'no' or check.lower() == 'n':
            raise Exception('Sorry, without models this code does not work!')
        else:
            print("Sorry, I didn't understand that.")
            continue
    return(output)


def load_models(hdf, group, model, weight, params, phase_low, phase_up, string):
    """loads the relevant evolutionary models for posterior calculations and performs some necessary pre-calculations.

    Parameters
    ----------
    hdf : HDF5 file object
        The HDF5 File Object necessary to load the correct models.
    group : list
        A list of the relevant groups to be loaded for the models.
    model : string
        A string containing the information of the current model metallicity.
    weight : float
        A float containing the value of the posterior weight for the current model metallicity.
    params : dictionairy
        A dictionairy containing information about user requirements.
    phase_low : float
        A float containing the lower evolutionary phase value for the models to be loaded.
    phase_up : float
        A float containing the upper evolutionary phase value for the models to be loaded.
    string : string
        A string containing information about the current evolutionary stage.

    Returns
    -------
    list
        A list containing pandas dataframes with the relevant models for posterior calculations.

    """
    list = []
    for i in range(len(group)):
        dataset = hdf[model+'/'+string+'/'+str(group[i][0])][()]
        attrs = hdf[model+'/'+string+'/' +
                    str(group[i][0])].attrs['header']
        df = pd.DataFrame(dataset, columns=attrs.split())
        if params['reverse'] == True:
            df['color'] = df[params['photometric_band_A']]-df[params['photometric_band_B']]
        else:
            df['color'] = df[params['photometric_band_B']]-df[params['photometric_band_A']]
        df['ABL'] = 10**(0.2*df[params['photometric_band_A']])
        df['met_weight'] = weight
        # use only 5 sigma range of color and ABL around star
        # (faster calculation and less ram use)

        df_load = df.loc[(df['phase'] >= phase_low) & (df['phase'] <= phase_up) & (df['color'] <= params['color_star'][0]+5*params['color_star'][1]) & (df['color'] >=
                                                                                                                                                        params['color_star'][0]-5*params['color_star'][0]) & (df['ABL'] >= params['ABL_star']-5*params['ABL_star_err']) & (df['ABL'] <= params['ABL_star']+5*params['ABL_star_err'])]
        list.append(sp_calc.get_mean_track(df_load))
    return pd.concat(list)


def write_outputfile(df, params, sol_type, probability):
    """Short summary.

    Parameters
    ----------
    df : pd.DataFrame
        The solution dataframe with all information about the models, model parameters, and their posterior weights.
    params : dictionairy
        A dictionairy containing information about user requirements.
    sol_type : string
        A string containing information about which evolutionary stage the solution dataframe corresponds to.
    probability : float
        A float providing the probability of the solutions dataframe with respect to the prior evolutionary stages.

    Returns
    -------
        None
    """
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
        if (10**result_list[1])-(10**result_list[0]) < 0.025 or (10**result_list[2])-(10**result_list[1]) < 0.025:
            print("WARNING: Uncertainty of mass for "+sol_type +
                  " solution is smaller than the available mass grid of the models. Uncertainty is probably not reliable.")
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
        if result_list[1]-result_list[0] < 0.025 or result_list[2]-result_list[1] < 0.025:
            print("WARNING: Uncertainty of mass for "+sol_type +
                  " solution is smaller than the available mass grid of the models. Uncertainty is probably not reliable.")
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
        if result_list[1]-result_list[0] < 0.025 or result_list[2]-result_list[1] < 0.025:
            print("WARNING: Uncertainty of mass for "+sol_type +
                  " solution is smaller than the available mass grid of the models. Uncertainty is probably not reliable.")
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
            # shallow copy of gloabl result variable
            result_list = copy.copy(sp_plots.get_global_list_uncertainty())
            result_list = result_list[0:18]
            result_list.extend([probability])
            if result_list[1]-result_list[0] < 0.025 or result_list[2]-result_list[1] < 0.025:
                print("WARNING: Uncertainty of mass for "+sol_type +
                      " solution is smaller than the available mass grid of the models. Uncertainty is probably not reliable.")
            np.savetxt(params['save_path']+params['object_name']+sol_type+'.dat', np.array(result_list).reshape((1, 19)),
                       header='Mass_neg_[Msun] Mass_mode_[Msun] Mass_pos_[Msun] Radius_neg_[Rsun] Radius_mode_[Rsun] Radius_pos_[Rsun] logg_neg_[cgs] logg_mode_[cgs] logg_pos_[cgs] log_Age_neg_[log_yr] log_Age_mode_[log_yr] log_Age_pos_[log_yr] Luminosity_neg_[Lsun] Luminosity_mode_[Lsun] Luminosity_pos_[Lsun] Temperature_neg_[K] Temperature_mode_[K] Temperature_pos_[K] Probability'+sol_type)

        else:
            for i in parameters:
                result_list.extend(corner.quantile(
                    df[i], [0.16, 0.5, 0.84], weights=df['posterior_weight']))
            result_list.extend([probability])
            if result_list[1]-result_list[0] < 0.025 or result_list[2]-result_list[1] < 0.025:
                print("WARNING: Uncertainty of mass for "+sol_type +
                      " solution is smaller than the available mass grid of the models. Uncertainty is probably not reliable.")
            np.savetxt(params['save_path']+params['object_name']+sol_type+'.dat', np.array(result_list).reshape((1, 19)),
                       header='Mass_q16_[Msun] Mass_q50_[Msun] Mass_q84_[Msun] Radius_q16_[Rsun] Radius_q50_[Rsun] Radius_q84_[Rsun] logg_q16_[cgs] logg_q50_[cgs] logg_q84_[cgs] log_Age_q16_[log_yr] log_Age_q50_[log_yr] log_Age_q84_[log_yr] Luminosity_q16_[Lsun] Luminosity_q50_[Lsun] Luminosity_q84_[Lsun] Temperature_q16_[K] Temperature_q50_[K] Temperature_q84_[K] Probability'+sol_type)

    pass
