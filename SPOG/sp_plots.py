import numpy as np
import pandas as pd
import corner
import matplotlib.pyplot as plt
import astropy.stats
import sp_calc
try:
    from scipy.ndimage import gaussian_filter
except ImportError:
    gaussian_filter = None
Rsun = 6.957e10
__author__ = "Stephan Stock @ ZAH, Landessternwarte Heidelberg"
__version__ = "0.92"
__license__ = "MIT"

global_list_uncertainty = []


def plot_cornerplot(df, params, sol_type):
    # define solar radius (cgs)

    if params['parameterization'] == 'log':
        df['logM'] = np.log10(df['mass_act'])
        df['log_age'] = np.log10(df['age'])
        df['logRR'] = np.log((10**df['logR'])/Rsun)  # transform from logR in cgs to logR in Rsun
        figure = corner.corner(df[['logM', 'logRR', 'logg_act', 'log_age', 'logL', 'logT']],
                               weights=df['posterior_weight'], quantiles=[0.16, 0.5, 0.84],
                               labels=[r'log($M [M_\odot])$', r'log($R [R_\odot])$', r'log($g[cm/s^2]$)', r'log($\tau$[yr])', r'log($L [L_\odot]$)', 'log(T[K])'], smooth=1.0, smooth1d=1.0, plot_contours=True, fill_contours=True,
                               plot_datapoints=False, show_titles=True)
        figure.savefig(params['save_path']+params['object_name']+sol_type+'_corner.pdf')
    elif params['parameterization'] == 'linear':
        df['g'] = 10**(df['logg_act'])
        df['L'] = 10**(df['logL'])
        df['T'] = 10**(df['logT'])
        df['Radius'] = (10**df['logR'])/Rsun
        df['age_Myr'] = df['age']/1e6
        figure = corner.corner(df[['mass_act', 'Radius', 'g', 'age_Myr', 'L', 'T']],
                               weights=df['posterior_weight'], quantiles=[0.16, 0.5, 0.84],
                               labels=[r'M $[M_\odot]$', r'R $[R_\odot]$', r'$g[cm/s^2]$', r'$\tau$[Myr]', r'$L [L_\odot]$', 'T[K]'], smooth=1.0, smooth1d=1.0, plot_contours=True, fill_contours=True,
                               plot_datapoints=False, show_titles=True)
        figure.savefig(params['save_path']+params['object_name']+sol_type+'_corner.pdf')
    elif params['parameterization'] == 'default2':
        df['Radius'] = (10**df['logR'])/Rsun
        df['L'] = 10**(df['logL'])
        df['T'] = 10**(df['logT'])
        df['log_age'] = np.log10(df['age'])
        df['mass_loss'] = pd.Series.abs(df['mass_ZAMS']-df['mass_act'])
        figure = corner.corner(df[['mass_act', 'mass_loss', 'Radius', 'logg_act', 'log_age', 'L', 'T', 'phase']],
                               weights=df['posterior_weight'], quantiles=[0.16, 0.5, 0.84],
                               labels=[r'M $[M_\odot]$', r'Mass loss $[M_\odot]$', r'R $[R_\odot]$', r'log($g[cm/s^2]$)', r'log($\tau$[yr])', r'$L [L_\odot]$', 'T[K]', 'Phase'], smooth=1.0, smooth1d=1.0, plot_contours=True, fill_contours=True,
                               plot_datapoints=False, show_titles=True)
        figure.savefig(params['save_path']+params['object_name']+sol_type+'_corner.pdf')
    else:
        df['Radius'] = (10**df['logR'])/Rsun
        df['L'] = 10**(df['logL'])
        df['T'] = 10**(df['logT'])
        df['log_age'] = np.log10(df['age'])
        if params['mode'] == 'classic':
            figure = corner.corner(df[['mass_act', 'Radius', 'logg_act', 'log_age', 'L', 'T']],
                                   weights=df['posterior_weight'],
                                   labels=[r'M $[M_\odot]$', r'R $[R_\odot]$', r'log($g[cm/s^2]$)', r'log($\tau$[yr])', r'$L [L_\odot]$', 'T[K]'], smooth=1.0, smooth1d=1.0, plot_contours=True, fill_contours=True,
                                   plot_datapoints=False, show_titles=False)
        else:
            figure = corner.corner(df[['mass_act', 'Radius', 'logg_act', 'log_age', 'L', 'T']],
                                   weights=df['posterior_weight'], quantiles=[0.16, 0.5, 0.84],
                                   labels=[r'M $[M_\odot]$', r'R $[R_\odot]$', r'log($g[cm/s^2]$)', r'log($\tau$[yr])', r'$L [L_\odot]$', 'T[K]'], smooth=1.0, smooth1d=1.0, plot_contours=True, fill_contours=True,
                                   plot_datapoints=False, show_titles=True)
        figure.savefig(params['save_path']+params['object_name']+sol_type+'_corner.pdf')


def plot_posterior(df, params, sol_type):

    if params['parameterization'] == 'log':
        df['logM'] = np.log10(df['mass_act'])
        df['log_age'] = np.log10(df['age'])
        df['logRR'] = np.log((10**df['logR'])/Rsun)  # transform from logR in cgs to logR in Rsun

        if params['posterior_fig_kwargs'] == {}:
            params['posterior_fig_kwargs'] = {'figsize': (6, 9)}
        if params['posterior_plot_kwargs'] == {}:
            params['posterior_plot_kwargs'] = {'color': 'black', 'linewidth': 1}
        recal_bin = False
        if params['posterior_bins'] == -1:
            recal_bin = True

        figure, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(
            3, 2, **params['posterior_fig_kwargs'])

        axes = [ax1, ax2, ax3, ax4, ax5, ax6]
        df_list = [df['logM'], df['logRR'], df['logg_act'], df['log_age'],
                   df['logL'], df['logT']]
        labels = [r'log($M[M_\odot]$)',  r'log($R[R_\odot]$)',
                  r'log($g[cm/s^2]$)', r'log($\tau$[yr])', r'$log(L[L_\odot]$)', 'log(T[K])']
        for ax, dataframe, label in zip(axes, df_list, labels):
            if recal_bin == True:
                params['posterior_bins'] = sp_calc.opt_bin(dataframe, df['posterior_weight'], 50.)
            n, b = astropy.stats.histogram(
                dataframe, bins=params['posterior_bins'], weights=df['posterior_weight'])
            n = gaussian_filter(n, params['smooth'])
            x0 = np.array(list(zip(b[:-1], b[1:]))).flatten()
            y0 = np.array(list(zip(n, n))).flatten()
            y0 = y0/(np.max(y0))  # normalize maximum to 1
            ax.plot(x0, y0, **params['posterior_plot_kwargs'])
            ax.set_xlabel(label, fontsize=12)

            q_16, q_50, q_84 = corner.quantile(dataframe, [0.16, 0.5, 0.84],
                                               weights=df['posterior_weight'])
            q_m, q_p = q_50-q_16, q_84-q_50
            for q in [q_16, q_50, q_84]:
                ax.axvline(q, ls="dashed", color='black')
            ax.set_title(str(label)+' = ${{{' +
                         str(format(round(q_50, 2), '.2f'))+'}}}_{{-{'+str(format(round(q_m, 2), '.2f'))+'}}}^{{+{'+str(format(round(q_p, 2), '.2f'))+'}}}$', fontsize=12)
            ax.set_ylim(bottom=0., top=1.04)
            ax.set_xlim(left=dataframe.min(), right=dataframe.max())

        plt.tight_layout()
        figure.savefig(params['save_path']+params['object_name']+sol_type+'_posterior.pdf')

    elif params['parameterization'] == 'linear':
        df['g'] = 10**(df['logg_act'])
        df['L'] = 10**(df['logL'])
        df['T'] = 10**(df['logT'])
        df['Radius'] = (10**df['logR'])/Rsun
        if params['posterior_fig_kwargs'] == {}:
            params['posterior_fig_kwargs'] = {'figsize': (6, 9)}
        if params['posterior_plot_kwargs'] == {}:
            params['posterior_plot_kwargs'] = {'color': 'black', 'linewidth': 1}
        recal_bin = False
        if params['posterior_bins'] == -1:
            recal_bin = True

        figure, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(
            3, 2, **params['posterior_fig_kwargs'])

        axes = [ax1, ax2, ax3, ax4, ax5, ax6]
        df_list = [df['mass_act'], df['Radius'], df['g'], df['age']/1e6,
                   df['L'], df['T']]
        labels = [r'Mass $[M_\odot]$',  r'Radius $[R_\odot]$',
                  r'Surf. Gravity $[cm/s^2]$', r'Age [Myr]', r'Luminosity $[L_\odot]$', 'Temperature [K]']
        for ax, dataframe, label in zip(axes, df_list, labels):
            if recal_bin == True:
                params['posterior_bins'] = sp_calc.opt_bin(dataframe, df['posterior_weight'], 50.)
            n, b = astropy.stats.histogram(
                dataframe, bins=params['posterior_bins'], weights=df['posterior_weight'])
            n = gaussian_filter(n, params['smooth'])
            x0 = np.array(list(zip(b[:-1], b[1:]))).flatten()
            y0 = np.array(list(zip(n, n))).flatten()
            y0 = y0/(np.max(y0))  # normalize maximum to 1
            ax.plot(x0, y0, **params['posterior_plot_kwargs'])
            ax.set_xlabel(label, fontsize=12)

            q_16, q_50, q_84 = corner.quantile(dataframe, [0.16, 0.5, 0.84],
                                               weights=df['posterior_weight'])
            q_m, q_p = q_50-q_16, q_84-q_50
            for q in [q_16, q_50, q_84]:
                ax.axvline(q, ls="dashed", color='black')
            ax.set_title(str(label)+' = ${{{' +
                         str(format(round(q_50, 2), '.2f'))+'}}}_{{-{'+str(format(round(q_m, 2), '.2f'))+'}}}^{{+{'+str(format(round(q_p, 2), '.2f'))+'}}}$', fontsize=12)
            ax.set_ylim(bottom=0., top=1.04)
            ax.set_xlim(left=dataframe.min(), right=dataframe.max())
        plt.tight_layout()
        figure.savefig(params['save_path']+params['object_name']+sol_type+'_posterior.pdf')

    elif params['parameterization'] == 'default2':
        if params['posterior_fig_kwargs'] == {}:
            params['posterior_fig_kwargs'] = {'figsize': (6, 12)}
        if params['posterior_plot_kwargs'] == {}:
            params['posterior_plot_kwargs'] = {'color': 'black', 'linewidth': 1}
        recal_bin = False
        if params['posterior_bins'] == -1:
            recal_bin = True

        figure, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(
            4, 2, **params['posterior_fig_kwargs'])

        axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]
        df_list = [df['mass_act'], df['Radius'], df['logg_act'], df['log_age'],
                   df['L'], df['T'], df['mass_loss'], df['phase']]
        labels = [r'Mass $[M_\odot]$',  r'Radius $[R_\odot]$',
                  r'log($g[cm/s^2]$)', r'log($\tau$[yr])', r'$L [L_\odot]$', 'T[K]', r'Mass loss $[M_\odot]$', 'Phase']
        for ax, dataframe, label in zip(axes, df_list, labels):
            if recal_bin == True:
                params['posterior_bins'] = sp_calc.opt_bin(dataframe, df['posterior_weight'], 50.)
            n, b = astropy.stats.histogram(
                dataframe, bins=params['posterior_bins'], weights=df['posterior_weight'])

            n = gaussian_filter(n, params['smooth'])

            x0 = np.array(list(zip(b[:-1], b[1:]))).flatten()
            y0 = np.array(list(zip(n, n))).flatten()
            y0 = y0/(np.max(y0))  # normalize maximum to 1
            ax.plot(x0, y0, **params['posterior_plot_kwargs'])

            ax.set_xlabel(label, fontsize=12)

            q_16, q_50, q_84 = corner.quantile(dataframe, [0.16, 0.5, 0.84],
                                               weights=df['posterior_weight'])
            q_m, q_p = q_50-q_16, q_84-q_50
            for q in [q_16, q_50, q_84]:
                ax.axvline(q, ls="dashed", color='black')
            ax.set_title(str(label)+' = ${{{' +
                         str(format(round(q_50, 2), '.2f'))+'}}}_{{-{'+str(format(round(q_m, 2), '.2f'))+'}}}^{{+{'+str(format(round(q_p, 2), '.2f'))+'}}}$', fontsize=12)
            ax.set_ylim(bottom=0., top=1.04)
            ax.set_xlim(left=dataframe.min(), right=dataframe.max())

        plt.tight_layout()
        figure.savefig(params['save_path']+params['object_name']+sol_type+'_posterior.pdf')
    else:
        if params['posterior_fig_kwargs'] == {}:
            params['posterior_fig_kwargs'] = {'figsize': (6, 9)}
        if params['posterior_plot_kwargs'] == {}:
            params['posterior_plot_kwargs'] = {'color': 'black', 'linewidth': 1}
        recal_bin = False
        if params['posterior_bins'] == -1:
            recal_bin = True

        figure, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(
            3, 2, **params['posterior_fig_kwargs'])

        axes = [ax1, ax2, ax3, ax4, ax5, ax6]
        df_list = [df['mass_act'], df['Radius'], df['logg_act'], df['log_age'],
                   df['L'], df['T']]
        labels = [r'Mass $[M_\odot]$',  r'Radius $[R_\odot]$',
                  r'log($g[cm/s^2]$)', r'log($\tau$[yr])', r'$L [L_\odot]$', 'T[K]']
        for ax, dataframe, label in zip(axes, df_list, labels):
            if recal_bin == True:
                params['posterior_bins'] = sp_calc.opt_bin(dataframe, df['posterior_weight'], 50.)
            n, b = astropy.stats.histogram(
                dataframe, bins=params['posterior_bins'], weights=df['posterior_weight'])
            n = gaussian_filter(n, params['smooth'])
            x0 = np.array(list(zip(b[:-1], b[1:]))).flatten()
            y0 = np.array(list(zip(n, n))).flatten()
            y0 = y0/(np.max(y0))  # normalize maximum to 1
            ax.plot(x0, y0, **params['posterior_plot_kwargs'])
            ax.set_xlabel(label, fontsize=12)
            if params['mode'] == 'classic':
                x_new, y_new = sp_calc.spline(b, n)
                mode, e_low, e_high = sp_calc.get_mode_uncertainty(x_new, y_new)
                ax.plot(x_new, y_new, color='red')

                ax.set_ylim(bottom=0., top=1.04)
                ax.set_xlim(left=x_new.min(), right=x_new.max())
                for q in [e_low, mode, e_high]:
                    ax.axvline(q, ls="dashed", color='black')
                global_list_uncertainty.extend([e_low, mode, e_high])
                e_m, e_p = mode-e_low, e_high-mode
                ax.set_title(str(label)+' = ${{{' +
                             str(format(round(mode, 2), '.2f'))+'}}}_{{-{'+str(format(round(e_m, 2), '.2f'))+'}}}^{{+{'+str(format(round(e_p, 2), '.2f'))+'}}}$', fontsize=12)

            else:
                q_16, q_50, q_84 = corner.quantile(dataframe, [0.16, 0.5, 0.84],
                                                   weights=df['posterior_weight'])
                q_m, q_p = q_50-q_16, q_84-q_50
                for q in [q_16, q_50, q_84]:
                    ax.axvline(q, ls="dashed", color='black')

                ax.set_title(str(label)+' = ${{{' +
                             str(format(round(q_50, 2), '.2f'))+'}}}_{{-{'+str(format(round(q_m, 2), '.2f'))+'}}}^{{+{'+str(format(round(q_p, 2), '.2f'))+'}}}$', fontsize=12)
                ax.set_ylim(bottom=0., top=1.04)
                ax.set_xlim(left=dataframe.min(), right=dataframe.max())

        plt.tight_layout()
        figure.savefig(params['save_path']+params['object_name']+sol_type+'_posterior.pdf')
    pass


def get_global_list_uncertainty():
    return global_list_uncertainty
