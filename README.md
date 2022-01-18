# *SPOG* (Stellar Parameters of Giants)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

*SPOG* is a Python script for uncomplicated determination of stellar parameters based on the original unpublished IDL code used and outlined in [Stock et al. (2018)](https://ui.adsabs.harvard.edu/abs/2018A%26A...616A..33S/abstract). The code uses Bayesian inference and provides publication quality plots, ascii parameter files readable by [*Topcat*](http://www.star.bris.ac.uk/~mbt/topcat/) and complete weighted posterior samples can be saved in the form of HDF5 files for direct analysis on the probability density functions of the derived stellar parameters.

## Installation
### From source (Only supported variant currently)
If you want to download the directory to a local path:
```bash
git clone https://github.com/StephanStock/SPOG.git
```

### Requirements
*SPOG* is written in Python 3 and should run with a standard [Anaconda distribution](https://www.anaconda.com/distribution/). The requirements are:
* astropy
* numpy
* matplotlib
* SciPy
* yaml
* pandas
* corner
* h5py


### Stellar Evolutionary Tracks
The code uses the stellar models based on the PAdova and TRieste Stellar Evolution  Code [Bressan et al. (2012)](https://ui.adsabs.harvard.edu/abs/2012MNRAS.427..127B/abstract) available under [this link](https://people.sissa.it/~sbressan/parsec.html).

However, the models require a particular preparation and certain modifications which are explained in more detail in Chapter 3 of [Stock et al. (2018)](https://ui.adsabs.harvard.edu/abs/2018A%26A...616A..33S/abstract) or further down in this readme. The prepared models can be downloaded from this [Link](https://heibox.uni-heidelberg.de/d/253b8d99e1324fa2b4f5/). Currently only sub-giant, red-giant and horizontal branch models are included. New prepared models will be available very soon allowing to also fit the main-sequence and pre-main sequence.


### Testing the installation
To test the functionality of *SPOG* after installing the program and the models go to the directory in which sp_main.py has been copied and run:
```bash
python sp_main.py param.yaml
```
Proceed through the program.

"param.yaml" is the path to a file in the "yaml"-format that contains the parameters with which the script is to be executed. All possible options are listed and explained in the example file:

```yaml
# Input file to define the parameters for the script.

#Essential parameters
object_name : 'Test_1p95msun'  # Name of the object, must be a python string, outputfiles will begin witht this name

save_path : '/home/sstock/Seafile/GJ251/stellar_param/' # Path were the results will be saved

model_path : '/media/sstock/Seagate Expansion Drive/Project_SPOG/Models_uncompressed_BV_new.h5'  # Path to the hdf5 file of the evolutionary models

photometric_band_A: 'V_johnson' #photometric band that defines the ordinate, currently supported kewords are V_johnson, B_johnson, I_johnson, J_2mass, H_2mass, Ks_2mass, G_gaia, G_BP_gaia, G_RP_gaia

photometric_band_B: 'B_johnson' #secondary photometric band, currently supported kewords are V_johnson, B_johnson, I_johnson, J_2mass, H_2mass, Ks_2mass, G_gaia, G_BP_gaia, G_RP_gaia

reverse: False   # if True: color equals photometric_band_A-photometric_band_B ; if False color equals photometric_band_B-photometric_band_A

mag_star : [7.18264, 0.010] # [magnitude of star, uncertainty] magnitude of star in photometric_band_A [mag]

color_star : [0.8788, 0.013] #[color of star, uncertainty] color in star [mag]

met_star : [0.05979, 0.01] #[metallicity, uncertainty] in Fe/H [dex]

par_star : [8.2719, 0.02] #[parallax, uncertainty] in [mas]

evolutionary_stage_prior : ['RGB', 'HB'] #[Include here a list of evolutionary stages that should be considered, you can mix any combination of subgiant or red-giant phase (RGB), red clump or horizontal branch phase (HB), main-sequence phase (MS), premain-sequence phase (PMS). ]


#Optional Parameters

par_err_dom : False #If the parallax error is assumed to dominate the uncertainty set to TRUE, if False slight BIAS is introduced since magnitude uncertainty is non-linear with respect to ABL

use_extinction : True  # Use extinction estimate

A_lambda : [0.21,0.1]      # [extinction A in photometric_band_A, uncertainty] in mag

E_color : [0.06,0.02]          # [reddening in color, uncertainty] in mag

parameterization: 'default'      # parametrization of derived paramers: default, default2 (as default but includes phase and mass loss), log, linear

mode: 'new'  # 'new': final parameters and uncertainties are sample quartiles [0.16, 0.5, 0.84] of weighted posterior (recommended and robust way as it is independent of binning!). 'classic': final parameters and uncertainties are weighted posterior modes (derived by spline interpolation) and uncertainties are calculated as in Eq. 11&12 of Stock et al., A&A 616, A33 (2018). This method may depend slightly on binning. This mode is only available in default parametrization.

smooth: 1.0 #Default should be 1. Posteriors are smoothed using a gaussian filter. This parameters sets the standard deviation for the Gaussian kernel. In classic mode it may effect the derived posterior modes.

model_sampling: 1           # Default is 1, higher integer value to use sparser sampling on the models by using every nth model in terms of metallicity, improving model load times but possibly affecting parameter precision and accuracy. Not recommended for final calculations

plot_corner: True #if True saves a corner plot

return_ascii: True # if True saves an ascii file (readable by Topcat)

plot_posterior: True #if True saves a posterior plot (without cornerplot)

posterior_bins: 20 #send the bin argument to astropy.stats.histogram (can be float or sting), if set to -1 the optimal number of bins will be calculated by using a self-implemented algorithm provided in Hogg (2008) (arkiv:0807.4820v1). Important note: The number of bins does not affect the results but only influences the visualised posterior plots!

posterior_fig_kwargs: {} #kwargs to sent to matplotlib.pyplot.subplot to influence the figure style, set empty {} for default

posterior_plot_kwargs: {} #kwargs to sent to matplotlib.pyplot.plot() to influence the plotting style, set empty {} for default

save_posterior: True    #if True saves all posterior results into a single hdf5 file, use pandas.read_hdf('filename.h5', 'key') where key is RGB or HB to reload the posterior into a Pandas dataframe

```

##A short description of applied changes to the evolutionary model by [Bressan et al. (2012)](https://ui.adsabs.harvard.edu/abs/2012MNRAS.427..127B/abstract)

Here you find some information about how the evolutionary models have been prepared and which parameter ranges are included in this routine.

1. A mass loss of $\eta=0.2$ was applied to the RGB models.  
2. The ages of the ZAHB models for stars experiencing a helium flash were corrected according to the applied mass loss in the RGB.
3. The evolutionary models have been interpolated to a much finder grid, each in mass and metallicity. For this interpolation the phase value apparent in each evolutionary track has been used. For details see [Stock et al. (2018)](https://ui.adsabs.harvard.edu/abs/2018A%26A...616A..33S/abstract). The following information is from [here](https://people.sissa.it/~sbressan/CAF09_V1.2S_M36_LT/readme.txt):

    1 PMS_BEG  Track begins here (Pre main sequence)\
    2 PMS_MIN  \
    3 PMS_END    PMS is near to end\
    4 NEAR_ZAM   This point is very near the ZAMS\
    5 MS_BEG     H burning fully active\
    6 POINT_B    Almost end of the H burning. Small contraction phase begins here for interm. & massive stars \  
    7 POINT_C    Small contraction ends here and star move toward RG\
    8 RG_BASE    RG base\
    9 RG_BMP1   RGB bump in Low Mass Stars (marked also for other masses)\
    10 RG_BMP2   RGB bump end in Low Mass Stars (marked also for other masses)\
    11 RG_TIP   Helium Flash or beginning of HELIUM Burning in intermediate and massive stars\
    12 Loop_A   Base of Red He burning, before possible loop\
    13 Loop_B   Bluest point of the loop (He burning\
    14 Loop_C   central He = 0  almost\
    15 TPAGB begins or c_burning begins (massive stars, generally if LC > 2.9)


    Note that the ZAHB of stars undergoing a helium flash have a different definition starting from 1 to 5.  

4. The following transformation from Z to [Fe/H] has been applied:

    Z_{Sun}=0.0152\
    Y_{Sun}=0.2485+1.78*Z_{Sun} ;Padova relation\
    X_{Sun}=1.-Y_{Sun}-Z_{Sun}\
    Y=0.2485+1.78*Z\
    X=1.-Y-Z\
    [Fe/H]=LOG10(Xsun/Zsun)+LOG10(Z/X)

5. Bolometric corrections were applied so that the evolution models contain absolute magnitudes of different photometric bands. The mass loss along the RGB was taken into account when calculating the absolute brightness in different photometric bands. The bands supported are V, B, I calculated based on [Worthey & Lee (2011)](https://ui.adsabs.harvard.edu/abs/2011ApJS..193....1W/abstract), J, H, Ks based on the 2MASS filters, and G, G_BP, G_RP based on the GAIA DR2 filters. The source for the latter correction is dustyAGB07, available at this [link](http://stev.oapd.inaf.it/cmd_3.6/photsys.html).

6. A cut on the evolutionary models was applied, which is based on the age of the universe. Only models with ages less than 13 Gyr are considered.

7. The metallicity range of the models available for this routine is from Z=0.0005 to Z=0.06 ([Fe/H] from ~ -1.51 to ~ +0.68).

8. The range of masses for this routine is from about 0.1 M_{Sun} to 12 M_{Sun}. Higher mass models are currently not supported due to a difference of how the phase value is defined for these models compared to models with smaller mass.



### Citation
If you are using *SPOG* for your research *please cite our paper* **[Stock et al. (2018)](https://ui.adsabs.harvard.edu/abs/2018A%26A...616A..33S/abstract)** and link this github repository.
The BibTeX entry is:
```
@ARTICLE{2018A&A...616A..33S,
       author = {{Stock}, Stephan and {Reffert}, Sabine and {Quirrenbach}, Andreas},
        title = "{Precise radial velocities of giant stars. X. Bayesian stellar parameters and evolutionary stages for 372 giant stars from the Lick planet search}",
      journal = {\aap},
     keywords = {stars: fundamental parameters, stars: late-type, stars: evolution, Hertzsprung-Russell and C-M diagrams, planetary systems, methods: statistical, Astrophysics - Solar and Stellar Astrophysics, Astrophysics - Earth and Planetary Astrophysics},
         year = 2018,
        month = aug,
       volume = {616},
          eid = {A33},
        pages = {A33},
          doi = {10.1051/0004-6361/201833111},
archivePrefix = {arXiv},
       eprint = {1805.04094},
 primaryClass = {astro-ph.SR},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2018A&A...616A..33S},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

```
Please make sure to also acknowledge the python packages required for *SPOG* as well as the stellar evolutionary models by [Bressan et al. (2012)](https://ui.adsabs.harvard.edu/abs/2012MNRAS.427..127B/abstract) and the source for the bolometric corrections.
