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

However, the models require a particular preparation and certain modifications which are explained in more detail in Chapter 3 of [Stock et al. (2018)](https://ui.adsabs.harvard.edu/abs/2018A%26A...616A..33S/abstract). The prepared models can be downloaded from this [Link](https://heibox.uni-heidelberg.de/d/253b8d99e1324fa2b4f5/). Currently only sub-giant, red-giant and horizontal branch models are included. New prepared models will be available soon allowing to also fit the main-sequence and pre-main sequence. 

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



## Citation
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
Please make sure to also acknowledge the python packages required for *SPOG* as well as the stellar evolutionary models by [Bressan et al. (2012)](https://ui.adsabs.harvard.edu/abs/2012MNRAS.427..127B/abstract).
