---
title: 'SPOG+: A Python script to determine Bayesian stellar parameters of stars'
tags:
  - Python
  - astronomy
  - stars: fundamental parameters
  - Hertzsprung-Russell and C-M diagrams
  - methods: statistical
authors:
  - name: Stephan Stock
    orcid: 0000-0002-1166-9338
    affiliation: "1" # (Multiple affiliations must be quoted)
affiliations:
 - name: Landessternwarte, Zentrum für Astronomie der Universität Heidelberg, Königstuhl 12, 69117 Heidelberg, Germany
   index: 1
date: 02 February 2022
bibliography: paper_spog.bib
---

# Summary
The mass and radius of a star are among its most important characteristics. However, the determination of these parameters for individual stars is extremely difficult, since a direct measurement is not possible in most cases. Nevertheless, the exact and precise determination of these parameters is enormously important for making statements about the formation and distribution of exoplanets, for example. The methodology underlying this code is based on the knowledge that stellar masses correlate with the brightness and temperature of the star. For this purpose, theoretical models, so-called evolutionary tracks, exist which correlate the time evolution of stars of a certain mass and element abundance, the so-called metallicity, with its spectral type and luminosity. These intrinsic parameters of the star can also be transformed into observable quantities, in particular into the absolute brightness in a specific wavelength range and the color of the star, by means of bolometric corrections and a distance determination of the star, e.g. in the form of a trigonometric parallax measurement. The difficulty in using the theoretical models is that many models of different evolutionary stages may be degenerate, leading to ambiguous solutions for a given set of observed stellar parameters.


``SPOG+``  is a Python script that requires only four measurable quantities of the star, the apparent magnitude in two different photometric bands, the distance in terms of a trigonometric parallax, and the measured metallicity. With this information ``SPOG+`` will determine the probability distribution of six essential stellar parameters: mass, radius, age, surface gravity, effective temperature, and luminosity. ``SPOG+`` uses a Bayesian approach by determining the maximum likelihood of a measurement in a predefined model grid of evolutionary tracks by @Bressan2012 and multiplies the derived likelihoods with a prior probability derived within the code based on physical relevant properties. In particular, the last point is essential in better distinguishing strongly degenerate models. As shown in @Stock2018 and independently confirmed in @Malla2020, the masses determined by SPOG in @Stock2018 for post-main sequence stars are in very good agreement with direct mass determinations of these stars using asteroseismology.

``SPOG+``  is the pythonic evolution of the unpublished IDL code SPOG used in @Stock2018 and is made freely available to all astronomers as a Python version to greatly simplify the determination of essential stellar parameters. ``SPOG+``, however, allows, to not only derive masses of post-main sequence stars as in @Stock2018, but also for main-sequence and pre-main sequence stars, provided one downloads the predefined and specifically for SPOG+ formatted models offered with this software package. Another feature of ``SPOG+``  is to determine the probability of a certain evolutionary stage of the star and to obtain the parameters based only on this evolutionary stage. If desired, it is even possible to obtain a probability distribution of the precisely defined evolutionary state in form of an integer value given in the evolutionary tacks.

``SPOG+`` is specifically written for Python 3 and should run with a standard Anaconda distribution with a few additional packages. The script is executed by passing a YAML-parameter file with important information about the star and the desired output. As a result one obtains the probability distributions, corner plots, ASCII output files and even the entire probability density functions in the form of HDF5 files of the six above mentioned stellar parameters.

``SPOG+``  is developed to be used by astronomers that are in need of accurate and precise stallar parameters, most importanly stellar masses of individual field stars that do not allow or are not feasible for a direct measurement of their mass. The stellar parameters derived by the IDL predecessor SPOG have already been used in several studies, e.g., @Stock2018, @Trifonov2019, @Quirrenbach2019, and @Tala2020 to name a few. ``SPOG+`` will allow other independent groups to make use of this code for a large number of potential targets, e.g., observed by the GAIA satellite [@Gaia], and therefore may potentially be used in many further publications. ``SPOG+`` is designed to be user-friendly and allows to produce high quality publication plots which will further encourage the usage of this script by other researchers. To our knowledge there is no other code freely available that allows to derive precise stellar masses in such a straight-forward way by using only two different measurements of the apparent magnitudes, a metallicity estimate, and a parallax measurement of the star.

# Acknowledgements

We acknowledge support from Sabine Reffert during the genesis of this project.  We would also like to thank our two referees. This code makes use of Astropy [@Astropy_comm] a community-developed core Python package for Astronomy [@Astropy2018], NumPy [@oliphant2006guide], tqdm [@tqdm], SciPy [@Scipy2019], matplotlib [@Hunter2007], pandas [@pandas], h5py [@h5py], corner [@corner] and pytable [@pytable].

# References
