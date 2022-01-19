from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='SPOG',
    version='1.0',
    description='SPOG is a Python script for uncomplicated determination of'
    'stellar parameters based on the method used and outlined'
    'in Stock et al. (2018).',
    license="MIT",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author='Stephan Stock @ ZAH, Landessternwarte Heidelberg',
    author_email='sstock@lsw.uni-heidelberg.de',
    url="https://github.com/StephanStock/SPOG",
    packages=setuptools.find_packages(),
    install_requires=['astropy', 'numpy', 'matplotlib',
                      'scipy', 'pandas', 'corner', 'h5py',
                      'pyyaml', 'tqdm', 'wget'],
    python_requires='>=3.6',
    entry_points={'console_scripts': ['SPOG=SPOG:main']},
)
