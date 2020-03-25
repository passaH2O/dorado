import os
from setuptools import setup

setup(
    name = 'particlerouting',
    version = '1.0',
    description = 'Lagrangian particle routing routine via weighted random walks',
    author = 'J. Hariharan, K. Wright, P. Passalacqua',
    author_email = 'jayaram.hariharan@utexas.edu',
    url = 'https://github.com/',
    packages = ['particlerouting'],
    long_description = 'See project webpage for details: https://github.com/',
    classifiers = ['Programming Language :: Python :: 2.7'],
    install_requires = ['numpy','matplotlib','scipy','netcdf4']
)
