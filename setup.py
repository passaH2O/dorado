import os
from setuptools import setup, find_packages

setup(
    name = 'pydorado',
    version = '2.5.4',
    license = 'MIT',
    description = 'dorado - Lagrangian particle routing routine via weighted random walks',
    author = 'J. Hariharan, K. Wright, P. Passalacqua',
    author_email = 'jayaram.hariharan@utexas.edu',
    url = 'https://github.com/passaH2O/dorado',
    packages = find_packages(exclude=['*.tests']),
    package_data = {'' : ['*.txt', '*.npz']},
    long_description = 'See project webpage for details: https://github.com/passaH2O/dorado',
    classifiers = [
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
        'Programming Language :: Python :: 3.14'],
    install_requires = ['numpy','matplotlib','scipy',
                        'future','tqdm'],
)
