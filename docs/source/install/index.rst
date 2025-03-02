.. _install:

=========================
Installation Instructions
=========================

dorado is compatible with Python 3.11+. There are only 4 dependencies: `numpy <https://numpy.org/install/>`_, `matplotlib <https://matplotlib.org/3.2.2/users/installing.html>`_, `scipy <https://www.scipy.org/install.html>`_, `future <https://python-future.org/>`_, and `tqdm <https://pypi.org/project/tqdm/>`_.

Installation via `pip`
----------------------

To `pip`-install this package, first ensure that you have the dependencies listed above installed, and then use the following command:
::

    $ pip install pydorado

Installation via `conda`
------------------------

This package is available via `conda-forge <https://anaconda.org/conda-forge>`_ and can be installed using `conda` (or `mamba`) using the command:
::

    $ conda install -c conda-forge pydorado


Installation from source
------------------------
1. Clone (or download) the repository
::

   $ git clone https://github.com/passaH2O/dorado

2. From the cloned (or extracted) folder, run the following in the command line:
::

   $ pip install .

to install the dorado package.

3. You can test your installation using `pytest`. To do this, first install pytest by running the following from the command line:
::

   $ pip install pytest

Then to run the unit tests, type the following from the cloned (or extracted) folder:
::

   $ pytest


Editable installation from source
---------------------------------
If you'd prefer an "editable" install (meaning that any modifications you make to the code will be used when you import and run scripts), run the following in the command line after cloning the repository (instead of following the above instructions):
::

   $ pip install -e .

The unit tests can be run (after installing `pytest` as indicated above) by typing:
::

   $ pytest

If you'd like your installation to include the dependencies needed for testing and locally building the documentation, you can install the package with the following command:
::

   $ pip install -e .[dev]

This should give you a "developers" installation of the package.
