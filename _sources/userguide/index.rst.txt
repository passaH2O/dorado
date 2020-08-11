.. _userguide:

==========
User Guide
==========

Overview
--------

The basic workflow when using `dorado` can be summarized in two steps:

1. Initialize parameters using the `particle_track.params` class and provide information about the domain, where and how many particles are going to be seeded, and adjust any parameters as desired.

2. Either use one of the high-level API procedures defined in `routines.py`, or define a particle using the `particle_track.Particle` and simulation particle evolution using the `particle_track.run_iteration` method.

Here we will describe some of the functionality built into the high-level API as well as the options available at the lower-level.


Defining the `params`
---------------------

Defining the parameters is require prior to doing any particle routing. Initiation of the parameter class can be achieved by typing the following:

.. doctest::

   >>> import dorado.particle_track as pt
   >>> params = pt.params()

After establishing the parameter class, information about the domain must be provided, as well as information about where the particles are going to be seeded and how many particles should be used. As `dorado` is a generic package designed to be used with a variety of input datasets or model outputs, we have tried to make the required input parameters as flexible as possible, however some values **must** be provided.

Required Parameters
^^^^^^^^^^^^^^^^^^^

A brief description of the required parameters is listed below, for a comprehensive list, please refer to the API: :ref:`particletrack`.

- x and y coordinates for the particles to be seeded in (`seed_xloc` and `seed_yloc`)
- the number of particles to simulate (`Np_tracer`)
- the length along a square cell face (`dx`)
- 2 of the following 3 arrays describing the domain:

   - the water depth (`depth`)
   - the water stage (`stage`)
   - the topography (`topography`)

- either the x and y components of water velocity (`u` and `v`) or the x and y components of the water discharge (`qx` and `qy`)

Optional Parameters
^^^^^^^^^^^^^^^^^^^

For a list of the optional parameters, please see :ref:`particletrack`, these parameters include values related to the random walk weighting scheme, travel time calculations, and even the distance and direction assumptions related to the grid.

The High-Level API
------------------
High-level functionality is provided in the `routines.py` script. Many of the example scripts take advantage of these methods and functions to route particles and make plots with minimal lines of code. A non-exhaustive list of the provided functionality is below.

**Functions to route particles:**

* Particle movement given a steady flow field with saving of the data and images after each particle iteration
* Particle movement in an unsteady flow field with automated data and plot saving
* Particle movement in a steady flow field with visualization of individual particle travel times

**Functions to plot/interpret output:**

* Query particle locations and travel times at a given iteration
* Query particle locations at a given travel time
* Plot the particle exposure time distributions
* Animate the output images of particle locations
* Plot the travel paths specified particles have taken
* Plot the particle positions for a specified iteration or travel time

For additional detail, either view the script itself (:download:`routines.py <../../../dorado/routines.py>`), or refer to the API: :ref:`routines`.

The Lower-Level API
-------------------
Lower-level functionality is provided in the `particle_track.py` script. At this level particle/domain parameters are assigned and more direct iterations with the particle routing functions are possible. A short list of the functionality present in this module is below.

* Parameter and domain information definition in the `params` class
* Particle definition and the ability to run a particle iteration in the `Particle` class
* Function for calculating exposure time (and residence time) of particles in a defined region of interest
* Assorted functions for transforming from real coordinate systems to the raster domain used for the particle routing

For additional detail, either view the script itself (:download:`particle_track.py <../../../dorado/particle_track.py>`), or refer to the API: :ref:`particletrack`.
