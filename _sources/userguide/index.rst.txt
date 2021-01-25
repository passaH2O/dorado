.. _userguide:

==========
User Guide
==========

Overview
--------

The basic workflow when using `dorado` can be summarized in three steps:

1. Initialize the parameters class (:obj:`dorado.particle_track.modelParams`) for parameters related to the grid (model domain). This includes information about flow fields, grid cell size, and tuning parameters for the random walk.

2. Initialize a :obj:`dorado.particle_track.Particles` class and generate a set of particles to be routed. A built in `particle_generator()` function is available to help with particle generation.

3. Either use one of the high-level API procedures defined in :download:`routines.py <../../../dorado/routines.py>`, or simulation particle evolution using the `run_iteration` method.

Here we will describe some of the functionality built into the high-level API as well as the options available at the lower-level.

Defining the `params`
---------------------

Defining the parameters is required prior to doing any particle routing. Initiation of the parameter class (:obj:`dorado.particle_track.modelParams`) can be achieved by typing the following:

.. doctest::

   >>> import dorado.particle_track as pt
   >>> params = pt.modelParams()

After establishing the parameter class, information about the domain must be provided. As `dorado` is a generic package designed to be used with a variety of input datasets or model outputs, we have tried to make the required input parameters as flexible as possible, however some values **must** be provided.

Required Parameters
^^^^^^^^^^^^^^^^^^^

A brief description of the required parameters is listed below, for a comprehensive list, please refer to the API documentation: :ref:`particletrack`.

- the length along a square cell face (`dx`)
- 2 of the following 3 arrays describing the domain:

   - the water depth (`depth`)
   - the water stage (`stage`)
   - the topography (`topography`)

- either the x and y components of water velocity (`u` and `v`) or the x and y components of the water discharge (`qx` and `qy`)

Optional Parameters
^^^^^^^^^^^^^^^^^^^

For a list of the optional parameters, please see :ref:`particletrack`, these parameters include values related to the random walk weighting scheme, travel time calculations, and even the distance and direction assumptions related to the grid.

Defining the `Particles`
------------------------

Defining a :obj:`dorado.particle_track.Particles` class is a key step in using `dorado` to perform particle routing. To define a set of particles, the model parameters must first be defined as described above. The `Particles` class is initialized using an instance of the model parameters. From there, particles can be generated and routed.

Particle Generation and Routing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Particles can be generated using the :obj:`dorado.particle_track.Particles.generate_particles` function, this allows for random and exact placement of particles

* Particles can be routed using either the High-Level or Low-Level API functionalities described below

The High-Level API
------------------
High-level functionality is provided in the :download:`routines.py <../../../dorado/routines.py>` script. Many of the example scripts take advantage of these methods and functions to route particles and make plots with minimal lines of code. A non-exhaustive list of the provided functionality is below.

Functions to route particles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Particle movement given a steady flow field with saving of the data and images after each particle iteration
* Particle movement in an unsteady flow field with automated data and plot saving
* Particle movement in a steady flow field with visualization of individual particle travel times

Functions to plot/interpret output
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Query particle locations and travel times at a given iteration
* Query particle locations at a given travel time
* Plot the particle exposure time distributions
* Animate the output images of particle locations
* Plot the travel paths specified particles have taken
* Plot the particle positions for a specified iteration or travel time

For additional detail, either view the script itself :download:`routines.py <../../../dorado/routines.py>`, or refer to the API documentation: :ref:`routines`.

The Lower-Level API
-------------------
Lower-level functionality is provided in the :download:`particle_track.py <../../../dorado/particle_track.py>` script. At this level particle/domain parameters are assigned and more direct iterations with the particle routing functions are possible. A short list of the functionality present in this module is below.

* Parameter and domain information definition in the `params` class
* Particle initialization and generation of particles
* The ability to run a particle iteration in the `Particle` class
* Function for calculating exposure time (and residence time) of particles in a defined region of interest
* Assorted functions for transforming from real coordinate systems to the raster domain used for the particle routing

For additional detail, either view the script itself :download:`particle_track.py <../../../dorado/particle_track.py>`, or refer to the API documentation: :ref:`particletrack`.
