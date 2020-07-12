.. _userguide:

==========
User Guide
==========

Overview
--------

The basic workflow when using `particlerouting`, can be summarized in two steps.

1. Initialize parameters using the `particle_track.params` class and provide information about the domain, where and how many particles are going to be seeded, and adjust any parameters as desired.

2. Either use one of the high-level API procedures defined in `routines.py`, or define a particle using the `particle_track.Particle` and simulation particle evolution using either the `particle_track.run_iteration` method.

Here we will describe some of the functionality built into the high-level API as well as the options available at the lower-level.


Defining the `params`
---------------------

Defining the parameters is require prior to doing any particle routing. Initiation of the parameter class can be achieved by typing the following:

.. doctest::

   >>> import particlerouting.particle_track as pt
   >>> params = pt.params()

After establishing the parameter class, information about the domain must be provided, as well as information about where the particles are going to be seeded and how many particles should be used. As `particlerouting` is a generic package designed to be used with a variety of input datasets or model outputs, we have tried to make the required input parameters as flexible as possible however some values **must** be provided.

Required Parameters
^^^^^^^^^^^^^^^^^^^

A brief description of the required parameters is listed below, for a comprehensive list, please refer to the :ref:`apiref`.

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

For a list of the optional parameters, please see the :ref:`apiref`, these parameters include values related to the random walk weighting scheme, travel time calculations, and even the distance and direction assumptions related to the grid.

The High-Level API
------------------
Brief descriptions of the high-level functionality provided in the `routines.py` script will be provided here. For additional detail, either view the script itself (:download:`routines.py <../../../particlerouting/routines.py>`), or refer to the :ref:`apiref`.

The Lower-Level API
-------------------
Here there will be some discussion about the lower-level API and how it might be accessed or modified to make case-specific changes. 
