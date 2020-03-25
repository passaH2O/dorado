.. _examples:

========
Examples
========

Examples contained in the ParticleRouting/examples folder include the following:

steady_anuga_particles
----------------------
add link to script

Example of particle routing on a simple flow field generated using ANUGA Hydro.

.. image:: images/steady_anuga.gif
    :width: 800px
    :align: center

steady_deltarcm_particles
-------------------------
add link to script

Example of particle routing on a numerically generated river delta created using DeltaRCM.

.. image:: images/steady_deltarcm.gif
    :width: 800px
    :align: center

animate_deltarcm_particles
--------------------------
add link to script

An example script showing how to use the animate_plots routine to generate videos of the particle movement like those shown on this page.

steepest_descent_deltarcm
-------------------------
add link to script

An example script showing how the randomness can be turned off and a 'steepest descent' routing of the particles can be done by only allowing particles to travel to the neighboring cell with the greatest weight. 

.. image:: images/steepest_descent.gif
    :width: 800px
    :align: center

timing_anuga_particles
----------------------
add link to script

An example script which shows how the individual particle travel times can be visualized by coloring the particles themselves.

.. image:: images/timing_anuga.gif
    :width: 800px
    :align: center

true_random_walk
----------------
add link to script

An example script showing that by making the depth uniform, and setting the water discharge and both random walk parameters to 0, the particles move in a true random walk again.

.. image:: images/random_walk.gif
    :width: 800px
    :align: center

parallel_comparison
-------------------
add link to script

This script is meant to provide an example of the local parallelization. This functionality has been tested on Ubuntu 18.04, and may not work on Windows due to the ways in which Python handles multiple processes and the global lock. In this script, 2 cores are compared to running the same process in serial and an approximate 2x speed-up was achieved in this small case by using the parallel option.
