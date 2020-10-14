.. _example04:

Example 4 - Toggling Steepest Descent
=====================================

In this example we will show you how to use the 'steepest_descent' toggle to create a deterministic simulation of the particle routing. When the 'steepest_descent' option is enabled, the random walk functionality is turned off and instead the particles move deterministically from cell to cell.

Full example script available :download:`here <../../../examples/steepest_descent_deltarcm.py>`.

After initializing the parameters and particles in the same way as in :ref:`example02`, we need to change one of the default settings for the parameters class to turn on 'steepest_descent'.

.. doctest::

   >>> params.steepest_descent = True

From here, we can continue as we did in :ref:`example02`.

.. doctest::

   >>> walk_data = steady_plots(particles, 50, 'steepest_descent_example')

.. image:: images/example04/steepest_descent.gif
