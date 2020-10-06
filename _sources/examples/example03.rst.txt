.. _example03:

Example 3 - Using the Built-In Animation Function
=================================================

In this example we will demonstrate how to use the built-in animation functionality. Animation in `dorado` is done using the `matplotlib <https://matplotlib.org/3.2.2/users/installing.html>`_ library (to keep dependencies at a minimum), however we note that many other animation libraries exist for Python that may provide greater functionality.

To use the `animate_plots()` function, you must have installed the `ffmpeg <https://matplotlib.org/3.2.2/users/installing.html>`_ image writer.

If you have run the previous example (:ref:`example02`), then you can run the following and generate an animation similar to the one shown below.

Full example script available :download:`here <../../../examples/animate_deltarcm_particles.py>`.

.. doctest::

   >>> from dorado.routines import animate_plots
   >>> animate_plots(0, 50, 'steady_deltarcm_example')

.. image:: images/example02/steady_deltarcm.gif
