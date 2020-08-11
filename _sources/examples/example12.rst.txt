.. _example12:

Example 12 - Unsteady Flow Fields
=================================

In this example we will revisit the `ANUGA <https://github.com/GeoscienceAustralia/anuga_core>`_ model domain we looked at in :ref:`example01`. This time, however, we will be routing our particles in an unsteady flow field, using the data contained in the `examples/example_data` subdirectory.

Full example script available :download:`here <../../../examples/unsteady_example.py>`.

First we have to load our modules and define the parameter items that are related to the cell size and the particle attributes (number, initial location).

.. doctest::

   >>> from dorado.routines import unsteady_plots
   >>> import dorado.particle_track as pt

   >>> # initialize a parameters object
   >>> params = pt.params()

   >>> # give params information not contained in the grid
   >>> params.dx = 5.
   >>> params.Np_tracer = 50
   >>> params.seed_xloc = list(range(5, 16))
   >>> params.seed_yloc = list(range(48, 53))

Now that the basic parameters have been established, we will apply the `unsteady_plots()` function to route our particles over unsteady flow data. The target travel time per timestep is set as 75 seconds.

.. doctest::

   >>> walk_data = unsteady_plots(params, 26, 75., 'unsteady_data',
   >>>                            'csv', 'unsteady_output')

If you run the example, data and figures will be placed in a subdirectory called `unsteady_output`. Below, there is an animation of the particles moving on the varying flow field. We can see that when the flow is inverted, our 'tide' pushes the particles back towards the inlet channel.

.. image:: images/example12/unsteady_example.gif
