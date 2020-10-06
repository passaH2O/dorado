.. _example12:

Example 12 - Unsteady Flow Fields
=================================

In this example we will revisit the `ANUGA <https://github.com/GeoscienceAustralia/anuga_core>`_ model domain we looked at in :ref:`example01`. This time, however, we will be routing our particles in an unsteady flow field, using the data contained in the `examples/example_data` subdirectory.

.. Note::
   This example must be run from the "examples" directory.

Full example script available :download:`here <../../../examples/unsteady_example.py>`.

First we have to load our modules and define the parameter items that are related to the cell size and the particle attributes (number, initial location).

.. doctest::

   >>> from dorado.routines import unsteady_plots

   >>> # define information that is not in the saved data
   >>> dx = 5.
   >>> Np_tracer = 50
   >>> seed_xloc = list(range(5, 16))
   >>> seed_yloc = list(range(48, 53))

Now that the basic parameters have been established, we will apply the `unsteady_plots()` function to perform particle routing over unsteady flow data. The target travel time per timestep is set as 75 seconds.

.. doctest::

   >>> walk_data = unsteady_plots(dx, Np_tracer, seed_xloc, seed_yloc,
   >>>                            26, 75., 'unsteady_data',
   >>>                            'csv', 'unsteady_output')

If you run the example, data and figures will be placed in a subdirectory called `unsteady_output`. Below, there is an animation of the particles moving on the varying flow field. We can see that when the flow is inverted, our 'tide' pushes the particles back towards the inlet channel.

.. image:: images/example12/unsteady_example.gif
