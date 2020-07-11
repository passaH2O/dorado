.. _quickstart:

==========
Quickstart
==========

Quick guide to using the `particlerouting` package using the 2 built-in sample datasets.


Quick Demo 1
------------
In this demo, we show how particles can be routed along a simulated river delta, the example output is from the delta simulation model, `DeltaRCM <https://github.com/DeltaRCM/pyDeltaRCM_WMT>`_.

First we load the sample parameters.

.. doctest::

    >>> import particlerouting as pr
    >>> rcmparams = pr.example_data.define_params.make_rcm_params()

We can visualize the water depth for this scenario from the parameters.

.. doctest::

    >>> plt.imshow(rcmparams.depth)
    >>> plt.title('Water Depth')
    >>> plt.show()

.. plot:: quickstart/demo1_depth.py

Now let's route 50 particles for 50 iterations.

.. doctest::

    >>> pr.routines.steady_plots(rcmparams, 50, 'demo-1')

The `steady_plots()` function saves plots of each iteration of the particle movement to the subfolder 'demo-1/figs'. The final plot is shown below; the initial particle locations are shown as blue dots, and the final particle locations are red dots.

.. plot:: quickstart/demo1.py


Quick Demo 2
------------
In this demo, we show how particles can be placed on a gridded flow field generated using the `ANUGA <https://github.com/GeoscienceAustralia/anuga_core>`_ software package. Instead of specifying the number of iterations for the particles to walk, we will be specifying a target travel time for the particles to travel in this demo.

First we load the sample parameters and define our `Particle` object.

.. doctest::

   >>> import particlerouting as pr
   >>> anugaparams = pr.example_data.define_params.make_anuga_params()
   >>> particle = pt.Particle(anugaparams)

We can visualize the flow discharge components for this scenario from the loaded parameters.

.. doctest::

   >>> plt.subplot(1, 2, 1)
   >>> plt.imshow(anugaparams.qx)
   >>> plt.title('x-component of water discharge')
   >>> plt.subplot(1, 2, 2)
   >>> plt.imshow(anugaparams.qy)
   >>> plt.imshow(anugaparams.qy)
   >>> plt.title('y-component of water discharge')
   >>> plt.show()

Now let's route 50 particles with a target travel time of 2100 seconds. We will return the output of the run 

.. doctest::

   >>>

The `steady_plots()` function saves plots of each iteration
