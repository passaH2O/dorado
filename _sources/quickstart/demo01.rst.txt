.. Quick Demo 1_:

Demo 1- Using the High-Level API
--------------------------------
In this demo, we show how particles can be routed along a simulated river delta, example output is from the delta simulation model, `DeltaRCM <https://github.com/DeltaRCM/pyDeltaRCM>`_.
To do this, we will be using one of the high-level API functions provided in :download:`routines.py <../../../dorado/routines.py>`.
A set of 50 particles will be seeded in the main channel of the synthetic delta.
Then using information about the water depths and flow rates in the delta, the particles will be allowed to travel downstream.
River deltas are distributary networks, and so the particles will not all travel along the same channels.
By using a tool like `dorado`, simulated particles can help us answer questions related to hydrological connectivity and nutrient residence times in deltaic systems.

First we load and create the sample particles object.

.. doctest::

    >>> import dorado
    >>> import matplotlib.pyplot as plt
    >>> from dorado.example_data import define_params as dp
    >>> rcmparticles = dp.make_rcm_particles()

We can visualize the water depth for this scenario from the `rcmparticles` object. If you'd like to download this portion of the demo as a standalone script, it is available :download:`here <../pyplots/quickstart/demo1_depth.py>`.

.. doctest::

    >>> plt.imshow(rcmparticles.depth)
    >>> plt.colorbar()
    >>> plt.title('Water Depth')
    >>> plt.show()

.. plot:: quickstart/demo1_depth.py

Now let's route 50 particles for 50 iterations.

.. doctest::

    >>> dorado.routines.steady_plots(rcmparticles, 50, 'demo-1')

The `steady_plots()` function saves plots of each iteration of the particle movement to the subfolder 'demo-1/figs'. The final plot is shown below; the initial particle locations are shown as blue dots, and the final particle locations are red dots. If you'd like to download this portion of the demo as a standalone script, it is available :download:`here <../pyplots/quickstart/demo1.py>`.

.. plot:: quickstart/demo1.py
