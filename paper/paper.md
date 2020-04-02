---
title: 'particlerouting: A Python package for simulating passive particle transport'
tags:
  - Python
  - Lagrangian transport
  - random walk
  - hydrodynamics
authors:
  - name: Jayaram Hariharan
    orcid: ?
    affiliation: 1
  - name: Kyle Wright
    orcid: ?
    affiliation: 1
  - name: Paola Passalacqua
    orcid: ?
    affiliation: 1
affiliations:
 - name: Department of Civil, Architectural and Environmental Engineering, The University of Texas at Austin
   index: 1
date: 4 April 2020
bibliography: paper.bib
---

# Summary

Hydrodynamic simulations of flow through different landscapes are performed to answer a variety of questions including those related to the transport of nutrients, chemicals, fish, and sediment through waterways [add citations]. In this paper, we present an open-source Python package, 'particlerouting', which is designed to provide a simple and accessible method for researchers to simulate Lagrangian transport of individual particles in flow fields. The majority of hydrodynamics solvers such as Anuga [@anuga], Delft3D [@delft3d], Frehd [cite] and others, provide flow fields as output. With 'particlerouting', a method to place particles in these flow fields and observe their movement through the domain is provided.

Existing particle tracking software contains complicated code structures [@freshkiss3d; @delft3d], is proprietary [@fluent], or was developed to solve problems occurring at oceanic and atmospheric scales [@Dagestad2018]. With 'particletracking', a weighted random walk is employed to simulate the movement of passive particles in a flow field. In doing so, this package fills the gap in the available methods for Lagrangian particle simulation by providing a simple and transparent framework which can be used to simulate particle movement on the output from any hydrodynamic model. 

# Background

Particle transport in 'particlerouting' is simulated through the use of a weighted random walk framework [@Pearson1905]. In particular, methods used in the numerical model DeltaRCM [@Liang2015a; @Liang2015] are applied to weight the random walk by the local water slopes and the water discharge components. The weights for the random walk are sensitive to 2 parameters, $\gamma$ and $\theta$. The $\gamma$ parameter controls the proportional dependence of water surface slope and water discharge on the downstream direction. A $\gamma$ of 0 means that the water slope is ignored when determining the downstream flow direction, and a $\gamma$ of 1 means that the discharge values (the flow field) is ignored. The second weighting parameter, $\theta$ modifies the routing weight for the random walk based on the local water depth value. The routing weights for each neighboring cell are calculated per [@Liang2015a], using $w_i = \frac{\frac{1}{R_i} \text{max} \left(0, F \cdot d_i \right)}{\Delta_i}$. The resistance value, $R_i$, is computed using the local water depth, $h_i$ and the parameter $theta$: $R_i = \frac{1}{h_i^\theta}$. When $\theta$ is 0, then the local depth values are not impacting the weighting factor associated with each cell, as $\theta$ gets larger, the local routing weights have an increasingly large dependence on the water depth of the neighboring cells. By altering the two parameters $\gamma$ and $\theta$, the user has control over the randomness with which the particles are routed. 

# Acknowledgements

We acknowledge contributions from

# References
