---
title: 'dorado: A Python package for simulating passive particle transport in shallow-water flows'
tags:
  - Python
  - Lagrangian transport
  - random walk
  - hydrodynamics
authors:
  - name: Jayaram Hariharan
    orcid: 0000-0002-1343-193X
    affiliation: 1
  - name: Kyle Wright
    orcid: 0000-0001-5142-1786
    affiliation: 1
  - name: Paola Passalacqua
    orcid: 0000-0002-4763-7231
    affiliation: 1
affiliations:
 - name: Department of Civil, Architectural and Environmental Engineering, The University of Texas at Austin
   index: 1
date: 8 August 2020
bibliography: paper.bib
---

# Summary

Hydrodynamic simulations of flow through landscapes allow scientists to answer questions related to the transport of water, nutrients, pollutants, biota, and sediment through waterways [@Duan2006; @Rynne2016; @Wild-Allen2016; @Lauzon2018; @Czuba2019]. In geophysical systems, such as rivers, estuaries, and deltas, hydrodynamic models typically solve the depth-integrated ``shallow water" equations in an Eulerian reference frame, which is concerned with fluxes through a given region of space -- examples of these solvers include ANUGA [@anuga], Delft3D [@delft3d], Frehd [@Hodges2014] and others. However, the spatial and temporal characteristics of the movement of material through a landscape are often better understood using a Lagrangian reference frame [@Doyle2009], which follows the movement of individual objects or parcels. In this paper, we present an open-source Python package, *dorado*, which provides a transparent and accessible method for researchers to simulate passive Lagrangian particle transport on top of Eulerian hydrodynamic solutions. This mixed Eulerian-Lagrangian methodology adapts the routing functionality from the popular numerical model DeltaRCM [@Liang2015a; @Liang2015] for use with the outputs of any shallow-water hydrodynamic solver.

# Statement of Need

Existing software for tracking water parcels or particles in a Lagrangian framework typically contains complicated code structures [@freshkiss3d; @delft3d] with steep learning-curves, are proprietary [@fluent], or were developed to solve problems at very different spatial scales [@Yeung1989; @Dagestad2018] than typical riverine applications. This package fills a gap in the available methods for Lagrangian particle simulation by providing a flexible, open, and transparent framework meant for rapid application in conjunction with any landscape-scale 2D hydrodynamic model. In addition, *dorado* comes with built-in pre-processing, analysis, and plotting functionality, including methods to compute the exposure time distribution of particles to specific sub-regions of the domain, which is often of interest in ecological applications [@Kadlec2008].

# Background

*dorado* makes use of the weighted random walk framework [@Pearson1905] to model particle transport as a stochastic Markovian process, using a grid-centric approach that takes into account local water inertial components and surface slopes, in a manner modeled after the numerical model DeltaRCM [@Liang2015; @Liang2015a]. The random walk routing weights are sensitive to two user-specified parameters, $\gamma$ and $\theta$. The $\gamma$ parameter controls the proportional importance of water surface gradient and water velocity in determining the downstream routing direction $F^{*}$, according to:
\begin{equation}
    F^{*} = \gamma F_{sfc} + (1-\gamma) F_{int}
\end{equation}
in which $F_{sfc}$ and $F_{int}$ represent the surface gradient and inertial components, respectively. Modifying $\gamma$ indirectly controls the diffusivity of the travel path. The second weighting parameter, $\theta$, modifies the routing weights, $w_i$, for the random walk based on the local water depth value, $h_i$. The routing weights for each neighboring cell are calculated per @Liang2015a:

\begin{equation}
    w_i = \frac{h_i^\theta \max \left(0, F \cdot d_i \right)}{\Delta_i}
\end{equation}
wherein $F$ are the routing directions, $d_i$ represents the unit vector from the current cell to its neighbor, and $\Delta_i$ is the D8 cellular distance. Here, the local flow depth acts as a resistance term -- as $\theta$ gets larger, the local routing weights have an increasingly large dependence on the water depth of the neighboring cells. By altering the two parameters $\gamma$ and $\theta$, the user has control over the randomness with which the particles are routed. For a full description of this methodology, see @Liang2015, @Liang2015a and the *dorado* documentation.

Particle travel times are back-calculated from their locations after routing by accounting for the flow velocity, flow direction, and grid size. The time elapsed during each particle step is assumed to be proportional to the distance traveled in the direction of the mean flow, $d_{eq}$, which is equal to the step distance $\Delta_i$ projected onto a unit vector oriented in the direction of the velocity field, $\phi$, according to $d_{eq} = \Delta_i \cdot \cos(\phi)$. During a model iteration, each particle is allowed to progress a different length of time, depending on the velocity of the flow in the vicinity of each particle and the orientation of each step in relation to the mean flow direction. The travel time of each step is then modified according to a user-specified parameter, $D_c$, which acts as a dispersion coefficient.

*dorado* provides functions with which users can choose to sync up the particle evolution by either number of step iterations or by individual particle travel times. In addition, *dorado* includes several post-processing functions for computing and plotting the fraction of time particles spend ``exposed" to a region of interest, otherwise known as the exposure time distribution (or residence time distribution, in steady flows) [@Kadlec2008; @Benjamin2013; @Hiatt2018]. This travel time methodology has been tested against analytical solutions for a plug-flow reactor with dispersion [@Benjamin2013], as well as against exposure time distributions from prior models of real systems [@Hiatt2018], and performed well at reproducing the observed travel times in both systems. As with all Lagrangian methods involving finite samples of particles, we expect the travel time computations to struggle to reproduce the heavy-tailed behavior observed in some systems [@Zhang2012], and expect that *dorado* will perform most accurately in advection-dominated flows.

# Functionality and Ease of Use

*dorado* is capable of handling both steady and unsteady hydrodynamic simulations. While the core functionality of *dorado* assumes that flow variables lie on a Cartesian grid of uniform grid size, the package can still be applied to unstructured hydrodynamic models through the use of built-in interpolation functions for gridding data. Due to the nature of the weighting scheme, model flow-fields must represent depth-averaged solutions (as vertical movement of particles is not considered). *dorado* supports current Python 3.x versions as well as the older 2.7 version, can be run on all major operating systems, and is dependent upon only a select few common scientific packages. To help new users, over a dozen examples are provided, illustrating uses of the package with DeltaRCM and ANUGA model simulations conducted on real and synthetic datasets. The package is organized into four modules to help users and future contributors easily access and interact with the source code:

  - `lagrangian_walker`: Contains internal functions and methods associated with the weighted random walk.
  - `particle_track`: Lower-level functionality in which the `modelParams` and `Particles` classes are defined, and parameters for the domain and particle evolution are set.
  - `routines`: Higher-level functions for common applications, designed to simplify use of the package.
  - `parallel_routing`: Parallel functionality which enables the user to distribute their particle routing commands across multiple CPU cores using Python's built-in multiprocessing library.

# Acknowledgements

We thank Nelson Tull for testing and providing feedback on the package. Thanks also to David Mohrig, Kathleen Wilson, Hima Hassenruck-Gudipati, Teresa Jarriel, Eric Prokocki, John Swartz, Shazzadur Rahman, and Andrew Moodie for their input and feedback during development. We also thank Mickey Lanning for suggesting the package name.

This work was supported in part by NSF EAR-1719670, the NSF GRFP under grant DGE-1610403, and the NASA Delta-X project, which is funded by the Science Mission Directorate’s Earth Science Division through the Earth Venture Suborbital-3 Program NNH17ZDA001N-EVS3.

# References
