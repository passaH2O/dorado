# -*- coding: utf-8 -*-
"""
dorado.spatial.py

Higher-level methods to quantify spatial exposure time.

Project Homepage: https://github.com/passaH2O/dorado
"""
import json
import numpy as np
import dorado.particle_track as pt
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import BoundaryNorm
from scipy.stats import gaussian_kde
from scipy.ndimage import gaussian_filter



def systemwide(walk_data, elevation, celltype, resolution_factor):
    """
    Compute localized exposure time by dividing the domain into fixed-size spatial chunks.
    Each chunk computes exposure time using only particles that started within the chunk,
    and using the celltype binary mask to define the region of interest.

    **Inputs** :
        walk_data : 'dict'
            Dictionary with 'xinds', 'yinds', and 'travel_times'.

        elevation : '2D float array'
            Terrain or bathymetry grid; used to define domain shape.

        celltype : '2D int array'
            Binary mask (1 = region of interest, 0 = ignore), same shape as elevation.

        resolution_factor : 'int'
            Number of cells per chunk edge (e.g., resolution = 5 -> 5x5 cells are in each chunk)

    **Outputs** :
        systemwide_exposure_times : 'list'
            Exposure time list for each chunk (NaN for empty chunks)
    """
    in_regions = celltype == 1
    rows, cols = elevation.shape

    # ROI
    regions = np.zeros_like(elevation, dtype="int")
    regions[in_regions] = 1

    n_chunks_y = rows // resolution_factor
    n_chunks_x = cols // resolution_factor

    exposure_time_data_chunks = []
    for r in range(n_chunks_y):
        for c in range(n_chunks_x):
            r0 = r * resolution_factor
            r1 = rows if r == n_chunks_y - 1 else (r + 1) * resolution_factor
            c0 = c * resolution_factor
            c1 = cols if c == n_chunks_x - 1 else (c + 1) * resolution_factor

            if np.any(in_regions[r0:r1, c0:c1]):
                chunk_walk = _find_common_indices_info(walk_data, (r0, r1), (c0, c1))
                exp_data = pt.exposure_time(chunk_walk, regions, verbose=False)
                exposure_time_data_chunks.append(exp_data)
            else:
                exposure_time_data_chunks.append(np.nan)

    return exposure_time_data_chunks





def _find_common_indices_info(walk_data, xind_range, yind_range):
    """ Extract particle trajectories whose initial positions fall within a specified grid chunk.
    
    This function is used to filter a subset of particles from `walk_data` based on their
    initial x and y index positions. Particle origin determines inclusion for
    chunk analysis. 

    **Inputs** :
        
        walk_data : 'dict'
            Dictionary with 'xinds', 'yinds', and 'travel_times'
        xind_range : 'tuple (int, int)'
            Minimum and maximum bounds of x-indices defining the chunk.
        yind_range : 'tuple (int, int)'
            Minimum and maximum bounds of y-indices defining the chunk.

    **Outputs** :
        
        common_indicies : 'dict'
            Dictionary containing the filtered subset of particle data.
            
    """
    common_indicies = {key: [] for key in walk_data}
    xinds = walk_data['xinds']
    yinds = walk_data['yinds']

    xind_min, xind_max = xind_range
    yind_min, yind_max = yind_range

    for idx in range(len(xinds)):
        if xinds[idx] and yinds[idx]:
            x0, y0 = xinds[idx][0], yinds[idx][0]
            if xind_min <= x0 <= xind_max and yind_min <= y0 <= yind_max:
                for key in walk_data:
                    common_indicies[key].append(walk_data[key][idx])

    return common_indicies




def localized(walk_data, elevation, celltype, resolution_factor):
    """
    Compute localized exposure time by dividing the domain into fixed-size spatial chunks.
    Each chunk computes exposure time using only particles that started within the chunk,
    and using the chunk itself as the region of interest.

    **Inputs** :
        walk_data : 'dict'
            Dictionary with 'xinds', 'yinds', and 'travel_times'

        elevation : '2D float array'
            Terrain or bathymetry grid defining the domain shape.

        celltype : '2D int array'
            Binary mask (1 = region of interest, 0 = ignore), same shape as elevation.

        resolution_factor : 'int'
            Number of cells per chunk edge (e.g., resolution = 5 -> 5x5 cells are in each chunk)

    **Outputs** :
        localized_exposure_times : 'list'
            Exposure time list for each chunk (NaN for empty chunks)
    """
    in_regions = celltype == 1
    rows, cols = elevation.shape
    n_chunks_y = rows // resolution_factor
    n_chunks_x = cols // resolution_factor

    # Binary region map (same as celltype)
    regions = np.zeros_like(elevation, dtype="int")
    regions[in_regions] = 1

    exposure_time_data_chunks = [] # Each chunk is ROI
    for r in range(n_chunks_y):
        for c in range(n_chunks_x):
            r0 = r * resolution_factor
            r1 = rows if r == n_chunks_y - 1 else (r + 1) * resolution_factor
            c0 = c * resolution_factor
            c1 = cols if c == n_chunks_x - 1 else (c + 1) * resolution_factor

            # Skip if chunk has no region of interest cells
            if np.any(in_regions[r0:r1, c0:c1]):
                region_chunk = np.zeros((rows, cols), dtype="int")
                region_chunk[r0:r1, c0:c1] = regions[r0:r1, c0:c1]

                chunk_walk = _find_common_indices_info(walk_data, (r0, r1), (c0, c1))
                exp_data = pt.exposure_time(chunk_walk, region_chunk)
                exposure_time_data_chunks.append(exp_data)
            else:
                exposure_time_data_chunks.append(np.nan)

    return exposure_time_data_chunks

        



def compute_thresholds(exposure_data,
                              elevation,
                              timedelta,
                              resolution_factor,
                              window_size=5,
                              threshold=3):
    """
    Compute smoothed exposure time percentiles (E50, E75, E90) on a chunked spatial domain.

    **Inputs** :

        exposure_data : 'list'
            Exposure time lists for each chunk.

        elevation_array : '2D float array'
            Grid of elevation or bathymetry used for shape.

        timedelta : 'float'
            Unit of time for time-axis of ETD plots, specified as time
            in seconds (e.g. 86400 = 1 day).

        resolution_factor : 'int'
            Number of chunks per axis.

        smoothing_window : 'int'
            Window size for local smoothing.

        smoothing_threshold : 'float'
            MAD-based filtering threshold to suppress outliers.

    **Outputs** :

        E50, E75, E90 : '2D arrays'
            Smoothed exposure time percentiles for each chunk.
    """

    def exposure_time_stats(exp_list):
        P50, P75, P90 = [], [], []
        for item in exp_list:
            if isinstance(item, list) and len(set(item)) > 0:
                filtered = [x for x in item if x != 0 and not np.isnan(x)]
                if len(set(filtered)) < 5:
                    P50.append(np.nan)
                    P75.append(np.nan)
                    P90.append(np.nan)
                    continue
                kde = gaussian_kde(filtered)
                x_grid = np.linspace(min(filtered), max(filtered), 1000)
                cdf = np.cumsum(kde(x_grid))
                cdf /= cdf[-1]
                def pval(p):
                    idx = np.searchsorted(cdf, p)
                    return x_grid[idx] if idx < len(x_grid) else np.nan
                P50.append(pval(0.50) / timedelta)
                P75.append(pval(0.75) / timedelta)
                P90.append(pval(0.90) / timedelta)
            else:
                P50.append(np.nan)
                P75.append(np.nan)
                P90.append(np.nan)
                
        shape = (elevation.shape[0] // resolution_factor, elevation.shape[1] // resolution_factor)
        return (np.array(P50).reshape(shape),
                np.array(P75).reshape(shape),
                np.array(P90).reshape(shape))

    def median_absolution_deviation_filter(data, window_size, threshold):
        smoothed = np.zeros_like(data, dtype='float')
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                if np.isnan(data[i, j]):
                    smoothed[i, j] = np.nan
                    continue
                i_min = max(i - window_size // 2, 0)
                i_max = min(i + window_size // 2 + 1, data.shape[0])
                j_min = max(j - window_size // 2, 0)
                j_max = min(j + window_size // 2 + 1, data.shape[1])
                window = data[i_min:i_max, j_min:j_max].flatten()
                window = window[~np.isnan(window)]
                if len(window) == 0:
                    smoothed[i, j] = np.nan
                    continue
                median = np.median(window)
                mad = np.median(np.abs(window - median))
                if mad == 0:
                    smoothed[i, j] = median
                    continue
                filtered = window[np.abs(window - median) <= (threshold * mad)]
                smoothed[i, j] = np.mean(filtered) if len(filtered) else np.nan
        return smoothed

    E50, E75, E90 = exposure_time_stats(exposure_data)
    E50 = median_absolution_deviation_filter(E50, window_size, threshold)
    E75 = median_absolution_deviation_filter(E75, window_size, threshold)
    E90 = median_absolution_deviation_filter(E90, window_size, threshold)

    return E50, E75, E90






def plot_spatial(exposure_maps,
                 elevation,
                 max_exposure,
                 cbar_levels,
                 cmap,
                 title):
    """
    Plot spatial exposure time maps for E50, E75, E90 percentiles using gridded chunks.

    **Inputs** :

        exposure_maps : 'list of 2D arrays'
            [E50, E75, E90] smoothed exposure time arrays. 

        elevation : '2D float array'
            Elevation or bathymetry grid used as shaded background.

        max_exposure : 'float'
            Upper limit of colorbar (e.g. 15 or 30 days).

        cbar_levels : 'int'
            Number of contour levels in the colorbar.

        cmap : 'colormap'
            Colormap for exposure values.

        title : 'str'
            Label for colorbar axis.

    """
    fig = plt.figure(figsize=(5.48, 1.4), dpi=600)
    gs = GridSpec(2, 4, figure=fig, hspace=0.06, wspace=0.06,
                  height_ratios=[0.53, 0.08], width_ratios=[0.05, 1, 1, 1])

    percentiles = [r'$E_{50}$', r'$E_{75}$', r'$E_{90}$']
    for j, label in enumerate(percentiles):
        fig.text(0.255 + 0.25 * j, 0.935, label, va='center', ha='center', fontsize=8)

    levels = np.linspace(0, max_exposure, cbar_levels)
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)
    

    for i, data in enumerate(exposure_maps):
        rows, cols = data.shape
        full_rows, full_cols = elevation.shape
      

        extent = [0, full_cols, full_rows, 0]  
        ax = fig.add_subplot(gs[0, i + 1])
        
        ax.imshow(elevation, cmap='bone', aspect='auto', alpha=0.75,
                  extent=extent, origin='upper')
    
    
        ax.imshow(data, cmap=cmap, norm=norm, aspect='auto',
                  interpolation='bilinear', extent=extent, origin='upper')
    

        ax.tick_params(axis='both', which='both', labelsize=8)
        ax.set_xticks([])
        ax.set_yticks([])


    cax = fig.add_axes([0.149, 0.13, 0.7525, 0.04])
    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),
                        cax=cax, orientation='horizontal')
    cbar.set_ticks(levels)
    cbar_labels = [str(int(t)) for t in levels]
    cbar_labels[-1] = f"{int(levels[-1])}+"
    cbar.ax.set_xticklabels(cbar_labels)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label(title, fontsize=8)
