# -*- coding: utf-8 -*-
"""
Created on Sat Dec 11 19:53:01 2021
CREATE COLORMAPS FOR PLOTS
@author: PMR
"""

import numpy as np
from matplotlib.colors import LinearSegmentedColormap

# colormap - 6 colors ----------------------------------------------------------------------------------------------------------
pmr_map_6 = [[223, 164, 201], [190,  71, 144], [148,  52, 111], [ 90,  32,  68],  # purple (from light to dark)
             [ 98, 176, 249], [  8, 122, 228], [  6,  86, 161], [  3,  42,  79],  # blue (from light to dark)
             [151, 237, 211], [ 43, 218, 165], [ 30, 170, 128], [ 18, 104,  78],  # green (from light to dark)
             [254, 229, 154], [254, 209,  72], [254, 197,  28], [223, 168,   1],  # yellow (from light to dark)
             [255, 202, 133], [255, 170,  59], [255, 149,  10], [205, 100,  19],  # orange (from light to dark)
             [248, 124, 119], [246,  63,  56], [238,  19,  11], [175,  14,   8]]  # red (from light to dark)

pmr_map_6 =  np.array(pmr_map_6)/255

n_bins = [4*6]   # Discretizes the interpolation into bins
cmap_name = 'pmr'
for n_bin in n_bins:
    # Create the colormap
    cmap_pmr_6 = LinearSegmentedColormap.from_list(cmap_name, pmr_map_6, N=n_bin)

# colormap - 5 colors ----------------------------------------------------------------------------------------------------------
pmr_map_5 = [[223, 164, 201], [190,  71, 144], [148,  52, 111], [ 90,  32,  68],  # purple (from light to dark)
             [ 98, 176, 249], [  8, 122, 228], [  6,  86, 161], [  3,  42,  79],  # blue (from light to dark)
             [151, 237, 211], [ 43, 218, 165], [ 30, 170, 128], [ 18, 104,  78],  # green (from light to dark)
             [254, 229, 154], [254, 209,  72], [254, 197,  28], [223, 168,   1],  # yellow (from light to dark)
             # [255, 202, 133], [255, 170,  59], [255, 149,  10], [205, 100,  19],  # orange (from light to dark)
             [248, 124, 119], [246,  63,  56], [238,  19,  11], [175,  14,   8]]  # red (from light to dark)

pmr_map_5 =  np.array(pmr_map_5)/255

n_bins = [4*5]   # Discretizes the interpolation into bins
cmap_name = 'pmr'
for n_bin in n_bins:
    # Create the colormap
    cmap_pmr_5 = LinearSegmentedColormap.from_list(cmap_name, pmr_map_5, N=n_bin)


# colormap - 4 colors ----------------------------------------------------------------------------------------------------------
pmr_map_4 = [
    # [223, 164, 201], [190,  71, 144], [148,  52, 111], [ 90,  32,  68],  # purple (from light to dark)
             [ 98, 176, 249], [  8, 122, 228], [  6,  86, 161], [  3,  42,  79],  # blue (from light to dark)
              [151, 237, 211], [ 43, 218, 165], [ 30, 170, 128], [ 18, 104,  78],  # green (from light to dark)
             [254, 229, 154], [254, 209,  72], [254, 197,  28], [223, 168,   1],  # yellow (from light to dark)
             # [255, 202, 133], [255, 170,  59], [255, 149,  10], [205, 100,  19],  # orange (from light to dark)
             [248, 124, 119], [246,  63,  56], [238,  19,  11], [175,  14,   8]]  # red (from light to dark)

pmr_map_4 =  np.array(pmr_map_4)/255

n_bins = [4*4]   # Discretizes the interpolation into bins
cmap_name = 'pmr'
for n_bin in n_bins:
    # Create the colormap
    cmap_pmr_4 = LinearSegmentedColormap.from_list(cmap_name, pmr_map_4, N=n_bin)

# colormap - 6 colors divergent-------------------------------------------------------------------------------------------------
pmr_map_div_6 = [[223, 164, 201], [190,  71, 144], [148,  52, 111], [ 90,  32,  68],  # purple (from light to dark)
                 [ 98, 176, 249], [  8, 122, 228], [  6,  86, 161], [  3,  42,  79],  # blue (from light to dark)
                 [ 18, 104,  78], [ 30, 170, 128], [ 43, 218, 165], [151, 237, 211],  # green (from light to dark)
                 [254, 229, 154], [254, 209,  72], [254, 197,  28], [223, 168,   1],  # yellow (from light to dark)
                 [255, 202, 133], [255, 170,  59], [255, 149,  10], [205, 100,  19],  # orange (from light to dark)
                 [248, 124, 119], [246,  63,  56], [238,  19,  11], [175,  14,   8]]  # red (from light to dark)

pmr_map_div_6 =  np.array(pmr_map_div_6)/255

n_bins = [4*6]   # Discretizes the interpolation into bins
cmap_name = 'pmr'
for n_bin in n_bins:
    # Create the colormap
    cmap_pmr_div_6 = LinearSegmentedColormap.from_list(cmap_name, pmr_map_div_6, N=n_bin)
    
# colormap - 5 colors divergent-------------------------------------------------------------------------------------------------
pmr_map_div_5 = [[223, 164, 201], [190,  71, 144], [148,  52, 111], [ 90,  32,  68],  # purple (from light to dark)
                 [ 98, 176, 249], [  8, 122, 228], [  6,  86, 161], [  3,  42,  79],  # blue (from light to dark)
                 [ 18, 104,  78], [ 30, 170, 128], [ 43, 218, 165], [151, 237, 211],  # green (from light to dark)               
                 [254, 229, 154], [254, 209,  72], [254, 197,  28], [223, 168,   1],  # yellow (from light to dark)
                 [255, 202, 133], [255, 170,  59], [255, 149,  10], [205, 100,  19],  # orange (from light to dark)
                 [248, 124, 119], [246,  63,  56], [238,  19,  11], [175,  14,   8]]  # red (from light to dark)

pmr_map_div_5 =  np.array(pmr_map_div_5)/255

n_bins = [4*5]   # Discretizes the interpolation into bins
cmap_name = 'pmr'
for n_bin in n_bins:
    # Create the colormap
    cmap_pmr_div_5 = LinearSegmentedColormap.from_list(cmap_name, pmr_map_div_5, N=n_bin)