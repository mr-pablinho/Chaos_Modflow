# -*- coding: utf-8 -*-
"""
Created on Tue Dec 28 21:04:07 2021
2D plots
@author: PMR
"""

# %% Import libraries

# import external libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# fonts and colors
import matplotlib
font = {'size'   : 10}
matplotlib.rc('font', **font)
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family'] = "sans-serif"

# import inner modules and variables
from f_colors import cmap_pmr_4, cmap_pmr_5, cmap_pmr_6, cmap_pmr_div_5, cmap_pmr_div_6
from set_sim import t, nrow, ncol, wave_names, nodes, weights, J_distro, polynomial_expansion, delr, delc, hk, duration


# %% Plot: UQ in the groundwater heads 2D

# import heads
E_gwh = np.zeros((len(wave_names), len(t), nrow, ncol))
S_gwh = np.zeros((len(wave_names), len(t), nrow, ncol))

for wave in range(len(wave_names)):
    E_gwh[wave,:,:,:] = np.load('./outputs/uq_heads/E_heads__%s.npy' % (wave_names[wave]))
    S_gwh[wave,:,:,:] = np.load('./outputs/uq_heads/S_heads__%s.npy' % (wave_names[wave]))


titles_wave = ['Sine waves',
               'Trapezoid waves',
               'Triangle waves',
               'Composed waves']

# display of the subplots
rowfig = 2
colfig = 3

# choose stress periods
time_eval = [84,115]
# time_eval = np.arange(0,30,3)

# set min and max for the colorbar
vmin_e, vmax_e = -1, 1
vmin_s, vmax_s = 0, 1

for ti in time_eval:
    fig, ax = plt.subplots(nrows=rowfig, ncols=colfig, figsize=[7,7], 
                            gridspec_kw={'width_ratios':[1,1,1],'height_ratios': [1,1]},
                            constrained_layout=True)
    
    fig.suptitle('SP: %d' % (ti), fontsize=20, horizontalalignment='left')
    
    for wave in range(len(wave_names)):
        mat_plot = E_gwh[wave][ti][:, :61]
        
        # correction for "-inf" and "nan"
        mat_plot = np.where(~np.isnan(mat_plot), mat_plot, -10)
        mat_plot = np.where(np.isfinite(mat_plot), mat_plot, -10)
        
        ime = ax[0][wave].imshow(mat_plot, 
                                 vmin=vmin_e, vmax=vmax_e, 
                                 cmap='Spectral')
        ax[0][wave].set_title(titles_wave[wave])
    ax[0][0].set_ylabel(r'$\mu$')
     
    # position for the colorbar (mu)
    cbaxes_e = fig.add_axes([1.04, 0.50, 0.025, 0.405]) 
    cb = plt.colorbar(ime, cax=cbaxes_e)
    
    
    
    for wave in range(len(wave_names)):
        ims = ax[1][wave].imshow(S_gwh[wave][ti][:, :61], 
                                 vmin=vmin_s, vmax=vmax_s, 
                                 cmap='gnuplot2')
    
    ax[1][0].set_ylabel(r'$\sigma$')
    
    # position for the colorbar (sigma)
    cbaxes_s = fig.add_axes([1.04, 0.04, 0.025, 0.405]) 
    cb = plt.colorbar(ims, cax=cbaxes_s)
    
    # plt.savefig('./images/heads/2d_gwh_uq_%d.png' % (ti), dpi=300, bbox_inches='tight')
    # plt.close(fig) 