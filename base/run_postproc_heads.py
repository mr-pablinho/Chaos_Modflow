# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 16:03:38 2021
Last update on Thu Nov 25 15:56:21 2021
GROUNDWATER HEADS RESULTS
"""


# %% Import libraries

# import external libraries
import numpy as np
import chaospy as cp
import matplotlib.pyplot as plt

# import inner modules and variables
from set_sim import t, nrow, ncol, wave_names, nodes, weights, J_distro, polynomial_expansion

# fonts and colors
import matplotlib
font = {'size'   : 10}
matplotlib.rc('font', **font)
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family'] = "sans-serif"

# %% Import results

heads_all = np.zeros((int(len(wave_names)), int(len(nodes.T)), int(len(t)), int(nrow), int(ncol)))

# create list for heads: [wave][node][time,row,col]
for wave in range(len(wave_names)):
    for node in range(len(nodes.T)):
        doo = np.load('./outputs/%s/heads_%d.npy' % (wave_names[wave], node))
        heads_all[wave,node,:,:,:] = doo
        print('Gathering groundwater heads: wave=%s, node=%d' % (wave_names[wave], node))


# %% Setup results 

# define evaluation points
# obs_point = [[10,0], [10,10], [10,20], [10,30], [10,40], [10,50], [10,60], [10,70]]
obs_point = [[10,10]]
obs_point_data = [None] * len(obs_point)

space = np.arange(0, ncol, 1)


for wave in range(len(wave_names)):
    plt.figure('Time v. head - ' + wave_names[wave])
    for node in range(len(nodes.T)):
        plt.plot(t, 
                 heads_all[wave,node,:,10,5], 
                 alpha=0.5, linewidth=0.5)
    plt.xlim(300000, 400000)
    plt.ylim(3,7)
    
    plt.figure('Space v. head - ' + wave_names[wave])
    for node in range(len(nodes.T)):
        plt.plot(space, 
                 heads_all[wave,node,200,10,:], 
                 alpha=0.5, linewidth=0.5)
    # plt.xlim(300000, 400000)
    plt.ylim(-2,2)


# %% Compute uncertainty

E_gwh = np.zeros((len(t), nrow, ncol))
S_gwh = np.zeros((len(t), nrow, ncol))

period = np.arange(0, len(t), 1)
# period = np.linspace(0, 250-1, 83).astype(int)
# period = np.arange(30, 50, 1)

for wave in range(len(wave_names)):
    for row in range(nrow):
        print('Computing polynomial expansions (wave:%s, row:%d)' % (wave_names[wave], row+1))
        for col in range(ncol):
            cell_eval = heads_all[wave, :, :, row, col]
            f_approx = cp.fit_quadrature(polynomial_expansion, nodes, weights, cell_eval)
            E_gwh[:, row, col] = cp.E(f_approx, J_distro)
            S_gwh[:, row, col] = cp.Std(f_approx, J_distro)

                
# %% Save arrays

    np.save('./outputs/uq_heads/E_heads__%s' % (wave_names[wave]), E_gwh)
    np.save('./outputs/uq_heads/S_heads__%s' % (wave_names[wave]), S_gwh)


