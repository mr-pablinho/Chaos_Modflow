# -*- coding: utf-8 -*-
"""
Created on Tue Dec 28 17:20:58 2021

@author: PMR
"""

# %% Import libraries

# import external libraries
import numpy as np
import chaospy as cp
import matplotlib.pyplot as plt

# fonts and colors
import matplotlib
font = {'size'   : 10}
matplotlib.rc('font', **font)
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family'] = "sans-serif"

# import inner modules and variables
import f_flow_field_eval as fpf
import f_wave_functions as fw
from set_sim import t, nrow, ncol, wave_names, nodes, weights, J_distro, polynomial_expansion, delr, delc, hk, duration, f


# %% Import results

heads_all = np.zeros((int(len(wave_names)), int(len(nodes.T)), int(len(t)), int(nrow), int(ncol)))

# create list for heads: [wave][node][time,row,col]
for wave in range(len(wave_names)):
    for node in range(len(nodes.T)):
        heads_all[wave,node,:,:,:] = np.load('./outputs/%s/heads_%d.npy' % (wave_names[wave], node))
        print('Gathering groundwater heads: wave=%s, node=%d' % (wave_names[wave], node))

f_approx_gwh = [None] * len(wave_names)
f_approx_okw = [None] * len(wave_names)



# %% Get expansion in one cell at all times

# choose specific cell for the analyisis
cell_rows = [50]
cell_cols = [5]

for wave in range(len(wave_names)):
        for row in cell_rows:
            print('Computing polynomial expansions (wave:%s, row:%d)' % (wave_names[wave], row+1))
            for col in cell_cols:                
                cell_eval_gwh = heads_all[wave, :, :, row, col]
                f_approx_gwh[wave] = cp.fit_quadrature(polynomial_expansion, nodes, weights, cell_eval_gwh)

del heads_all

                       
# %% Plot: Groundwater heads uncertainty PDFs
 
w_labels = ['Sine wave', 'Trapezoid wave', 'Triangular wave', 
            # 'Complex wave'
            ]

x = np.arange(-2, 2, 0.01)
time_plot = np.arange(0,1,1)
fp = f * 3600  # in hours
t_c = np.linspace(0, (len(t)/240*3), len(t))
tau = 1/fp
dp = (duration * 24)  # in hours
x_ticks = np.linspace(0, dp/tau, len(t))

coolors = ['royalblue', 'darkorange', 'tomato', 'green']

for ti in time_plot:
    print('Plotting t: %d...' % (ti))
    fig = plt.figure(ti, figsize=[3.5,3.5])
    plt.title(r'$\tau = %d/f$ (row:%d, column:%d)' % (ti, cell_rows[0], cell_cols[0]), loc='left')
    for wave in range(len(wave_names)):
        qoi_dist = cp.QoI_Dist(f_approx_gwh[wave][ti], J_distro, sample=10000)
        
        a = qoi_dist.pdf(x)
        a = a / np.sum(a)
        
        plt.plot(x, a, alpha=0.99, 
                 label=w_labels[wave], 
                 linewidth=1.0, color=coolors[wave])
        
    plt.legend(loc='upper right', borderpad=1, edgecolor='black', fancybox=0, fontsize='x-small', framealpha=1)
    plt.xlim(-2,2)
    plt.ylim(0, 0.07)
    # plt.ylim(0, 6)
    plt.xlabel('Groundwater head')
    plt.ylabel('P[h]')
    plt.tight_layout()
    plt.savefig('./results/uq_heads/pdf_ts%s.png' % (ti))
    plt.close(fig) 

