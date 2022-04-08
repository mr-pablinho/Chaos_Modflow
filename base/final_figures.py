# -*- coding: utf-8 -*-
"""
Created on Tue Dec 28 13:53:06 2021
FINAL PLOTS - STOCHASTIC WAVES v2
@author: PMR
"""

# %% Import libraries

# import external libraries
import numpy as np
import matplotlib.pyplot as plt

# fonts and colors
import matplotlib
font = {'size'   : 10}
matplotlib.rc('font', **font)
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family'] = "sans-serif"

# import inner modules and variables
from set_sim import t, nrow, ncol, wave_names, duration, f


# %% Import results

E_gwh = np.zeros((len(wave_names), 240, nrow, ncol))
S_gwh = np.zeros((len(wave_names), 240, nrow, ncol))
E_okw = np.zeros((len(wave_names), 240, nrow, ncol))
S_okw = np.zeros((len(wave_names), 240, nrow, ncol))

for wave in range(len(wave_names)):
    E_gwh[wave,:,:,:] = np.load('./outputs/uq_heads/E_heads__%s.npy' % (wave_names[wave]))[:240,:,:]
    S_gwh[wave,:,:,:] = np.load('./outputs/uq_heads/S_heads__%s.npy' % (wave_names[wave]))[:240,:,:]
    E_okw[wave,:,:,:] = np.load('./outputs/uq_topo/E_okw__%s.npy' % (wave_names[wave]))[:240,:,:]
    S_okw[wave,:,:,:] = np.load('./outputs/uq_topo/S_okw__%s.npy' % (wave_names[wave]))[:240,:,:]



# %% Labels
w_labels = ['Sine wave', 'Trapezoid wave', 'Triangular wave', 
            # 'Complex wave'
            ]


# %% Fig: UQ in groundwater heads and Okubo-Weiss


# create ticks to plot
fp = f * 3600  # in hours
t_c = np.linspace(0, (len(t)/240*3), len(t))
tau = 1/fp
dp = (duration * 24)  # in hours
x_ticks = np.linspace(0, dp/tau, len(t))
n_tau = 3

# choose row and column to evaluate
choose_row = 50
choose_col = 5

# plot
fig, axs = plt.subplots(2, len(w_labels))
from_tau, to_tau = 0,240


for wave in range(len(w_labels)):
    
    # UQ groundwater heads
    lineE_gwh = np.mean(E_gwh[wave, from_tau:to_tau, :, choose_col], axis=1)
    lineS_gwh = np.mean(S_gwh[wave, from_tau:to_tau, :, choose_col], axis=1)
    
    axs[0,wave].set_title(w_labels[wave], loc='left')
    axs[0,wave].fill_between(x_ticks[from_tau:to_tau], 
                             lineE_gwh - lineS_gwh, 
                             lineE_gwh + lineS_gwh,
                             color='red', alpha=0.3, edgecolor=None,
                             label=r'[$\mu_H - \sigma_H, \mu_H + \sigma_H$]')
    
    axs[0,wave].plot(x_ticks[from_tau:to_tau], 
                     lineE_gwh,
                     linewidth=1,
                     color='red',
                     label=r'$E_H$')
    
    ymin, ymax = -1.8, 1.8
    
    axs[0,wave].set_xlim(0, n_tau)
    axs[0,wave].set_ylim(ymin, ymax)
    axs[0,wave].vlines(x=np.arange(1,n_tau,1), ymin=ymin, ymax=ymax, 
                       linestyle='--', alpha=.5, color='black', linewidth=1.0)
    axs[0,wave].hlines(xmin=0, xmax=n_tau, y=0, 
                       linestyle='-', alpha=.99, color='black', linewidth=.70)

    # UQ Okubo-Weiss
    lineE_okw = np.mean(E_okw[wave, from_tau:to_tau, :, choose_col], axis=1)
    lineS_okw = np.mean(S_okw[wave, from_tau:to_tau, :, choose_col], axis=1)
    
    axs[1,wave].fill_between(x_ticks[from_tau:to_tau], 
                             lineE_okw - lineS_okw, 
                             lineE_okw + lineS_okw,
                             color='blue', alpha=0.3, edgecolor=None,
                             label=r'[$\mu_O - \sigma_O, \mu_O + \sigma_O$]')   
    
    axs[1,wave].plot(x_ticks[from_tau:to_tau], 
                     lineE_okw,
                     linewidth=1,
                     color='blue',
                     label=r'$E_O$')
    
    ymin, ymax = -1e-6, 3e-6
    
    axs[1,wave].set_xlim(0, n_tau)
    axs[1,wave].set_ylim(ymin, ymax)
    axs[1,wave].set_xlabel(r'$\tau$')
    axs[1,wave].set_yticks([-1*1e-6, 0, 1*1e-6, 2*1e-6, 3*1e-6])
    
    axs[1,wave].vlines(x=np.arange(1,n_tau,1), ymin=ymin, ymax=ymax, 
                       linestyle='--', alpha=.5, color='black', linewidth=1.0)
    axs[1,wave].hlines(xmin=0, xmax=n_tau, y=0, 
                       linestyle='-', alpha=.99, color='black', linewidth=0.70)  

    
axs[0,0].legend(loc='upper left', borderpad=1, edgecolor='black', fancybox=0, fontsize='x-small', framealpha=1)
axs[1,0].legend(loc='upper left', borderpad=1, edgecolor='black', fancybox=0, fontsize='x-small', framealpha=1)
axs[0,0].set_ylabel('Groundwater head')
axs[1,0].set_ylabel('Okubo-Weiss')

fig.tight_layout()
fig.set_figheight(5.5)
fig.set_figwidth(8)
