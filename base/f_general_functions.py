# -*- coding: utf-8 -*-
"""
Created on Mon Nov 15 14:28:41 2021
GENERAL FUNCTIONS
@author: PMR
"""


import matplotlib.pyplot as plt


def model_summary(Lx, Ly, ztop, zbot, ncol, nrow, nlay, sy, ss, hk_mu, hk_sigma):
    print('\n')
    print('*' * 50)
    print('Model setup summary:')
    print('• Lenght in x-axis: ............... %d' % (Lx))
    print('• Lenght in y-axis: ............... %d' % (Ly))
    print('• Surface level: .................. %d' % (ztop))
    print('• Aquifer bottom: ................. %d' % (zbot))
    print('• Number of columns: .............. %d' % (ncol))
    print('• Number of rows: ................. %d' % (nrow))
    print('• Number of layers: ............... %d' % (nlay))
    print('• Specific yield: ................. %.3f' % (sy))
    print('• Specific storage: ............... %.3e' % (ss))
    print('• Hyd. conductivity (mean): ....... %.3e' % (10**hk_mu))
    print('• Hyd. conductivity (std): ........ %.3e' % (10**hk_sigma))
    print('*' * 50)
    
      
def plot_condField(hk_mu, hk_sigma, hk_log):
    plt.figure('Random hydraulic conductivity field')
    plt.title(('Random hydraulic conductivity field ($\mu$=%.2f, $\sigma$=%.2f)' % (hk_mu, hk_sigma)))
    plt.imshow(hk_log, cmap='rainbow')
    plt.colorbar(orientation='horizontal', shrink=0.75)


def plot_waveFunctions_sep(wave_names, p_heads_all_waves, t):
    for i in range(len(p_heads_all_waves)):
        plt.figure('Specified heads - Type: %s' % (wave_names[i]))
        plt.title('Specified heads - Type: %s' % (wave_names[i]))
        plt.plot(t, p_heads_all_waves[i], linewidth=.75)
        plt.xlim(0,200000)
        plt.ylim(3,7)
     

def plot_waveFunctions_all(wave_names, p_heads_all_waves, t):
    plt.figure('Specified heads - Type: All')
    plt.title('Specified heads - Type: All')
    for i in range(len(p_heads_all_waves)):
        plt.plot(t, p_heads_all_waves[i], linewidth=1., alpha=0.4, label=wave_names[i])
    plt.xlim(0,200000)
    plt.ylim(3,7)
    plt.legend()
