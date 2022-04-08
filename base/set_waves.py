# -*- coding: utf-8 -*-
"""
Created on Mon Nov 15 14:50:44 2021
CREATE CHD DATA (STOCHASTIC)
@author: PMR
"""

# %% Import libraries

# import external libraries
import numpy as np
import matplotlib.pyplot as plt

# import functions
import f_wave_functions as nf

# import inner modules
from set_sim import t, nper, nrow, nodes, wave_names, f

# fonts and colors
import matplotlib
font = {'size'   : 10}
matplotlib.rc('font', **font)
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family'] = "sans-serif"


# %% Compute stochastic waves

# create empty storage lists
p_heads_sine = [None] * len(nodes.T)
p_heads_trapezoid = [None] * len(nodes.T)
p_heads_triangle = [None] * len(nodes.T)
p_heads_composed = [None] * len(nodes.T)
p_heads_square = [None] * len(nodes.T)
p_heads_sawtooth_r0 = [None] * len(nodes.T)
p_heads_sawtooth_r1 = [None] * len(nodes.T)

 
for i in range(len(nodes.T)):
        
    # value for realization
    A = nodes[0,i]
    phase_e = nodes[1,i]
    
    # deterministic parameters
    p = 1 / (f * 3600)
    masl = 0  # meters above sea level
    
    # initial phase + random noise
    phase = 0 + phase_e
    triangular_phase = (p * 3600 / 4) +  phase_e
    saw_phase = (p * 3600 / 2) + phase_e

    # compute the wave functions        
    p_heads_sine[i] = nf.sine_wave(f, t, A, phase) + masl 
    p_heads_trapezoid[i] = nf.trapezoid_wave(f, t, A, phase, ext=6) + masl
    p_heads_triangle[i] = nf.triangular_wave(f, t, A, phase) + masl
    # p_heads_composed[i] = ((nf.sine_wave(f*3, t, A=(2*A)/5, phase=phase*3)) + (nf.sine_wave(f*1, t, A=A, phase=phase))) + masl
    
    # p_heads_square[i] = nf.square_wave(f, t, A, phase=phase) + masl
    # p_heads_sawtooth_r0[i] = nf.sawtooth_wave(f, t, A, phase=saw_phase, ramp=0) + masl
    # p_heads_sawtooth_r1[i] = nf.sawtooth_wave(f, t, A, phase=saw_phase, ramp=1) + masl
    

# gather all in one list
p_heads_all_waves = [p_heads_sine, p_heads_trapezoid, p_heads_triangle, p_heads_composed, p_heads_square, p_heads_sawtooth_r0, p_heads_sawtooth_r1]

# choose only active waves
p_heads_all_waves = p_heads_all_waves[0:len(wave_names)]


# %% Create synthetic CHD data

waves_chd = [[None] * len(nodes.T)] * len(wave_names)

for wave in range(len(wave_names)):
    
    for node in range(len(nodes.T)):
    
        # create CHD empty array (empty at every loop evaluation)
        chd_array = np.zeros(((nper*nrow), 6))
        chd_row = np.arange(0, nrow, 1)
        chd_col = [0]
        chd_data = {}
        
        # choose specific wave for the realization
        p_heads = p_heads_all_waves[wave][node]

        # create array for the boundary conditions
        i = 0
        for si in range(nper):
            for ri in chd_row:
                for ci in chd_col:
                    chd_array[i,0] = si
                    chd_array[i,1] = 0
                    chd_array[i,2] = ri
                    chd_array[i,3] = ci
                    chd_array[i,4] = p_heads[si]
                    chd_array[i,5] = p_heads[si]
                    i = i + 1
      
        # create CHD package objects 
        for i in range(nper):
            ii = i * nrow
            iii = ii + nrow
            chd_data[i] = chd_array[ii:iii,1:]
        
        # create list with all simulations
        waves_chd[wave][node] = chd_data


# %% Plots

for wave in range(len(wave_names)):
    plt.figure('Evaluations (%s)' % (wave_names[wave]))
    plt.title('Evaluations (%s)' % (wave_names[wave]), loc='left')
    for node in range(len(nodes.T)):
        plt.plot(t, p_heads_all_waves[wave][node], linewidth=0.5, color='blue', alpha=.1)
    plt.ylabel('CHD heads')
    plt.xlabel('Time')
