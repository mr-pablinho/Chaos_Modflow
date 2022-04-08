# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 14:18:57 2021
Last update on Mon Nov 29 14:18:57 2021
UQ IN PERIODIC FUNCTIONS
@author: PMR
"""

# %% Import libraries

# import external libraries
import numpy as np
import chaospy as cp
import matplotlib.pyplot as plt

# import inner modules and variables
import f_flow_field_eval as fpf
from set_sim import t, delr, delc, hk, nodes, weights, nrow, ncol, wave_names, polynomial_expansion, J_distro
from set_waves import p_heads_all_waves


# %% Import results and variables

E_per = np.zeros((len(wave_names), len(t)))
S_per = np.zeros((len(wave_names), len(t)))

# groundwater heads

for wave in range(len(wave_names)):
    
    eval_perFunction = np.array(p_heads_all_waves[wave])
    f_approx_per = cp.fit_quadrature(polynomial_expansion, nodes, weights, eval_perFunction)
    E_per[wave] = cp.E(f_approx_per, J_distro)
    S_per[wave] = cp.Std(f_approx_per, J_distro)
    
    for node in range(len(nodes.T)):
        plt.figure('Evaluations %s' % (wave_names[wave]))
        plt.title('Evaluations %s' % (wave_names[wave]))
        plt.plot(t, p_heads_all_waves[wave][node][:], color='purple', alpha=.3)


