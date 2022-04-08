# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 10:05:05 2021
Last update on Thu Nov 25 15:56:21 2021
FLOW FIELD EVALUATION
"""

# %% Import libraries

# import external libraries
import numpy as np
import chaospy as cp

# import inner modules and variables
import f_flow_field_eval as fpf
from set_sim import t, delr, delc, hk, nodes, weights, nrow, ncol, wave_names, polynomial_expansion, J_distro



# %% Compute flow field and metrics

# create storage arrays
qx = np.zeros((len(nodes.T), len(t), nrow, ncol))
qy = np.zeros((len(nodes.T), len(t), nrow, ncol))
ff_metrics = np.zeros((len(nodes.T), 4, len(t), nrow, ncol))
okuboweiss = np.zeros((len(nodes.T), len(t), nrow, ncol))
# vorticity = np.zeros((len(nodes.T), len(t), nrow, ncol))
# shear = np.zeros((len(nodes.T), len(t), nrow, ncol))
# stretch = np.zeros((len(nodes.T), len(t), nrow, ncol))

for wave in range(len(wave_names)):
    print('Computing topology metrics: %s wave...' % (wave_names[wave]))
    
    for node in range(len(nodes.T)): 
        
        # choose head array to evaluate
        heads_eval = np.load('.\outputs\%s\heads_%d.npy' % (wave_names[wave], node))
        
        # compute Darcy's velocity
        qx[node, :, :, :] = fpf.darcy_velocity(heads_eval, hk_field=hk, delta_x=delc, delta_y=delr)[0]
        qy[node, :, :, :] = fpf.darcy_velocity(heads_eval, hk_field=hk, delta_x=delc, delta_y=delr)[1]
        
        # compute topological
        ff_metrics[node, :, :, :, :] = fpf.okubo_weiss_2d(qx[node, :, :, :],
                                                          qy[node, :, :, :],
                                                          heads_eval, 
                                                          delta_x=delc,
                                                          delta_y=delr)
        
        # separate metrics 
        okuboweiss[node, :, :, :] = ff_metrics[node, 0, :, :, :]
        # vorticity[wave, node, :, :, :] = ff_metrics[node, 1, :, :, :]
        # shear[wave, node, :, :, :] = ff_metrics[node, 2, :, :, :]
        # stretch[wave, node, :, :, :] = ff_metrics[node, 3, :, :, :]


    np.save('./outputs/uq_topo/evals_okw__%s' % (wave_names[wave]), okuboweiss)
    # np.save('./outputs/uq_topo/evals_vor__%s' % (wave_names[wave]), vorticity)
    # np.save('./outputs/uq_topo/evals_she__%s' % (wave_names[wave]), shear)
    # np.save('./outputs/uq_topo/evals_str__%s' % (wave_names[wave]), stretch)


# %% Compute uncertainty with PCE

    E_okw, S_okw = np.zeros((len(t), nrow, ncol)), np.zeros((len(t), nrow, ncol))
    # E_vor, S_vor = np.zeros((len(wave_names), len(t), nrow, ncol)), np.zeros((len(wave_names), len(t), nrow, ncol))
    # E_she, S_she = np.zeros((len(wave_names), len(t), nrow, ncol)), np.zeros((len(wave_names), len(t), nrow, ncol))
    # E_str, S_str = np.zeros((len(wave_names), len(t), nrow, ncol)), np.zeros((len(wave_names), len(t), nrow, ncol))
    
    period = np.arange(0, len(t), 1)
    
    # period = np.linspace(0, 250-1, 83).astype(int)

    print('Computing polynomial chaos expansions - Topology metrics')


    for row in range(nrow):
        print('Computing polynomial expansions (wave:%s, row:%d)' % (wave_names[wave], row+1))
        for col in range(ncol):

            okw_eval = okuboweiss[:, :, row, col]
            # vor_eval = vorticity[wave, :, :, row, col]
            # she_eval = shear[wave, :, :, row, col]
            # str_eval = stretch[wave, :, :, row, col]
            
            f_approx_okw = cp.fit_quadrature(polynomial_expansion, nodes, weights, okw_eval)
            # f_approx_vor = cp.fit_quadrature(polynomial_expansion, nodes, weights, vor_eval)
            # f_approx_she = cp.fit_quadrature(polynomial_expansion, nodes, weights, she_eval)
            # f_approx_str = cp.fit_quadrature(polynomial_expansion, nodes, weights, str_eval)
            
            E_okw[:, row, col], S_okw[:, row, col] = cp.E(f_approx_okw, J_distro), cp.Std(f_approx_okw, J_distro)
            # E_vor[wave, :, row, col], S_vor[wave, :, row, col] = cp.E(f_approx_vor, J_distro), cp.Std(f_approx_vor, J_distro)
            # E_she[wave, :, row, col], S_she[wave, :, row, col] = cp.E(f_approx_she, J_distro), cp.Std(f_approx_she, J_distro)
            # E_str[wave, :, row, col], S_str[wave, :, row, col] = cp.E(f_approx_str, J_distro), cp.Std(f_approx_str, J_distro)      


# %% Save arrays

    np.save('./outputs/uq_topo/E_okw__%s' % (wave_names[wave]), E_okw)
    # np.save('./outputs/uq_topo/E_vor', E_vor)
    # np.save('./outputs/uq_topo/E_she', E_she)
    # np.save('./outputs/uq_topo/E_str', E_str)
    
    np.save('./outputs/uq_topo/S_okw__%s' % (wave_names[wave]), S_okw)
    # np.save('./outputs/uq_topo/S_vor', S_vor)
    # np.save('./outputs/uq_topo/S_she', S_she)
    # np.save('./outputs/uq_topo/S_str', S_str)


