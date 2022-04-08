# -*- coding: utf-8 -*-
"""
Created on Sat Dec 18 20:56:29 2021
Uncertainty propagtion in the CHD  -- Monte Carlo test
@author: PMR
"""


# %% Import libraries

# import external libraries
import numpy as np
import chaospy as cp
import matplotlib.pyplot as plt

# import functions
import f_wave_functions as nf

# fonts and colors
import matplotlib
font = {'size'   : 10}
matplotlib.rc('font', **font)
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family'] = "sans-serif"


# %% Setup stochastic problem

# time discretization
duration = 3  # in days
inputs_day = 240  # inputs per unit of duration
nper = duration * inputs_day  # number of stress periods
delta_t_h = 24/inputs_day  # in hours
delta_t = delta_t_h * 3600  # in seconds
perlen = [delta_t] * nper  # simulation period length
nstp = [1] * nper  # time steps per stress periods
nstp[0] = 1  # to edit in detail the number of time steps in a specific stress period
steady = [False] * nper  # steady or transient state

# sampling rate, or number of measurements per second for wave functions
ms = nper  # time coordinates
last_second = duration * 24 * 3600
t = np.linspace(0, last_second, ms, endpoint=False)  # time coordinates


# amplitude
a_A = 0.90
b_A = 1.10
A_distro = cp.Uniform(a_A, b_A)

# frequency [cycles per second]
# mu_err = 0
# std_err = 1e-6
# err_distro = cp.Normal(mu_err, std_err)

a_phase = 2 * (-3600)
b_phase = 2 * (3600)
phase_distro = cp.Uniform(a_phase, b_phase)


# Monte Carlo test
nSamples = 1e4
# rule = 'Random'
rule = 'Halton'
# rule = 'Sobol'

# joint distribution
# J_distro = cp.J(A_distro, err_distro)
J_distro = cp.J(A_distro, phase_distro)
samples = J_distro.sample(nSamples, rule=rule)

# nSamples = 4
# samples = np.array(([1, 1, 1, 1],
#                     [1/(7*3600), 1/(7.25*3600), 1/(7.5*3600), 1/(9*3600)]))


# create empty storage lists
p_heads_sine = [None] * int(nSamples)
p_heads_trapzoid = [None] * int(nSamples)
p_heads_triangle = [None] * int(nSamples)
p_heads_composed = [None] * int(nSamples)
 
for i in range(int(nSamples)):
    
    # value for realization
    A = samples[0,i]
    phase_e = samples[1,i]
    f = 1 / (8 * 3600) 
    
    phase=0+phase_e
    p = 1 / (f * 3600)
    triangular_phase = (p * 3600 / 4) +  phase_e
    
    # err = samples[1,i]
    # f = 1 / (8 * 3600) + err
    
    # computation steps
    omega = nf.angular_frequency(f)
    
    
    # compute the wave functions
    
    masl = 0  # meters above sea level
    
    p_heads_sine[i] = nf.sine_wave(omega, t, A, phase=phase) + masl   
    p_heads_trapzoid[i] = nf.trapzoid_wave(f, t, A, phase=triangular_phase, ext=6) + masl
    p_heads_triangle[i] = nf.triangle_wave(f, t, A, phase=triangular_phase) + masl
    p_heads_composed[i] = ((nf.sine_wave(omega*3, t, A=A/3, phase=phase)) + (nf.sine_wave(omega*1, t, A=A, phase=phase))) + masl


p_heads_sine = np.array(p_heads_sine)
p_heads_trapzoid = np.array(p_heads_trapzoid)
p_heads_triangle = np.array(p_heads_triangle)
p_heads_composed = np.array(p_heads_composed)

t_c = np.arange(0,len(t),1)
for i in range(len(p_heads_sine)):
    plt.plot(t_c[:], p_heads_trapzoid[i,:], alpha=0.1, linewidth=1, color='purple')
plt.ylim(-1.5,1.5)
plt.xlim(0, t_c[-1])

# %% Compute UQ metrics

E_gwh_sine = np.mean(p_heads_sine, axis=0)
S_gwh_sine = np.std(p_heads_sine, axis=0)

E_gwh_trapzoid = np.mean(p_heads_trapzoid, axis=0)
S_gwh_trapzoid = np.std(p_heads_trapzoid, axis=0)

E_gwh_triangle = np.mean(p_heads_triangle, axis=0)
S_gwh_triangle = np.std(p_heads_triangle, axis=0)

E_gwh_composed = np.mean(p_heads_composed, axis=0)
S_gwh_composed = np.std(p_heads_composed, axis=0)



# %%

def plot_uq(E_gwh, S_gwh, t, w_title, N):
    
    comb = ' (N=%d, rule=%s) ' % (N, rule)
    t_c = np.arange(0,len(t),1)
    
    # plt.figure('E' + w_title, figsize=(5,4))
    # plt.plot(t_c, E_gwh, 
    #           linewidth=1, color='blue')
    # plt.title(w_title + comb, loc='left')
    # plt.ylabel('Expected value - Heads')
    # plt.xlabel('Stress period')
    # plt.xlim(0,t_c[-1])
    # plt.ylim(-1.3,1.3)
    # plt.savefig('./new_figs/MC/E__' + w_title + comb + '.png')
    # # plt.close() 
    
    # plt.figure('S' + w_title, figsize=(5,4))
    # plt.plot(t_c, S_gwh,
    #           linewidth=1, color='red')
    # plt.title(w_title + comb, loc='left')
    # plt.ylabel('Standard deviation - Heads')
    # plt.xlabel('Stress period')
    # plt.xlim(0,t_c[-1])
    # plt.ylim(0,1)
    # plt.savefig('./new_figs/MC/S__' + w_title + comb + '.png')
    # # plt.close() 
    
    plt.figure('UQ' + w_title, figsize=(5,4))
    plt.fill_between(t_c, E_gwh - 1*S_gwh, E_gwh + 1*S_gwh, alpha=0.25, color='red')
    plt.plot(t_c, E_gwh, 
             linewidth=1, color='black')
    plt.title(w_title + comb, loc='left')
    plt.ylabel('UQ - Heads')
    plt.xlabel('Stress period')
    plt.xlim(0,t_c[-1])
    plt.ylim(-1.5,1.5)
    plt.savefig('./new_figs/MC/UQ__' + w_title + comb + '.png')
    # plt.close() 


plot_uq(E_gwh_sine, S_gwh_sine, t, 'Sine', nSamples)
plot_uq(E_gwh_trapzoid, S_gwh_trapzoid, t, 'Trapzoid', nSamples)
plot_uq(E_gwh_triangle, S_gwh_triangle, t, 'Triangle', nSamples)
plot_uq(E_gwh_composed, S_gwh_composed, t, 'Composed', nSamples)