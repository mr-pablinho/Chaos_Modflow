# -*- coding: utf-8 -*-
"""
Created on Tue Dec 28 11:22:21 2021
Uncertainty propagtion in the CHD  -- Monte Carlo test - Complex wave
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


# deterministic parameters
f = 1 / (8 * 3600)
p = 1 / (f * 3600)
omega = nf.angular_frequency(f)
masl = 0  # meters above sea level

# amplitude (big wave)
a_Ab = 0.9
b_Ab = 1.1
Ab_distro = cp.Uniform(a_Ab, b_Ab)

# amplitude (small wave)
a_As = 0.9 / 3
b_As = 1.1 / 3
As_distro = cp.Uniform(a_As, b_As)

# phase (big wave) [cycles per second]
a_phase_b = -1/(8*f)  # tau/8*t
b_phase_b =  1/(8*f)  # tau/8*t
phase_distro_b = cp.Uniform(a_phase_b, b_phase_b)

# phase (small wave) [cycles per second]
a_phase_s = -1/(24*f)  # tau/24*t
b_phase_s =  1/(24*f)  # tau/24*t
phase_distro_s = cp.Uniform(a_phase_s, b_phase_s)


# Monte Carlo test
nSamples = 1e2
# rule = 'Random'
rule = 'Halton'
# rule = 'Sobol'

# joint distribution
J_distro = cp.J(Ab_distro, As_distro, phase_distro_b, phase_distro_s)
samples = J_distro.sample(nSamples, rule=rule)

p_heads_composed = [None] * int(nSamples)
 
for i in range(int(nSamples)):
    
    # value for realization
    Ab = samples[0,i]
    As = samples[1,i]
    phase_e_b = samples[2,i]
    phase_e_s = samples[3,i]
    
    omega = nf.angular_frequency(f)
    
    # initial phase + random noise
    phase_b = 0 + phase_e_b
    phase_s = 0 + phase_e_s

    # compute the wave functions
    p_heads_composed[i] = ((nf.sine_wave(omega*3, t, A=As, phase=phase_s)) + (nf.sine_wave(omega*1, t, A=Ab, phase=phase_b))) + masl
    
p_heads_composed = np.array(p_heads_composed)


# %% Compute UQ metrics

E_gwh_composed = np.mean(p_heads_composed, axis=0)
S_gwh_composed = np.std(p_heads_composed, axis=0)



# %%

def plot_uq(E_gwh, S_gwh, t, w_title, N):
    
    comb = ' (N=%d, rule=%s) ' % (N, rule)
    t_c = np.arange(0,len(t),1)
    
    plt.figure('E' + w_title, figsize=(5,4))
    plt.plot(t_c, E_gwh, 
              linewidth=1, color='blue')
    plt.title(w_title + comb, loc='left')
    plt.ylabel('Expected value - Heads')
    plt.xlabel('Stress period')
    plt.xlim(0,t_c[-1])
    plt.ylim(-1.3,1.3)
    plt.savefig('./new_figs/MC/E__' + w_title + comb + '.png')
    # # plt.close() 
    
    plt.figure('S' + w_title, figsize=(5,4))
    plt.plot(t_c, S_gwh,
              linewidth=1, color='red')
    plt.title(w_title + comb, loc='left')
    plt.ylabel('Standard deviation - Heads')
    plt.xlabel('Stress period')
    plt.xlim(0,t_c[-1])
    plt.ylim(0,1)
    plt.savefig('./new_figs/MC/S__' + w_title + comb + '.png')
    # # plt.close() 
    
    plt.figure('UQ' + w_title, figsize=(5,4))
    plt.fill_between(t_c, E_gwh - 1*S_gwh, E_gwh + 1*S_gwh, alpha=0.25, color='red')
    plt.plot(t_c, E_gwh, 
             linewidth=1, color='black')
    plt.title(w_title + comb, loc='left')
    plt.ylabel('UQ - Heads')
    plt.xlabel('Stress period')
    plt.xlim(0,240)
    plt.ylim(-1.5,1.5)
    plt.savefig('./new_figs/MC/UQ__' + w_title + comb + '.png')
    # plt.close() 


plot_uq(E_gwh_composed, S_gwh_composed, t, 'Composed', nSamples)