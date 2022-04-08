# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 16:39:42 2021
INITIAL SETUP
@author: PMR
"""

# %% Import libraries

# import external libraries
import numpy as np
import chaospy as cp
import matplotlib.pyplot as plt


# %% Model discretization

# model extent
Lx = 100.0  # x extension [m]
Ly = 100.0  # y extension [m]
ztop = 10.0  # surface level [m a.s.l.]
zbot = -10.0  # bottom of the aquifer [m a.s.l.]

# grid resolution
factor_col = 1  # factor to increase the spatial discretization
factor_row = 1  # factor to increase the spatial discretization
nlay = 1  # number of layers
nrow = int(Ly * factor_col)  # number of rows
ncol = int(Lx * factor_col)  # number of columns 
delr = Lx / ncol  # cell size in x
delc = Ly / nrow  # cell size in y 
delv = (ztop - zbot) / nlay  # cell size in z
botm = np.linspace(ztop, zbot, nlay + 1)  # layer height 

# variables for the BAS package
ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
init_heads = 0  # initial conditions [m a.s.l.]
strt = init_heads * np.ones((nlay, nrow, ncol), dtype=np.float32)

# time discretization
duration = 3 # in days
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


# %% Hydraulic parameters

# hydraulic parameters
laytyp = 1  # layer type: 0 confined or > 0 convertible
sy = 0.25  # specific yield and specific storage [1/m]
ss = 6e-05  # specific yield and specific storage [-]    

# choose hydraulic conductivity field [m/s]
hk = np.loadtxt('./inputs/hk_s1.5_333.txt')

hk = hk[0:nrow, 0:ncol]


# %% Setup stochastic problem

# frequency
f = 1 / (8 * 3600)

# amplitude
a_A = 0.9
b_A = 1.1
A_distro = cp.Uniform(a_A, b_A)

# phase [seconds]
a_phase = -1/(8*f)
b_phase =  1/(8*f)
phase_distro = cp.Uniform(a_phase, b_phase)

# joint distribution
J_distro = cp.J(A_distro, phase_distro)

# polynomial chaos expansions
quadOrder = 9
polyOrder = 3
rule_nodes = 'Gaussian'

nodes, weights = cp.generate_quadrature(order=quadOrder, 
                                        dist=J_distro, 
                                        sparse=0, 
                                        rule=rule_nodes)

polynomial_expansion = cp.expansion.stieltjes(polyOrder, 
                                              J_distro, 
                                              normed=True)

# %% Plot nodes

plt.figure('Nodes and weights', figsize=(4,4))
plt.scatter(*nodes, s=weights*4e3, alpha=0.3, color='purple')
plt.title('Nodes and weights', loc='left')
plt.xlabel('Amplitude')
plt.ylabel('Phase shift')


# %% Choose periodic functions to evaluate

wave_names = ['sine', 'trapezoid',  'triangle']
# wave_names = ['sine']
