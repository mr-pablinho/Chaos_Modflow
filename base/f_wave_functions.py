# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 13:36:38 2021
UQ WAVES - FUNCTIONS
@author: PMR
"""

# %% Import libraries
import numpy as np
# from scipy import signal
# import matplotlib.pyplot as plt


# %% Wave functions


def sine_wave(f, t, A=1, phase=0):
    y = A * np.sin(2 * np.pi * f * (t + phase))
    return y


def square_wave(f, t, A=1, phase=0):
    y = A * np.sign(np.sin((2 * np.pi * f * t) + phase)) 
    return y


def trapezoid_wave(f, t, A, phase, ext):
    A_l = A * ext
    y = A_l * triangular_wave(f, t, A, phase)
    y[y >  A] = A
    y[y < -A] = -A
    return y


def triangular_wave(f, t, A, phase):
    # per = 1/f
    # y = 4 * (A/per) * np.abs(((((t - phase) - per/4 ) % per) + per ) % per - per/2) - A
    y = (4 * A * f) * np.abs(((((t - phase) - (1/(4 * f)) ) % (1/f)) + (1/f) ) % (1/f) - 1/(2*f)) - A
    return y
  

'''
def sine_wave(f, t, A=1, phase=0):
    y = A * np.sin((2 * np.pi * f * t) + phase)
    return y


def square_wave(f, t, A=1, phase=0):
    y = A * np.sign(np.sin((2 * np.pi * f * t) + phase)) 
    return y


def triangle_wave(f, t, A=1, phase=0):
    y = A * signal.sawtooth((2 * np.pi * f * t) + phase, width=0.5)
    return y


def sawtooth_wave(f, t, A=1, phase=0, ramp=1):
    y = A * signal.sawtooth((2 * np.pi * f * t) + phase, width=ramp)
    return y


def square_wave_fourier(f, t, A=1, phase=0, expansion_terms=10):
    terms = np.arange(1, expansion_terms + 1, 2)
    expansion = 0
    for i in terms:
        sin_term = (1/i) * (np.sin((2 * f * np.pi * i * t) + phase))
        expansion = expansion + sin_term
        expansion[400]
    y = A * (4/np.pi) * expansion
    return y


def trapezoid_wave(f, t, A=1, phase=0.25, ext=5):
    A_l = A * ext
    y = A_l * signal.sawtooth((2 * np.pi * f * t) + phase, width=0.5)
    y[y >  A] = A
    y[y < -A] = -A
    return y
  
'''


