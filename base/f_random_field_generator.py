# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 13:57:40 2021
RANDOM FIELD GENERATOR 
@author: Pablo Merchán-Rivera (2021)
"""

import chaospy as cp
import numpy as np



def random_field(Nx, Ny, mu, sigma, dist_type=None, corr_x=1, corr_y=1, cov_function='exponential', verbose=False):
    
    '''
    Example:
    mu, sigma = 0, 1
    dist_type = 'Normal'
    Nx = 100
    Ny = 20
    corr_x = 5
    corr_y = 1
    cov_function = 'both'
    field = random_field(Nx, Ny, mu, sigma, dist_type, corr_x, corr_y, cov_function, verbose=False)
    plt.figure('0')
    plt.imshow(field[0], cmap='rainbow')
    plt.colorbar()
    plt.figure('1')
    plt.imshow(field[1], cmap='rainbow')
    plt.colorbar()
    '''
    
    print('Creating random field (%dx%d)... ' % (Nx, Ny))
    
    x_min, x_max, y_min, y_max = 0.0, 1*Nx/corr_x, 0.0, 1*Ny/corr_y
    mesh_size_x, mesh_size_y = Nx, Ny
    
    start_x, stop_x, step_x = x_min + x_max/(2*mesh_size_x), x_max, (x_max - x_min)/mesh_size_x
    start_y, stop_y, step_y = y_min + y_max/(2*mesh_size_y), y_max, (y_max - y_min)/mesh_size_y
    
    x_coord = np.arange(start_x, stop_x, step_x)
    y_coord = np.arange(start_y, stop_y, step_y)
    
    if verbose == True:
        message_random_field(start_x, stop_x, step_x, x_coord, start_y, stop_y, step_y, y_coord)
    
    mesh_coord = [None] * (Nx*Ny)
    ii = 0
    for i in range(mesh_size_x):
        for j in range(mesh_size_y):
            mesh_coord[ii] = x_coord[i], y_coord[j]
            ii += 1
    
    C1 = np.zeros(shape=(Nx*Ny, Nx*Ny))
    C2 = np.zeros(shape=(Nx*Ny, Nx*Ny))
    m = np.zeros(shape=Nx*Ny)
    
    for i in range(Nx*Ny):
        for j in range(Nx*Ny):
            C1[i, j] = c1(mesh_coord[i], mesh_coord[j])
            C2[i, j] = c2(mesh_coord[i], mesh_coord[j])
        C2[i, i] += 1E-12
        m[i] = 0.0
    
    
    L1 = np.linalg.cholesky(C1)
    L2 = np.linalg.cholesky(C2)
      
    
    if dist_type == 'normal' or dist_type == 'Normal' or dist_type == 'gaussian' or dist_type == 'Gaussian':
        print('Creating gaussian random field (mu=%.3f, sigma=%.3f)' % (mu, sigma))
        dist = cp.Normal(0, sigma).sample((Nx*Ny))
        np.random.shuffle(dist)
        G1 = mu + np.matmul(L1, dist)
        G2 = mu + np.matmul(L2, dist)

    
    if cov_function == 'exponential' or cov_function == 'Exponential':
        foo = np.reshape(G1, (Nx,Ny)).T
    elif cov_function == 'squared exponential' or cov_function == 'Squared exponential':
        foo = np.reshape(G2, (Nx,Ny)).T
    elif cov_function == 'both' or cov_function == 'Both':
        foo = [np.reshape(G1, (Nx,Ny)).T, np.reshape(G2, (Nx,Ny)).T]
    else:
        print('The covariance function was not defined. Using exponential covariance function by default.')
        foo = np.reshape(G1, (Nx,Ny)).T

    return foo



# %% Covariance functions


# exponential 
def c1(x, y):
    d = (x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1])
    return np.exp(-(np.sqrt(d)))


# squared exponential
def c2(x, y):
    d = (x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1])
    return np.exp(-d*0.5)



# %% Verbose message


def message_random_field(start_x, stop_x, step_x, x_coord, start_y, stop_y, step_y, y_coord):
    print('Mesh for random field:')
    print('x-coord -->  start: %.3f, stop: %.3f, step: %.3f, length:%.3f' % (start_x, stop_x, step_x, len(x_coord)))
    print('y-coord -->  start: %.3f, stop: %.3f, step: %.3f, length:%.3f' % (start_y, stop_y, step_y, len(y_coord)))
    print('---'*10)
    print('x-coord array --> %s' % (x_coord))
    print('y-coord array --> %s' % (y_coord))
    
