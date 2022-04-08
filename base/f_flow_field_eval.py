# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 18:20:29 2021
Last update on on Wed Nov 24 16:18:43 2021
DARCY VELOCITY & OKUBO-WEISS 
@author: PMR
"""

import numpy as np
import matplotlib.pyplot as plt



def darcy_velocity(output_heads, hk_field, delta_x, delta_y):
    
    # extend the array to compute the gradient   
    x_head_ex = np.hstack((output_heads, np.tile(output_heads[:, [-1]], 1)))
    y_head_ex = np.vstack((output_heads, np.tile(output_heads[[-1], :], 1)))
    
    # difference of heads
    xoo_h = x_head_ex[:,1:] - x_head_ex[:,0:-1]
    yoo_h = y_head_ex[1:,:] - y_head_ex[0:-1,:]
    
    # evaluate if delta of the distance was provide as an array or number
    if type(delta_x) == np.ndarray:
        delta_x_array = delta_x
    else:
        delta_x_array = np.full(output_heads.shape, delta_x)
        
    if type(delta_y) == np.ndarray:
        delta_y_array = delta_y
    else:
        delta_y_array = np.full(output_heads.shape, delta_y)
        
    # compute Darcy velocity
    qx = -(np.multiply(hk_field, (np.divide(xoo_h, delta_x_array))))
    qy = -(np.multiply(hk_field, (np.divide(yoo_h, delta_y_array))))
    
    return qx, qy


def okubo_weiss_2d(qx, qy, output_heads, delta_x, delta_y):
    
    # extend the velocity array to compute the derivatives   
    vx = np.hstack((qx, np.tile(qx[:, [-1]], 1)))
    vx = np.vstack((vx, np.tile(vx[[-1], :], 1)))
    
    vy = np.hstack((qy, np.tile(qy[:, [-1]], 1)))
    vy = np.vstack((vy, np.tile(vy[[-1], :], 1)))
    
    
    # difference of velocities
    dvx_x = (vx[:,1:] - vx[:,0:-1])[0:-1,:]
    dvx_y = (vx[1:,:] - vx[0:-1,:])[:,0:-1]
    dvy_x = (vy[:,1:] - vy[:,0:-1])[0:-1,:]
    dvy_y = (vy[1:,:] - vy[0:-1,:])[:,0:-1]
           
    # evaluate if delta of the distance was provide as an array or number
    if type(delta_x) == np.ndarray:
        delta_x_array = delta_x
    else:
        delta_x_array = np.full(dvx_x.shape, delta_x)
        
    if type(delta_y) == np.ndarray:
        delta_y_array = delta_y
    else:
        delta_y_array = np.full(dvy_y.shape, delta_y)
   
    
    # compute derivatives
    dvx_dx = np.divide(dvx_x, delta_x_array)
    dvx_dy = np.divide(dvx_y, delta_y_array)
    dvy_dy = np.divide(dvy_y, delta_y_array)
    dvy_dx = np.divide(dvy_x, delta_x_array)
    
    # curl
    curl2d =  dvy_dx - dvx_dy
    
    # shear
    shear2d =  dvy_dx + dvx_dy
    
    # streatch
    stretch2d = 2 * dvx_dx
    
    # Okubo-Weiss
    okuboweiss2d = (np.power(stretch2d, 2) + np.power(shear2d, 2)) - np.power(curl2d, 2)
    
    return okuboweiss2d, curl2d, shear2d, stretch2d
