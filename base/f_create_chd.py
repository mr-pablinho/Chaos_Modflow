# -*- coding: utf-8 -*-
"""
Created on Wed Nov 17 18:25:26 2021
CREATE CHD DICTIONARY
@author: PMR
"""

import numpy as np


def chd_dictionary(p_heads, nper, nrow):
    
    # create CHD empty array (empty at every loop evaluation)
    chd_array = np.zeros(((nper*nrow), 6))
    chd_row = np.arange(0, nrow, 1)
    chd_col = [0]
    chd_data = {}
    
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
        
    return chd_data