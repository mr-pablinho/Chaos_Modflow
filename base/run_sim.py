# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 14:47:24 2021
FLOPY
@author: PMR
"""

# %% Import libraries

# import external libraries
import numpy as np
import flopy
import flopy.utils.binaryfile as bf
import time as time

# import functions
import f_create_chd as cdic

# import inner modules and variables
import set_sim as ss
import set_waves as sw

# fonts and colors
import matplotlib
font = {'size'   : 10}
matplotlib.rc('font', **font)
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family'] = "sans-serif"


# %% Create and run MODFLOW model

if __name__ == "__main__":
    
    # running time (start)
    start_main = time.time()


    # %% Define settings for modules and packages
    
    # create Time-Invariant Flopy Objects
    modelname = "synth_model"
    mf = flopy.modflow.Modflow(modelname, exe_name="mf2005dbl")
    
    dis = flopy.modflow.ModflowDis(mf, 
                                   ss.nlay, ss.nrow, ss.ncol, 
                                   delr=ss.delr, delc=ss.delc, top=ss.ztop, botm=ss.botm[1:], 
                                   nper=ss.nper, perlen=ss.perlen, nstp=ss.nstp, steady=ss.steady)
    
    bas = flopy.modflow.ModflowBas(mf, ibound=ss.ibound, strt=ss.strt)
    pcg = flopy.modflow.ModflowPcg(mf)
    
    # create Output Control objects
    stress_period_data = {}
    # only_sp = np.arange(0,10,1)
    
    for kper in range(ss.nper):
    # for kper in only_sp:
        for kstp in range(ss.nstp[kper]):
            stress_period_data[(kper, kstp)] = ["save head", "print head"]
            
    oc = flopy.modflow.ModflowOc(mf, stress_period_data=stress_period_data, compact=True)
    
    # create Layer-Property Flow package 
    lpf = flopy.modflow.ModflowLpf(mf, hk=ss.hk, sy=ss.sy, ss=ss.ss, laytyp=ss.laytyp, ipakcb=53)
    
    
    # %% Bring CHD data, write files and run model
    
    for wave in range(len(ss.wave_names)):
        
        for node in range(len(ss.nodes.T)):
            
            
            p_heads = sw.p_heads_all_waves[wave][node]
            
            chd_data = cdic.chd_dictionary(p_heads, ss.nper, ss.nrow)

            chd = flopy.modflow.ModflowChd(mf, stress_period_data=chd_data)
        
            # write the model input files
            mf.write_input()
    
            # run the model
            success, mfoutput = mf.run_model(silent=True, pause=False)
            if not success:
                raise Exception("MODFLOW did not terminate normally.")
            else:
                print('Modflow simulation %d/%d (wave: %s)' % (node+1, len(ss.nodes.T), ss.wave_names[wave]))
             
                
            # %% Extract results
            
            # create the headfile and budget file objects
            headobj = bf.HeadFile(f"{modelname}.hds")
            times = headobj.get_times()
            
            mw_output = []
            for i in range(ss.nper):
                goo = headobj.get_data(totim=times[i])[0,:,:]
                mw_output.append(goo)
            
            foo = np.array(mw_output)
                       
            np.save('./outputs/%s/heads_%d.npy' % (ss.wave_names[wave], node), foo)
            
            
        # %% Close and running time (stop)  
        
        end_main = time.time()
        print('\n')
        print('*' * 50)
        print('Running time: %.4f minutes (wave: %s)' % ((end_main - start_main)/60, ss.wave_names[wave]))
        print('*' * 50)
        print('\n')


    
