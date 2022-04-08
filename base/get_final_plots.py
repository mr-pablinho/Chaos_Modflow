# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 17:07:10 2021
FINAL PLOTS - STOCHASTIC WAVES
@author: PMR
"""

# %% Import libraries

# import external libraries
import numpy as np
import chaospy as cp
import matplotlib.pyplot as plt

# fonts and colors
import matplotlib
font = {'size'   : 10}
matplotlib.rc('font', **font)
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family'] = "sans-serif"

# import inner modules and variables
import f_flow_field_eval as fpf
import f_wave_functions as fw
from set_sim import t, nrow, ncol, wave_names, nodes, weights, J_distro, polynomial_expansion, delr, delc, hk, duration


# %% Import results

heads_all = np.zeros((int(len(wave_names)), int(len(nodes.T)), int(len(t)), int(nrow), int(ncol)))

# create list for heads: [wave][node][time,row,col]
for wave in range(len(wave_names)):
    for node in range(len(nodes.T)):
        heads_all[wave,node,:,:,:] = np.load('./outputs/%s/heads_%d.npy' % (wave_names[wave], node))
        print('Gathering groundwater heads: wave=%s, node=%d' % (wave_names[wave], node))



# %% 

# times = [1223]
# # times = np.arange(1220,1225,1)

# okw_wave = [None] * len(wave_names)
# vor_wave = [None] * len(wave_names)
# she_wave = [None] * len(wave_names)
# str_wave = [None] * len(wave_names)

# for wave in range(len(wave_names)):
#     print(wave)
#     for step in times:
#         okw_wave[wave] = np.load('./outputs/uq_topo/evals_okw.npy')[wave,:,step,:,:]


        
# # %% Compute uncertainty

# E_gwh = np.zeros((len(wave_names), len(t), nrow, ncol))
# S_gwh = np.zeros((len(wave_names), len(t), nrow, ncol))
f_approx_gwh = [None] * len(wave_names)
f_approx_okw = [None] * len(wave_names)



# %% Get expansion in one cell at all times

# choose specific cell for the analyisis
cell_rows = [50]
cell_cols = [5]

for wave in range(len(wave_names)):
        for row in cell_rows:
            print('Computing polynomial expansions (wave:%s, row:%d)' % (wave_names[wave], row+1))
            for col in cell_cols:                
                cell_eval_gwh = heads_all[wave, :, :, row, col]
                f_approx_gwh[wave] = cp.fit_quadrature(polynomial_expansion, nodes, weights, cell_eval_gwh)

                       

# %% Plot: Groundwater heads uncertainty PDFs 

x = np.arange(-2, 2, 0.01)
time_plot = np.arange(0,20,5)

for ti in time_plot:
    fig = plt.figure(ti, figsize=[4,2.5])
    plt.title('Time step: ' + str(ti), loc='left')
    for wave in range(len(wave_names)):
        qoi_dist = cp.QoI_Dist(f_approx_gwh[wave][ti], J_distro, sample=10000)
        
        a = qoi_dist.pdf(x)
        a = a / np.sum(a)
        
        plt.plot(x, a, alpha=0.99, label=wave_names[wave], linewidth=1.5)
        
    plt.legend()
    plt.xlim(-2,2)
    plt.ylim(0, 0.1)
    # plt.ylim(0, 6)
    plt.xlabel('Groundwater head')
    plt.ylabel('Frequency')
    # plt.savefig('./results/uq_heads/pdf_ts%s.png' % (ti))
    # plt.close(fig) 



# %% Plot: Hydraulic conductivity field

# import final hyd. conductivity field
plt.figure('Hyd. field')
plt.imshow(np.log10(hk), cmap='viridis')
plt.colorbar()
# plt.clim(-4.5,-1.5)



# %% Plot: Simple waves

from fractions import Fraction

# amplitude
A = 1
# frequency
p = 8  # in hours
f = 1 / (p * 3600)  # in cycles per second

triangular_phase = p * 3600 / 4

# waves = np.load('./inputs/simple_waves.npy')
waves = [fw.sine_wave(f, t, A, phase=0),
          fw.trapezoid_wave(f, t, A, phase=0, ext=6),
          fw.triangular_wave(f, t, A, phase=0),
          fw.sine_wave(f, t, A, phase=0) + fw.sine_wave(f*3, t, A=0.4*A, phase=0)
          ]
waves = np.array(waves)
t_coord = np.arange(0,len(t),1)


x_ticks = np.linspace(0, 480, 480)


X_tick = np.array([])
xx = np.linspace(0, 9, int((8*9)/8), endpoint=0)
# xx = np.hstack((xx, 9))

for item in xx:
    X_tick = np.append(X_tick,Fraction(item).limit_denominator())
# plt.xticks(np.unique(xx), X_tick, fontsize=9)
# plt.xticks(xx, fontsize=9)


# import final hyd. conductivity field
fig, axs = plt.subplots(4, 1)
notation_wave = ['s', 'z', 't', 'c']
wave_title = ['Sine wave', 'Trapezoid wave', 'Triangle wave', 'Complex wave']
for i in range(4):
    axs[i].plot(x_ticks, waves[i,:480], color='red')
    axs[i].set_xlim(0,430)
    axs[i].set_ylim(-1.5,1.5)
    axs[i].set_yticks([-1,0,1])
    axs[i].set_xticks(np.arange(0,480,80))
    
    axs[i].set_title(wave_title[i], loc='left')
    
    axs[i].spines['right'].set_visible(False)
    axs[i].spines['top'].set_visible(False)
    axs[i].yaxis.set_ticks_position('left')
    axs[i].xaxis.set_ticks_position('bottom')
    
    axs[i].set_ylabel(r'$\mathrm{y_\widetilde %s}$' % (notation_wave[i]))
    
    axs[i].plot((1), (-1.5), ls="", marker=">", ms=5, color="k",
            transform=axs[i].get_yaxis_transform(), clip_on=False)
    axs[i].hlines(xmin=0, xmax=480, y=0, 
                  linestyle=':', alpha=.99, color='black', linewidth=1.0)

axs[-1].set_xlabel(r'$\mathrm{\tau/T}$')

fig.tight_layout()
fig.set_figheight(5)
fig.set_figwidth(3)



# %% Plot: UQ in the groundwater heads 2D

# import external libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# import heads
E_gwh = np.zeros((len(wave_names), len(t), nrow, ncol))
S_gwh = np.zeros((len(wave_names), len(t), nrow, ncol))

for wave in range(len(wave_names)):
    E_gwh[wave,:,:,:] = np.load('./outputs/uq_heads/E_heads__%s.npy' % (wave_names[wave]))
    S_gwh[wave,:,:,:] = np.load('./outputs/uq_heads/S_heads__%s.npy' % (wave_names[wave]))


titles_wave = ['Sine waves',
               'Trapzoid waves',
               'Triangle waves',
               'Composed waves']

# display of the subplots
rowfig = 2
colfig = 3

# choose stress periods
time_eval = [65,115]
# time_eval = np.arange(0,30,3)

# set min and max for the colorbar
vmin_e, vmax_e = -1, 1
vmin_s, vmax_s = 0, 1

for ti in time_eval:
    fig, ax = plt.subplots(nrows=rowfig, ncols=colfig, figsize=[7,7], 
                            gridspec_kw={'width_ratios':[1,1,1],'height_ratios': [1,1]},
                            constrained_layout=True)
    
    fig.suptitle('SP: %d' % (ti), fontsize=20, horizontalalignment='left')
    
    for wave in range(len(wave_names)):
        mat_plot = E_gwh[wave][ti][:, :61]
        
        # correction for "-inf" and "nan"
        mat_plot = np.where(~np.isnan(mat_plot), mat_plot, -10)
        mat_plot = np.where(np.isfinite(mat_plot), mat_plot, -10)
        
        ime = ax[0][wave].imshow(mat_plot, 
                                 vmin=vmin_e, vmax=vmax_e, 
                                 cmap='Spectral')
        ax[0][wave].set_title(titles_wave[wave])
    ax[0][0].set_ylabel(r'$\mu$')
     
    # position for the colorbar (mu)
    cbaxes_e = fig.add_axes([1.04, 0.50, 0.025, 0.405]) 
    cb = plt.colorbar(ime, cax=cbaxes_e)
    
    
    
    for wave in range(len(wave_names)):
        ims = ax[1][wave].imshow(S_gwh[wave][ti][:, :61], 
                                 vmin=vmin_s, vmax=vmax_s, 
                                 cmap='gnuplot2')
    
    ax[1][0].set_ylabel(r'$\sigma$')
    
    # position for the colorbar (sigma)
    cbaxes_s = fig.add_axes([1.04, 0.04, 0.025, 0.405]) 
    cb = plt.colorbar(ims, cax=cbaxes_s)
    
    # plt.savefig('./images/heads/2d_gwh_uq_%d.png' % (ti), dpi=300, bbox_inches='tight')
    # plt.close(fig) 



# %% Plot: UQ in the Okubo-Weiss 2D

# import external libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from f_colors import cmap_pmr_4, cmap_pmr_5, cmap_pmr_6, cmap_pmr_div_5, cmap_pmr_div_6


E_okw = np.zeros((len(wave_names), len(t), nrow, ncol))
S_okw = np.zeros((len(wave_names), len(t), nrow, ncol))

# import Okubo-Weiss 2D
for wave in range(len(wave_names)):
    E_okw[wave,:,:,:] = np.load('./outputs/uq_topo/E_okw__%s.npy' % (wave_names[wave]))
    S_okw[wave,:,:,:] = np.load('./outputs/uq_topo/S_okw__%s.npy' % (wave_names[wave]))


titles_wave = ['Sine waves',
               'Trapezoid waves',
               'Triangle waves',
               'Composed waves']

# display of the subplots
rowfig = 2
colfig = 3

# choose stress periods
time_eval = [65]
# time_eval = np.arange(0,30,3)

# set min and max for the colorbar
vmin_e, vmax_e = -1e-6, 1e-6
vmin_s, vmax_s = 0, 1e-5

for ti in time_eval:
    fig, ax = plt.subplots(nrows=rowfig, ncols=colfig, figsize=[7,7], 
                            gridspec_kw={'width_ratios':[1,1,1],'height_ratios': [1,1]},
                            constrained_layout=True)
    
    fig.suptitle('SP: %d' % (ti), fontsize=20, horizontalalignment='left')
    
    for wave in range(len(wave_names)):
        mat_plot = E_okw[wave][ti][:, :61]
        
        # correction for "-inf" and "nan"
        mat_plot = np.where(~np.isnan(mat_plot), mat_plot, -10)
        mat_plot = np.where(np.isfinite(mat_plot), mat_plot, -10)
        
        ime = ax[0][wave].imshow(mat_plot, 
                                 vmin=vmin_e, vmax=vmax_e, 
                                 cmap='Spectral')
        ax[0][wave].set_title(titles_wave[wave])
    ax[0][0].set_ylabel(r'$\mu$')
     
    # position for the colorbar (mu)
    cbaxes_e = fig.add_axes([1.04, 0.50, 0.025, 0.405]) 
    cb = plt.colorbar(ime, cax=cbaxes_e)
    
    
    
    for wave in range(len(wave_names)):
        ims = ax[1][wave].imshow((S_okw[wave][ti][:, :61]), 
                                   vmin=vmin_s, vmax=vmax_s, 
                                 cmap='gnuplot2',  # cmap_pmr_5
                                 )
        
    
    ax[1][0].set_ylabel(r'$\sigma$')
    
    # position for the colorbar (sigma)
    cbaxes_s = fig.add_axes([1.04, 0.04, 0.025, 0.405]) 
    cb = plt.colorbar(ims, cax=cbaxes_s)
    
    plt.savefig('./images/topology/2d_okw_uq_%d.png' % (ti), dpi=300, bbox_inches='tight')
    # plt.close(fig) 



# %% Plot: UQ in the groundwater heads

# import heads
E_gwh = np.zeros((len(wave_names), len(t), nrow, ncol))
S_gwh = np.zeros((len(wave_names), len(t), nrow, ncol))

for wave in range(len(wave_names)):
    E_gwh[wave,:,:,:] = np.load('./outputs/uq_heads/E_heads__%s.npy' % (wave_names[wave]))
    S_gwh[wave,:,:,:] = np.load('./outputs/uq_heads/S_heads__%s.npy' % (wave_names[wave]))

# choose stress periods
burn = 0
up_to = 720
t_coord = np.arange(burn, up_to, 1)

fig, axs = plt.subplots(3, 3)

choose_row = 10
choose_col = 5

for wave in range(3):
    
    lineE = np.mean(E_gwh[wave, burn:up_to, :, choose_col], axis=1)

    axs[0,wave].plot(t_coord, lineE,
                   linewidth=1,
                   color='black')

    axs[0,wave].set_xlim(burn, up_to)
    axs[0,wave].set_ylim(-1.5,1.5)
    axs[0,wave].set_yticks([-1,0,1])
    # axs[0,wave].set_xticks([0,400,800,1200,1600])
    
    lineS = np.mean(S_gwh[wave, burn:up_to, :, choose_col], axis=1)

    axs[1,wave].plot(t_coord, lineS,
                   linewidth=1,
                   color='red')
    axs[1,wave].set_xlim(burn, up_to)
    axs[1,wave].set_ylim(0,0.9)
    # axs[wave].set_yticks([3,5,7])
    # axs[1,wave].set_xticks([0,400,800,1200,1600])
    
    axs[2,wave].fill_between(t_coord, lineE - lineS, lineE + lineS,
                   color='red', alpha=0.5, edgecolor=None)
    axs[2,wave].fill_between(t_coord, lineE - 2*lineS, lineE + 2*lineS,
                   color='red', alpha=0.25, edgecolor=None)
    axs[2,wave].plot(t_coord, lineE,
                   linewidth=1,
                   color='black')
    axs[2,wave].set_xlim(burn, up_to)
    axs[2,wave].set_ylim(-2,2)
    # axs[wave].set_yticks([3,5,7])
    # axs[1,wave].set_xticks([0,400,800,1200,1600])

fig.tight_layout()
fig.set_figheight(10)
fig.set_figwidth(10)



# %% Plot: UQ in the Okubo-Weiss 

fp = f * 3600  # in hours
t_c = np.linspace(0, (len(t)/240*3), len(t))
tau = 1/fp
dp = (duration * 24)  # in hours
x_ticks = np.linspace(0, dp/tau, len(t))


# import Okubo-Weiss results
E_okw = np.zeros((len(wave_names), len(t), nrow, ncol))
S_okw = np.zeros((len(wave_names), len(t), nrow, ncol))

# import Okubo-Weiss 2D
for wave in range(len(wave_names)):
    E_okw[wave,:,:,:] = np.load('./outputs/uq_topo/E_okw__%s.npy' % (wave_names[wave]))
    S_okw[wave,:,:,:] = np.load('./outputs/uq_topo/S_okw__%s.npy' % (wave_names[wave]))

# choose stress periods
burn = 0
n_tau = 3
up_to = 720
to_show = int(720/n_tau)
t_coord = np.arange(burn,up_to,1)

fig, axs = plt.subplots(3, 3)

choose_row = 10
choose_col = 5

for wave in range(3):
    
    lineE = np.mean(E_okw[wave, burn:up_to, :, choose_col], axis=1)
    # lineE = E_okw[wave, burn:up_to, choose_row, choose_col]

    axs[0,wave].plot(x_ticks, lineE,
                   linewidth=0.75,
                   color='black')

    axs[0,wave].set_xlim(0, n_tau)
    # axs[wave].set_ylim(3,7)
    # axs[wave].set_yticks([3,5,7])
    # axs[0,wave].set_xticks([400,800,1200,1600])
    
    lineS = np.mean(S_okw[wave, burn:up_to, :, choose_col], axis=1)
    # lineS = S_okw[wave, burn:up_to, choose_row, choose_col]

    axs[1,wave].plot(x_ticks, lineS,
                   linewidth=0.75,
                   color='red')
    axs[1,wave].set_xlim(0, n_tau)
    # axs[wave].set_ylim(3,7)
    # axs[wave].set_yticks([3,5,7])
    # axs[1,wave].set_xticks([400,800,1200,1600])
    
    axs[2,wave].fill_between(x_ticks, lineE - lineS, lineE + lineS,
                   color='red', alpha=0.5, edgecolor=None)
    axs[2,wave].fill_between(x_ticks, lineE - 2*lineS, lineE + 2*lineS,
                   color='red', alpha=0.0, edgecolor=None)
    axs[2,wave].plot(x_ticks, lineE,
                   linewidth=1,
                   color='black')
    axs[2,wave].set_xlim(0, n_tau)
    # axs[2,wave].set_ylim(0,1.5e-8)
    
    '''
    Me quedé aquí. Quiero incluir el plot de las sombras en el okubo-weiss
    axs[2,wave].fill_between(t_coord, lineE - lineS, lineE + lineS,
                   color='red', alpha=0.5, edgecolor=None)
    axs[2,wave].fill_between(t_coord, lineE - 2*lineS, lineE + 2*lineS,
                   color='red', alpha=0.25, edgecolor=None)
    axs[2,wave].plot(t_coord, lineE,
                   linewidth=1,
                   color='black')
    axs[2,wave].set_xlim(burn, up_to)
    axs[2,wave].set_ylim(-2,2)
    '''

fig.tight_layout()
fig.set_figheight(6)
fig.set_figwidth(10)


