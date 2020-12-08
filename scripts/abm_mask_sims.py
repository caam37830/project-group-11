import os
dir_path = os.getcwd()

from sir.abm import *

nrow, ncol = 100, 100
pop = new_pop(15, nrow, ncol)
time_plot(pop, 3, 1/(nrow*ncol), 0.1, 365, 0, 0, 0.2, 0.7, dir_path+'/doc/final/figures/no_mask1.png')

pop = new_pop(15, nrow, ncol)
time_plot(pop, 3, 1/(nrow*ncol), 0.1, 365, 0.25, 0.05, 0.2, 0.7, dir_path+'/doc/final/figures/mask1.png')

pop = new_pop(200, nrow, ncol)
masked_pop = update_masks(pop, 0.2, 0.7, 0.67, 0)
time_plot(masked_pop, 3, 1/(nrow*ncol), 0.1, 365, 0.25, 0.05, 0.2, 0.7, dir_path+'/doc/final/figures/mask2.png')

pop = new_pop(200, nrow, ncol)
time_plot(pop, 3, 1/(nrow*ncol), 0.1, 365, 0.25, 0.05, 0.2, 0.7, dir_path+'/doc/final/figures/no_mask2.png')

pop = new_pop(15, nrow, ncol)
time_plot(pop, 3, 1/(nrow*ncol), 0.1, 365, 0.01, 1/365, 1, 0.05, dir_path+'/doc/final/figures/vaccine.png')

abm_mask_phase(np.linspace(0,1,20), np.linspace(0,1,20), 0.2, 0.7, save_path = dir_path+'doc/final/figures/mask_phase1.png',
nrow = 100, ncol = 100)

abm_mask_phase(0.2, 0.05, np.linspace(0,1,50), np.linspace(0,1,50), save_path=dir_path+'doc/final/figures/mask_phase2.png',
nrow=100, ncol=100)