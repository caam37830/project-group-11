import sys
sys.path.append("../")
from sir.abm import abm_phase, time_plot, new_pop
import matplotlib.pyplot as plt
import numpy as np

# Create 3 time plots at various b, k values
pop = new_pop(1, 100, 100)
time_plot(pop, 3, 0, .01, 365, '../doc/checkpoint/figures/abm1.png')
pop = new_pop(1, 100, 100)
time_plot(pop, 3, 0, .1, 365, '../doc/checkpoint/figures/abm2.png')
pop = new_pop(1, 100, 100)
time_plot(pop, 1, 0, .01, 365, '../doc/checkpoint/figures/abm3.png')

# Code for phase plots

bs=np.arange(1, 11, dtype=np.int64)
ks=np.linspace(0.01, .5,25)

abm_phase(100, 100, 1, t=5, ks=ks, bs=bs, save_path='../doc/checkpoint/figures/abm_phase_diagram1.png')
abm_phase(100, 100, 1, t=10, ks=ks, bs=bs, save_path='../doc/checkpoint/figures/abm_phase_diagram2.png')
abm_phase(100, 100, 1, t=50, ks=ks, bs=bs, save_path='../doc/checkpoint/figures/abm_phase_diagram3.png')

