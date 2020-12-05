import sys
sys.path.append("../")
from sir.abm import *
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable


# Create 3 new time plots with outside interaction chance.
pop = new_pop(1, 100, 100)
time_plot(pop, 3, 0.01, .01, 50,'../doc/checkpoint/figures/abm3.png')
pop = new_pop(1, 100, 100)
time_plot(pop, 3, 0.01, .1, 50,'../doc/checkpoint/figures/abm3.png')
pop = new_pop(1, 100, 100)
time_plot(pop, 1, 0.01, .01, 50,'../doc/checkpoint/figures/abm3.png')

# Create 3 new phase plots to check ob with b
bs=np.arange(1, 11, dtype=np.int64)
obs=np.logspace(-3, 0,10)

abm_phase(100, 100, 1, t=5, bs=bs, obs=obs, obs_phase=True, save_path='../doc/checkpoint/figures/abm_phase_diagram1.png')
abm_phase(100, 100, 1, t=10, bs=bs, obs=obs, obs_phase=True, save_path='../doc/checkpoint/figures/abm_phase_diagram2.png')
abm_phase(100, 100, 1, t=50, bs=bs, obs=obs, obs_phase=True, save_path='../doc/checkpoint/figures/abm_phase_diagram3.png')


# Make a heatmap
I_list = np.empty(shape=1)
ncol = 90
nrow = 90
num_sims = 50
t = 20


# Simulate a number of epidemics with b=3, ob=.01, k=.01 
for i in tqdm(range(num_sims)):
    pop = new_pop(1, nrow, ncol)
    pop, S, I, R = abm_pop_sim(pop, 3, .01, .01, t)
    pop_f = pop.flatten()
    
    # Record the agents that remain susceptible at day t
    I_list = np.append(I_list, get_indices(pop_f, 'S'))

# Calculate the average rate of infection per agent
infected_rates = np.empty(shape=(ncol*nrow))
for i in range(nrow*ncol):
    suscept_instances = np.array(np.where(I_list == i))
    susceptible_rate = suscept_instances.shape[1] / num_sims
    infected_rates[i] = 1 - susceptible_rate

# Plot an infection heat map
infected_rates = np.reshape(infected_rates, (nrow, ncol))
ax = plt.subplot(111)
plt.title('Infection Rate Heatmap(b=3,ob=.01,k=.01)')
im = ax.imshow(infected_rates, cmap='gist_heat')

# Add a colorbar
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
plt.savefig('../doc/checkpoint/figures/Var1heatmap.png')
plt.show()

# Make a heatmap
I_list = np.empty(shape=1)
ncol = 90
nrow = 90
num_sims = 100
t = 20

# Simulate a number of epidemics with b=3, ob=.01, k=.01 
for i in tqdm(range(num_sims)):
    pop = new_pop(1, nrow, ncol)
    pop, S, I, R = abm_pop_sim(pop, 3, .05, .01, t, strict_cohort=True, cohort_size=9)

    # Record the agents that remain susceptible at day t
    I_list = np.append(I_list, get_indices(pop, 'S'))

# Calculate the average rate of infection per agent
infected_rates = np.empty(shape=(ncol*nrow))
for i in range(nrow*ncol):
    suscept_instances = np.array(np.where(I_list == i))
    susceptible_rate = suscept_instances.shape[1] / num_sims
    infected_rates[i] = 1 - susceptible_rate

# Plot an infection heat map
infected_rates = np.reshape(infected_rates, (nrow, ncol))
ax = plt.subplot(111)
plt.title('Strict Cohort Infection Rate Heatmap(b=3,ob=.01,k=.01)')
im = ax.imshow(infected_rates, cmap='gist_heat')

# Add a colorbar
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)
plt.savefig('../doc/checkpoint/figures/Var1heatmap.png')
plt.show()
