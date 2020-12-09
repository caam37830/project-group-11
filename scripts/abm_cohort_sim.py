import sys
sys.path.append("../")
from sir.abm import *
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable


# Create a new phase plot to check ob with b
b = 3
obs = np.logspace(-3, 0,15)
ts = np.linspace(0, 50, 11)
cts = np.zeros((len(obs), len(ts)))
for i, ob in tqdm(enumerate(obs)):
    pop = new_pop(1, 100, 100)
    pop, S, I, R = abm_pop_sim(pop, b, ob, .1, 50)
    for j, t in enumerate(ts):
        t = np.int(t)
        cts[i,j] = np.true_divide(I[t-1], 100*100)

plt.figure(figsize=(30,10))
plt.imshow(cts, extent=[np.min(ts), np.max(ts), np.max(obs), np.min(obs)], aspect=5, cmap="Reds")
plt.colorbar()
plt.ylabel('Non-Local Interaction Chance')
plt.yscale('log')
plt.xlabel('# Days')
plt.title('Rate of Infection by Chance of Non-Local Interaction')
plt.savefig('../doc/final/figures/NonLocalPhase.png')
plt.show()



# Make a heatmap
I_list = np.empty(shape=1)
ncol = 100
nrow = 100
num_sims = 100
t = 15

# Simulate a number of epidemics with b=3, ob=.01, k=.01 
for i in tqdm(range(num_sims)):
    pop = new_pop(4, nrow, ncol)
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
plt.savefig('../doc/final/figures/Var1heatmap.png')
plt.show()


# Make a strict cohort heatmap
I_list = np.empty(shape=1)
ncol = 90
nrow = 90
num_sims = 100
t = 15

# Simulate a number of epidemics with b=3, ob=.01, k=.01 
for i in tqdm(range(num_sims)):
    pop = new_pop(4, nrow, ncol)
    pop, S, I, R = abm_pop_sim(pop, 3, .01, .01, t, strict_cohort=True)
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
plt.savefig('../doc/final/figures/StrictCohort Heatmap.png')
plt.show()



# Create a new phase plot to check ob with b this time with strict cohorts
b = 3
obs = np.logspace(-3, 0,15)
ts = np.linspace(0, 50, 11)
cts = np.zeros((len(obs), len(ts)))
for i, ob in tqdm(enumerate(obs)):
    pop = new_pop(1, 90, 90)
    pop, S, I, R = abm_pop_sim(pop, b, ob, .1, 50, strict_cohort=True)
    for j, t in enumerate(ts):
        t = np.int(t)
        cts[i,j] = np.true_divide(I[t-1], 90*90)

plt.figure(figsize=(30,10))
plt.imshow(cts, extent=[np.min(ts), np.max(ts), np.max(obs), np.min(obs)], aspect=5, cmap="Reds")
plt.colorbar()
plt.ylabel('Non-Local Interaction Chance')
plt.yscale('log')
plt.xlabel('# Days')
plt.title('Rate of Infection by Chance of Non-Local Interaction')
plt.savefig('../doc/final/figures/CohortPhase.png')
plt.show()