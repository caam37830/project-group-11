import sys
sys.path.append("../")
from sir.abm import *
import matplotlib.pyplot as plt
import numpy as np

<<<<<<< HEAD
# Create 3 time plots at various b, k values
pop = new_pop(10_000, 1)
time_plot(pop, 3, .01, 365, '../doc/checkpoint/figures/abm1.png')
pop = new_pop(10_000, 1)
time_plot(pop, 3, .1, 365, '../doc/checkpoint/figures/abm2.png')
pop = new_pop(10_000, 1)
time_plot(pop, 1, .01, 365, '../doc/checkpoint/figures/abm3.png')

# Code for phase plots

bs=np.arange(1, 11, dtype=np.int64)
ks=np.linspace(0.01, .5,10)

abm_phase(10_000, 1, t=5, ks=ks, bs=bs, save_path='../doc/checkpoint/figures/abm_phase_diagram1.png')
abm_phase(10_000, 1, t=10, ks=ks, bs=bs, save_path='../doc/checkpoint/figures/abm_phase_diagram2.png')
abm_phase(10_000, 1, t=50, ks=ks, bs=bs, save_path='../doc/checkpoint/figures/abm_phase_diagram3.png')
=======
# Create base pop with 0.1% infected rate

def newpop(N, I0):
    """
    Creates a population of size N with I0 initial infected cases
    """
    pop = [Person() for i in range(N)]
    for i in range(I0):
        pop[i].infect()
    return pop

# Pre-compile with a small population to enable jit functionality
pop = newpop(100, 2)
abm_pop_sim(pop, b=3, k=0.1, t=10)

# Set population and disease paramters, get pop sim data
N = 1000
I0 = 10
b = 3
k = 0.1

pop = newpop(N, I0)
pop, S, I, R = abm_pop_sim(pop, b, k, 100)

s = np.true_divide(S, N)
i = np.true_divide(I, N)
r = np.true_divide(R, N)

# plot pop sim data
t = np.linspace(1,100,100)
plt.plot(t, s, label='Susceptible', color='green')
plt.plot(t, i, label='Infectious', color='red')
plt.plot(t, r, label='Removed', color='blue')
plt.title("SIR ABM Simulation, b={},k={}".format(b, k))
plt.ylabel("ratio")
plt.xlabel("day")
plt.legend()
plt.show()



# Code for phase plot
bs = np.arange(1, 10, 0.5)
ks = np.linspace(0.01,.5, 50)

cts = np.zeros((len(bs), len(ks)))
for i, b in enumerate(bs):
    for j, k in enumerate(ks):
        pop = newpop(N, I0)
        pop, S, I, R = abm_pop_sim(pop, np.int(b), k, 25)
        cts[i,j] = np.true_divide(len(I), N)

plt.figure(figsize=(10,5))
plt.imshow(cts, extent=[np.min(ks), np.max(ks), np.max(bs), np.min(bs)])
plt.colorbar()
plt.xlabel('k')
plt.ylabel('b')

plt.show()
>>>>>>> 18bccdcd55913c402ba788f299873f76dcfc5479
