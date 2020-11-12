from sir.abm import *
import matplotlib.pyplot as plt
import numpy as np

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