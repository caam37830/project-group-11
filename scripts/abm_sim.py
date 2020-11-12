from sir.abm import *
import matplotlib.pyplot as plt
import numpy as np

# Create base pop with 0.1% infected rate
N = 100_000
def newpop():
    pop = [Person() for i in range(N)]
    for i in range(100):
        pop[i].infect()
    return pop

# Set b and k variables, get pop sim data
b = 3
k = 0.1

pop = newpop()
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
bs = np.linspace(1, 49, 50)
ks = np.linspace(0.01,.5, 10)

cts = np.zeros((len(bs), len(ks)))
for i, b in enumerate(bs):
    for j, k in enumerate(ks):
        pop = newpop()
        pop, S, I, R = abm_pop_sim(pop, np.int(b), k, 25)
        cts[i,j] = np.true_divide(len(I), N)

plt.figure(figsize=(10,5))
plt.imshow(cts, extent=[np.min(ks), np.max(ks), np.max(bs), np.min(bs)])
plt.colorbar()
plt.xlabel('k')
plt.ylabel('b')

plt.show()