from sir.abm import *
import matplotlib.pyplot as plt
import numpy as np

# Create base pop with 0.1% infected rate
N = 10_000
pop = [Person() for i in range(N)]
for i in range(10):
    pop[i].infect()

# Set b and k variables, get pop sim data
b = 5
k = 0.1
pop, S, I, R = abm_pop_sim(pop, 2, .1, 50)

# plot pop sim data
t = np.linspace(1,50,50)
plt.plot(t, S, label='Susceptible', color='green')
plt.plot(t, I, label='Infected', color='red')
plt.plot(t, R, label='Removed', color='blue')
plt.title("SIR ABM Simulation, b={},k={}".format(b, k))
plt.ylabel("count")
plt.xlabel("day")
plt.legend()
plt.show()
