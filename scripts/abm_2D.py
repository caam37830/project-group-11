import sys
sys.path.append("../")
from sir.abm import *
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm



q = np.sqrt(3 / 10000 / np.pi)
p = .3
num_sims = 3
t = 15
center = np.zeros(t)
corners = np.zeros(t)
randoms = np.zeros(t)

for i in tqdm(range(num_sims)):
    # Centered infected start
    pop = new_pop(0,50,50)
    for i in range(4):
        pop[i, 0].infect()
        pop[i, 0].pos = (0.5,0.5)
    pop, S, I, R = abm_2D_sim(pop, .3, q=q, k=0.1, t=15)
    center = np.add(center, np.true_divide(S, 50*50))

    # Corner infected start
    pop = new_pop(0,50,50)
    for i in range(4):
        pop[i, 0].infect()
        pop[i, 0].pos = (0,0)
    pop, S, I, R = abm_2D_sim(pop, .3, q=q, k=0.1, t=15)
    corners = np.add(corners, np.true_divide(S, 50*50))

    pop = new_pop(4, 50, 50)
    pop, S, I, R = abm_2D_sim(pop, .3, q=q, k=0.1, t=15)
    randoms = np.add(randoms, np.true_divide(S, 50*50))


t = np.linspace(1,t,t)
plt.plot(t, corners/num_sims, label='Center', color='gold')
plt.plot(t, center/num_sims, label='Corner', color='purple')
plt.plot(t, randoms/num_sims, label='Random', color='black')
plt.legend('Infection Start')
plt.title("Infection Progress in Spatial Model, q={}, p={}".format(q, p))
plt.ylabel("ratio")
plt.xlabel("day")

plt.savefig('../doc/checkpoint/figures/spatialComp.png')