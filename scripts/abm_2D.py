import sys
sys.path.append("../")
from sir.abm import *
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

q = np.sqrt(3 / (50*50) / np.pi)
k = .1

# Compare step value impact on infection rates
ps = [.02, 0.05, .1, .25, .4, .6, .75]
cts = np.zeros((1, len(ps)))
for j, p in enumerate(ps):
    pop = new_pop(0, 50, 50)
    for i in range(2):
        pop[i, 0].infect()
        pop[i, 0].pos = (0.5,0.5)
    time_plot(pop, 0, 0, .1, 100 , 0, 0, 0, 0, '../doc/final/figures/abmStepSize%.2f.png' % p, q=q, p=p)



t = 100
p = .02

# Centered infected start
pop = new_pop(0,50,50)
for i in range(4):
    pop[i, 0].infect()
    pop[i, 0].pos = (0.5,0.5)
time_plot(pop, 0, 0, .1, t, 0, 0, 0, 0, '../doc/final/figures/abm2dCenter.png', q=q, p=p)

# Corner infected start
pop = new_pop(0,50,50)
for i in range(4):
    pop[i, 0].infect()
    pop[i, 0].pos = (0,0)
time_plot(pop, 0, 0, .1, t, 0, 0, 0, 0, '../doc/final/figures/abm2dCorner.png', q=q, p=p)

# Random start
pop = new_pop(4, 50, 50)
time_plot(pop, 0, 0, .1, t, 0, 0, 0, 0, '../doc/final/figures/abm2dRandom.png', q=q, p=p)
