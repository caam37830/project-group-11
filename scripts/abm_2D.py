import sys
sys.path.append("../")
from sir.abm import *
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm


pop = new_pop(0, 20, 20)
pop[0,0].pos = (0.5,0.5)
pop[0,0].infect()
pop = pop.flatten()

q = np.sqrt(3 / 10000 / np.pi)
pop, S, I, R = abm_2D_sim(pop, .3, q=q, k=0.01, t=15)