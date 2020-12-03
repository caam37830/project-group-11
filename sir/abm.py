import numpy as np
import random
import matplotlib.pyplot as plt
from tqdm import tqdm
from numba import njit

class Person():
    """
    An agent representing an individual person.

    The state attribute covers the agent's relationship with the disease.
    """

    def __init__(self, state='S'):
        """
        Creates an agent with a default susceptible status
        """

        # Check for valid status if not default
        assert state in {'S', 'I', 'R'}, 'State must be S, I, or R'

        self.state = state

    def cur_state(self):
        """
        Returns the agent's current status
        """

        return self.state

    def __repr__(self):
        return f"Person({self.state})"
    
    def infect(self):
        """
        Agent catches the disease
        Only possible for suscepitble agents
        """
        if self.state == 'S':
            self.state = 'I'
    
    def remove(self):
        """
        Agent recovers and becomes immune (or succumbs to disease).
        Only possible for infected agents
        """

        if self.state == 'I':
            self.state = 'R'

def abm_pop_sim(pop, b, ob, k, t):
    """
    Simulate the spread of a disease on a given population over t days.
    pop designates the starting state, b and k control the spread and recovery rate
    """
    S = []
    I = []
    R = []

    for i in range(t):
        pop = infect_pop(pop, b, ob)
        pop = remove_pop(pop, k)
        S.append(len(get_indices(pop, 'S')))
        I.append(len(get_indices(pop, 'I')))
        R.append(len(get_indices(pop, 'R')))

    return pop, S, I, R


def remove_pop(pop, k):
    """
    Remove k proportion of the infected population 
    """
    nrow, ncol = pop.shape

    # Remove population
    infected = get_indices(pop, 'I')
    
    # Determine the rounded number of infected to remove
    removals = np.int(np.floor(len(infected) * k))

    # Remove k proportion of infected
    for rem in random.sample(infected, removals):
        i = rem // nrow
        j = rem % ncol
        pop[i, j].remove()
    
    return pop


def infect_pop(pop, b, ob):
    """
    Have each infected agent interact with b others
    """
    # Find infected agents
    infected = get_indices(pop, 'I')

    nrow, ncol = pop.shape

    # Infect population 
    for inf in infected:

        # Find grid index of infected
        i = inf // nrow
        j = inf % ncol
        nbrs = neighbors(i, j, nrow, ncol)

        # Identify the b nbr interactions
        # Only use b if there are b neighbors to interact with!
        nei_interactions = min(b, len(nbrs))
        receivers = random.sample(nbrs, nei_interactions)

        # Infect neighbors
        for rec in receivers:
            pop[rec].infect()

        # Check for interacting outside of cohort
        if np.random.rand(1) < ob:
            continue

        # Infect a random other agent
        out_i = np.random.choice(nrow)
        out_j = np.random.choice(ncol)
        pop[out_i,out_j].infect()
        

    return pop

def neighbors(i, j, m, n):
    """
    Returns the neighboring indices for a given index
    in a matrix  of m x n shape
    """

    # Sets default delta indices
    # Adjust delta indices for given index on edge of matrix
    inbrs = [-1, 0, 1]
    if i == 0:
        inbrs = [0, 1]
    if i == m-1:
        inbrs = [-1, 0]
    jnbrs = [-1, 0, 1]
    if j == 0:
        jnbrs = [0, 1]
    if j == n-1:
        jnbrs = [-1, 0]

    nbrs = []
    # Applies deltas and yields neighboring indices
    for delta_i in inbrs:
        for delta_j in jnbrs:
            if delta_i == delta_j == 0:
                continue
            else:
                nbrs.append((i+delta_i, j+delta_j))

    return nbrs

                

def get_indices(pop, state):
    """
    Finds the population agents with a given state.
    Flattens the population for easier index finding.
    """
    # Convert from grid pop to flattened pop
    pop_flat = pop.flatten()

    indices = [i for i, agent in enumerate(pop_flat) if agent.cur_state() == state]
    return indices

def abm_phase(nrow, ncol, infected, t, bs=np.arange(1, 11, dtype=np.int64), ks=np.linspace(0.01, .5, 10), save_path=None):
    """
    plot phase diagram
    :param N: starting population
    :param infected: starting infected count
    :param t: time
    :param bs: discrete b
    :param ks: discrete k
    :param save_path:
    :return:
    """
    # store initial state of pop for future use
    cts = np.zeros((len(bs), len(ks)))
    N = nrow * ncol
    ob = 0

    for i, b in tqdm(enumerate(bs)):
        # ensure b is an int, not a float
        for j, k in enumerate(ks):
            pop = new_pop(infected, nrow, ncol)
            pop, S, I, R = abm_pop_sim(pop, b, ob, k, t)
            cts[i, j] = np.true_divide(I[-1], N)


    # Create a phase plot
    fig, ax = plt.subplots()
    axcontour = ax.contour(ks, bs, cts)
    fig.colorbar(axcontour)
    ax.set_title('phase diagram [I, t={}]'.format(t))
    ax.set_xlabel('k')
    ax.set_ylabel('b')

    # Save plot
    if save_path is not None:
        fig.savefig(save_path)
    plt.show()

def new_pop(start_inf, nrow, ncol):
    """
    Create an initial state population
    :param N: population size
    :param start_inf: count of infected initial state
    """
    # Create a population of susceptibles
    pop = [Person() for i in range(nrow * ncol)]
    pop_size = len(pop)

    infected = random.sample(range(0,pop_size), start_inf)
    # Infect a given number of the population
    for i in infected:
        pop[i].infect()
    
    # Move population into a grid
    pop_grid = np.reshape(pop, (nrow, ncol))
   
    return pop_grid


def time_plot(pop, b, ob, k, t, save_path):
    """
    Runs a simulation with the given parameters and 
    Outputs a plot of the three state ratios over time t
    """
    N = pop.size
    pop, S, I, R = abm_pop_sim(pop, b, ob, k, t)

    s = np.true_divide(S, N)
    i = np.true_divide(I, N)
    r = np.true_divide(R, N)

    # plot pop sim data
    t = np.linspace(1,t,t)
    plt.plot(t, s, label='Susceptible', color='green')
    plt.plot(t, i, label='Infectious', color='red')
    plt.plot(t, r, label='Removed', color='blue')
    plt.title("SIR ABM Simulation, b={},k={}".format(b, k))
    plt.ylabel("ratio")
    plt.xlabel("day")
    plt.legend()
    if save_path is not None:
            plt.savefig(save_path)
    plt.show()
