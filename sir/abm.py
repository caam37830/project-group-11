import numpy as np
import random
import matplotlib.pyplot as plt
from tqdm import tqdm
from numba import jit

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

@jit
def abm_pop_sim(pop, b, k, t):
    """
    Simulate the spread of a disease on a given population over t days.
    pop designates the starting state, b and k control the spread and recovery rate
    """
    S = []
    I = []
    R = []

    for i in range(t):
        pop = remove_pop(pop, k)
        pop = infect_pop(pop, b)
        S.append(len(get_indices(pop, 'S')))
        I.append(len(get_indices(pop, 'I')))
        R.append(len(get_indices(pop, 'R')))

    return pop, S, I, R

@jit
def remove_pop(pop, k):
    """
    Remove k proportion of the infected population 
    """
    # Remove population
    infected = get_indices(pop, 'I')

    # Determine the rounded number of infected to remove
    removals = np.int(np.floor(len(infected) * k))

    # Remove k proportion of infected
    for i in random.sample(infected, removals):
        pop[i].remove()
    
    return pop

@jit
def infect_pop(pop, b):
    """
    Have each infected agent interact with b others
    """
    # Find infected agents
    infected = get_indices(pop, 'I')
    pop_size = len(pop)

    # Infect population 
    for i in infected:

        # Identify the b interactions for this agent
        receivers = random.sample(range(0,pop_size), b)

        for j in receivers:
            pop[j].infect()

    return pop


def get_indices(pop, state):
    """
    Finds the population agents with a given state
    """
    indices = [i for i, agent in enumerate(pop) if agent.cur_state() == state]
    return indices

def abm_phase(N, infected, t, bs=np.arange(1, 11, dtype=np.int64), ks=np.linspace(0.01, .5, 10), save_path=None):
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

    for i, b in tqdm(enumerate(bs)):
        # ensure b is an int, not a float
        for j, k in enumerate(ks):
            pop = new_pop(N, infected)
            pop, S, I, R = abm_pop_sim(pop, b, k, t)
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

def new_pop(N, infected):
    """
    Create an initial state population
    :param N: population size
    :param infected: count of infected initial state
    """
    # Create a population of susceptibles
    pop = [Person() for i in range(N)]

    # Infect a given number of the population
    for i in range(infected):
        pop[i].infect()
    return pop


def time_plot(pop, b, k, t, save_path):
    """
    Runs a simulation with the given parameters and 
    Outputs a plot of the three state ratios over time t
    """
    N = len(pop)
    pop, S, I, R = abm_pop_sim(pop, b, k, t)

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
