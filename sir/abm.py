import numpy as np
import random
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



def get_indices(pop, state):
    """
    Finds the population agents with a given state
    """
    indices = [i for i, agent in enumerate(pop) if agent.cur_state() == state]
    return indices