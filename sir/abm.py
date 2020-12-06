import numpy as np
import random
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.spatial import KDTree

class Person():
    """
    An agent representing an individual person.

    The state attribute covers the agent's relationship with the disease.
    """

    def __init__(self, state='S', cohort='No Cohort', pos=np.random.uniform(0,1,2)):
        """
        Creates an agent with a default susceptible status
        """

        # Check for valid status if not default
        assert state in {'S', 'I', 'R'}, 'State must be S, I, or R'

        self.state = state
        self.cohort = cohort
        self.pos = pos

        # Start with no mask
        self.masked = False

    def cur_state(self):
        """
        Returns the agent's current status
        """

        return self.state

    def __repr__(self):
        return f"Person({self.state}, C: {self.cohort})"
    
    def infect(self):
        """
        Agent catches the disease
        Only possible for susceptible agents
        """
        if self.state == 'S':
            if self.masked and np.random.rand(1) < self.risk or not self.masked:
                self.state = 'I'

    def remove(self):
        """
        Agent recovers and becomes immune (or succumbs to disease).
        Only possible for infected agents
        """

        if self.state == 'I':
            self.state = 'R'
    
    def group(self, cohort):
        """
        Assign agent to a cohort
        """
        self.cohort = cohort
    
    def move(self, p):
        """
        Move an agent by p distance
        """
        dpos = np.random.randn(2)
        dpos = dpos / np.linalg.norm(dpos)

        new = self.pos + (dpos * p)
        if ((new[0] >= 0 and new[0] <= 1) and (new[1] >= 0 and new[1] <= 1)):
            self.pos = new

    def mask(self, infectivity, risk):
        """
        Grants an individual the mask attribute which can be used to modify their infectivity/ risk of infection depending on state
        """
        self.masked = True
        self.infectivity = infectivity
        self.risk = risk

    def unmask(self):
        """
        Change an individual's mask status to False
        """
        if self.masked == True:
            self.masked = False




def abm_pop_sim(pop, b, ob, k, t, strict_cohort=False, cohort_size=9):
    """
    Simulate the spread of a disease on a given population over t days.
    pop designates the starting state, b and k control the spread and recovery rate
    """
    S = []
    I = []
    R = []

    if strict_cohort == True:
        pop = make_cohorts(pop, size=cohort_size)

    for i in range(t):
        pop = infect_pop(pop, b, ob, strict_cohort)
        pop = remove_pop_grid(pop, k)
        pop_flat = pop.flatten()
        S.append(len(get_indices(pop_flat, 'S')))
        I.append(len(get_indices(pop_flat, 'I')))
        R.append(len(get_indices(pop_flat, 'R')))

    return pop, S, I, R


def remove_pop_grid(pop, k):
    """
    Remove k proportion of the infected population 
    Uses a grid population
    """
    nrow, ncol = pop.shape

    # Find current infected population
    infected = get_indices(pop.flatten(), 'I')
    
    # Determine the rounded number of infected to remove
    removals = np.int(np.floor(len(infected) * k))

    # Remove k proportion of infected
    for rem in random.sample(infected, removals):
        i = rem // nrow
        j = rem % ncol
        pop[i, j].remove()
    
    return pop


def remove_pop(pop, k):
    """
    Remove k proportion of the infected population 
    Uses a flattened population
    """
    # Find current infected population
    infected = get_indices(pop, 'I')
    
    # Determine the rounded number of infected to remove
    removals = np.int(np.floor(len(infected) * k))

    # Remove k proportion of infected
    for rem in random.sample(infected, removals):
        pop[rem].remove()
    return pop

def update_masks(pop, infectivity, risk, m, u):
    """
    Grants masks to m proportion of the unmasked susceptible and infected population
    Take away masks from u proportion of the masked susceptible and infected population    
    """
    sus_m = get_indices(pop.flatten(), 'S', True)
    inf_m = get_indices(pop.flatten(), 'I', True)
    sus_u = get_indices(pop.flatten(), 'S', False)
    inf_u = get_indices(pop.flatten(), 'I', False)
    masks = sus_u + inf_u
    unmasks = sus_m + inf_m
    mask_num = np.int(np.floor(len(masks) * m))
    unmask_num = np.int(np.floor(len(unmasks)*u))

    for unmasked in random.sample(masks, mask_num):
        pop.flatten()[unmasked].mask(infectivity, risk)

    for masked in random.sample(unmasks, unmask_num):
        pop.flatten()[masked].unmask()
    
    return pop

def infect_pop(pop, b, ob, strict_cohort):
    """
    Have each infected agent interact with b of their neighbors.
    Each infected agent also has a chance 'ob' to interact with any other agent.
    """
    # Find flattened index of infected agents
    
    pop_flat = pop.flatten()
    infected = get_indices(pop_flat, 'I')
    nrow, ncol = pop.shape
    
    # Infect population 
    for inf in infected:

        # Get indices on 2D grid
        i = inf // ncol
        j = inf % nrow

        # Determine if one of the interactions will be out of cohort
        if np.random.rand(1) < ob:
            outside_int = True
            inside_int = b-1
        else:
            outside_int = False
            inside_int = b
        
        # Strict Cohort Policy -- Find members of static cohort
        if strict_cohort == True:
            # Identify infected agent's cohort and other members of this cohort
            inf_cohort = pop_flat[inf].cohort
            cohort_mems = [i for i, agent in enumerate(pop_flat) if agent.cohort == inf_cohort]
            
            # Identify the b cohort interactions
            # If fewer than b cohort members, just do all members
            coh_interactions = min(inside_int, len(cohort_mems))
            receivers = random.sample(cohort_mems, coh_interactions)
            
            for rec in receivers:
                if pop[i, j].masked and np.random.rand(1) < pop[i, j].infectivity or not pop_flat[inf].masked:
                    pop_flat[rec].infect()
        # Non-strict Cohort Policy -- Find dynamic neighbors
        else:
            # Find neighbors of current agent
            nbrs = neighbors(i, j, nrow, ncol)

            # Identify the b nbr interactions
            # If fewer than b neighbors, just do all neighbors
            nei_interactions = min(inside_int, len(nbrs))
            receivers = random.sample(nbrs, nei_interactions)

            # Infect neighbors
            for rec in receivers:
                if pop[i, j].masked and np.random.rand(1) < pop[i, j].risk or not pop[i, j].masked:
                    pop[rec].infect()

        # Interact outside of cohort if check passed
        if outside_int:
            # Infect a random other agent
            out_i = np.random.choice(nrow)
            out_j = np.random.choice(ncol)
            other = pop[out_i,out_j] 
            
            if pop[i, j].masked and np.random.rand(1) < pop[i, j].infectivity or not pop[i, j].masked:
                other.infect()

    return pop

def neighbors(i, j, m, n):
    """
    Returns the neighboring indices for a given index
    in a matrix of m x n shape
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

            # Ignore 0,0 neighbor (the agent itself)
            if delta_i == delta_j == 0:
                continue
            else:
                nbrs.append((i+delta_i, j+delta_j))

    return nbrs


def get_indices(pop, state, masked=None):
    """
    Finds the population agents with a given state.
    Flattens the population for easier index finding.
    """
    # use flattened pop
    pop = pop.flatten()
    if masked is None:
        indices = [i for i, agent in enumerate(pop) if agent.cur_state() == state]
    else:
        assert isinstance(masked, bool), 'masked input must be a boolean'
        indices = [i for i, agent in enumerate(pop) if agent.cur_state() == state and agent.masked == masked]
    return indices

def abm_phase(nrow, ncol, infected, t, bs=np.arange(1, 11, dtype=np.int64), ks=np.linspace(0.01, .5, 10), save_path=None, obs=None, obs_phase=False):
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

    for i, b in tqdm(enumerate(bs)):
        # ensure b is an int, not a float

        # loop through obs or ks depending on which we want to plot
        if obs_phase == True:
            for j, ob in enumerate(obs):
                pop = new_pop(infected, nrow, ncol)
                pop, S, I, R = abm_pop_sim(pop, b, ob, .01, t)
                cts[i, j] = np.true_divide(I[-1], N)
        else:
            for j, k in enumerate(ks):
                pop = new_pop(infected, nrow, ncol)
                pop, S, I, R = abm_pop_sim(pop, b, 0, k, t)
                cts[i, j] = np.true_divide(I[-1], N)


    # Create a phase plot with variation depending on if using k or ob as variable
    fig, ax = plt.subplots()
    if obs_phase == True:
        axcontour = ax.contour(obs, bs, cts)
    else:
        axcontour = ax.contour(ks, bs, cts)
    fig.colorbar(axcontour)
    ax.set_title('phase diagram [I, t={}]'.format(t))
    if obs_phase == True:
        ax.set_xlabel('Chance of a Non-local Interaction')
        ax.set_ylabel('Number of Interactions')
    else:
        ax.set_xlabel('k')
        ax.set_ylabel('b')
    

    # Save plot
    if save_path is not None:
        fig.savefig(save_path)
    plt.show()

def new_pop(start_inf, nrow, ncol):
    """
    Create an initial state population of size nrow * ncol which lives on a grid
    :param nrow: number of rows in grid
    :param ncol: number of columns in grid
    :param start_inf: count of infected initial state
    """
    # Create a population of susceptibles
    pop = [Person() for i in range(nrow * ncol)]
    pop_size = len(pop)

    infected = random.sample(range(pop_size), start_inf)
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

def make_cohorts(pop, size=9):
    """
    Splits the population into cohorts of a given size
    """
    nrow, ncol = pop.shape

    # Ensure cohort size is a perfect square, and a factor of population size
    assert pop.size % size == 0 , 'Cohort size is not a factor of population size'
    root_size = np.int(np.sqrt(size))
    assert root_size == np.floor(root_size), 'Cohort size is not a perfect square'

    # Split population grid into chunks of columns
    rsplit = np.split(pop, nrow/root_size)
    cohort_num = 0

    for r in range(len(rsplit)):

        # Splits column chunks into cohorts
        csplit = np.hsplit(rsplit[r], ncol/root_size)
        for c in range(len(csplit)):
            cohort_num = cohort_num + 1

            # Assigns group number to every agent within cohort
            for i in range(root_size):
                for j in range(root_size):
                    csplit[c][i,j].group(cohort_num)

    return pop


def abm_2D_sim(pop, p, q, k, t):
    """
    Simulate the spread of a disease on a given population over t days.
    pop designates the starting state, b and k control the spread and recovery rate
    """
    S = []
    I = []
    R = []

    for i in range(t):
        pop = move_pop(pop, p, q)
        pop = remove_pop(pop, k)
        S.append(len(get_indices(pop, 'S')))
        I.append(len(get_indices(pop, 'I')))
        R.append(len(get_indices(pop, 'R')))

    return pop, S, I, R


def move_pop(pop, p, q):
    """
    Have each agent move p distance and interact with other agents
    in a q radius.
    Ignores recovered agents as they have no meaningful interactions.
    """
    # Use flattened pop
    pop = pop.flatten()
 
    # Make a matrix of current positions
    X = np.zeros((pop.size, 2))
    for i in range(pop.size):
        X[i,0] = pop[i].pos[0]
        X[i,1] = pop[i].pos[1]

    tree = KDTree(X)

    # Find and move susceptible agents
    susceptible = get_indices(pop, 'S')
    for s in susceptible:
        pop[s].move(p)
        inds = tree.query_ball_point(pop[s].pos, q)

        for j in inds:
            if pop[j].state == 'I':
                pop[s].infect()
                break

    # Find and move infected agents
    infected = get_indices(pop, 'I')
    for i in infected:
        pop[i].move(p)
        inds = tree.query_ball_point(pop[i].pos, q)
        for j in inds:
            pop[j].infect()

    return pop