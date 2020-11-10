from scipy.integrate import solve_ivp
import numpy as np
import sympy as sym
import matplotlib.pyplot as plt

class covid():
    """
    A population in which Covid is present

    Methods for projecting and visualizing the trajectory given certain parameters are implemented
    """
    def __init__(self, S, I, R, b, k, **kwargs):
        """
        initialize class with initial values for each compartment

        b is number of interactions per individual per day
        k is fraction of infectious period which recovers each day (0 < k < 1)
        """
        # todo 1. here I introduce **kwargs for further variations
        assert type(S) == int, 'S should be int'
        assert type(I) == int, 'I should be int'
        assert type(R) == int, 'R should be int'
        assert b > 0, 'b must be a positive number'
        assert (k > 0 and k < 1), 'k must be between 0 and 1' 

        self.S = S
        self.I = I
        self.R = R
        self.N = S + I + R

        self.s = S / self.N
        self.i = I / self.N
        self.r = R / self.N

        self.k = k
        self.b = b

        self.sol = None
        # below, we should add more init for further variations.

    def __call__(self, t):
        """
        return proportion of each group at particular time
        """
        # raise error if the solve method has not ran
        if self.sol is None:
            raise AssertionError('Covid class is callable only after solve has been run')

        return self.sol.sol(t)

    def __repr__(self):
        return "Covid({}, {}, {})(b={}, k={})".format(self.s, self.i, self.r, self.b, self.k)


    def solve(self, t_bound, fun=None, h=1, **kwargs):
        """
        define ode solve for this system
        input t_bound to simulation for t days
        default function is self.fun, default step size is h=1
        """
        if fun is None:
            def fun(t, y):
                """
                define initial function, the basic sir differential equations
                """
                s_ = -self.b * y[0] * y[1]
                i_ = self.b * y[0] * y[1] - self.k * y[1]
                r_ = self.k * y[1]
                return np.array([s_, i_, r_])

        y0 = np.array([self.s, self.i, self.r])
        t_span = (0, t_bound)
        t_eval = np.arange(0, t_bound, h)

        self.sol = solve_ivp(fun, t_span, y0=y0, t_eval=t_eval, dense_output=True, **kwargs)

        return self.sol

    def plot(self):
        """
        integrate plot related lines
        """
        plt.plot(self.sol.t, self.sol.y[0], label='Susceptible', color='green')
        plt.plot(self.sol.t, self.sol.y[1], label='Infectious', color='red')
        plt.plot(self.sol.t, self.sol.y[2], label='Removed', color='blue')
        plt.title("SIR ODE simulation, b={},k={}".format(self.b, self.k))
        plt.ylabel("ratio")
        plt.xlabel("day")
        plt.legend()
        plt.show()

    def event_solve(self):
        """
        under some condition, the society take some actions, then solve the problem.
        """
        pass
