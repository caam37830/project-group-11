from scipy.integrate import solve_ivp
import numpy as np
import sympy as sym
import matplotlib.pyplot as plt

class population:
    def __init__(self, S, I, R):
        assert type(S) == int, 'S should be int'
        assert type(I) == int, 'I should be int'
        assert type(R) == int, 'R should be int'

        self.S = S
        self.I = I
        self.R = R
        self.N = S + I + R

        self.s = S / self.N
        self.i = I / self.N
        self.r = R / self.N

class covid(population):
    def __init__(self, S, I ,R, k, b, **kwargs):
        # todo 1. here I introduce **kwargs for further variations
        super(covid, self).__init__(S, I, R)

        self.k = k
        self.b = b
        # below, we should add more init for further variations.

    def fun(self, t, y):
        s_ = -self.b * y[0] * y[1]
        i_ = self.b * y[0] * y[1] - self.k * y[1]
        r_ = self.k * y[1]
        return np.array([s_, i_, r_])

    def solve(self, t_span, fun=None, **kwargs):
        if fun is None:
            fun = self.fun
        y0 = np.array([self.s, self.i, self.r])
        self.sol = solve_ivp(fun, t_span, y0=y0, **kwargs)
        return self.sol

    def plot(self, *args, **kwargs):
        plt.plot(self.sol.t, self.sol.y[0], label='Susceptible', col='green')
        plt.plot(self.sol.t, self.sol.y[1], label='Infectious', col='red')
        plt.plot(self.sol.t, self.sol.y[2], label='Removed', col='blue')
        plt.title("SIR simulation")
        plt.ylabel("ratio")
        plt.xlabel("day")
        plt.legend()
        plt.show()

    def event_solve(self):
        """
        under some condition, the society take some actions, then solve the problem.
        """
        pass


