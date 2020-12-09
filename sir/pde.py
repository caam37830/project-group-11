import scipy.sparse as sparse
from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

# create matrix A to apply forward difference scheme
def forward_diff_matrix(n):
    data = []
    i = []
    j = []
    for k in range(n - 1):
        i.append(k)
        j.append(k)
        data.append(-1)

        i.append(k)
        j.append(k+1)
        data.append(1)

    return sparse.coo_matrix((data, (i, j)), shape=(n, n)).tocsr()

# delta matrix
def delta(shape: tuple):
    m, n = shape
    d = forward_diff_matrix(n)
    DD = -d.T @ d
    k = np.eye(m)
    DDx = sparse.kron(k, DD)
    DDy = sparse.kron(DD, k)
    return DDx + DDy

class covid2D():
    """
    2D spatial SIR model

    Methods for projecting and visualizing the trajectory given certain parameters are implemented
    """
    def __init__(self, b=3, k=0.1, p=1, M=200, sourcepoint=1, source_loc=None, **kwargs):
        """
        initialize class with initial values for each compartment
        """
        # universal parameters, child class will have them
        self.S = kwargs.get('S', 999_995)
        self.I = kwargs.get('I', 5)
        self.R = kwargs.get('R', 0)
        self.N = kwargs.get('N', self.S + self.I + self.R)

        assert b > 0, 'b must be a positive number'
        assert (k > 0 and k < 1), 'k must between 0 and 1'

        self.b = b
        self.k = k
        self.p = p

        self.sol = None

        # init flattened s,i,r matrix
        self.i = np.zeros(M * M)

        if isinstance(sourcepoint, int):
            if source_loc is None:
                source_index = np.random.randint(0, M*M, sourcepoint)
            else:
                assert sourcepoint == len(source_loc[0]), '# of source point should = len(source_loc[0])'
                source_index = np.ravel_multi_index(source_loc, dims=(M, M))
        else:
            source_index = np.arange(M*M)
        self.i[source_index] = self.I / self.N
        self.s = 1 - self.i
        self.r = np.ones(M * M) * self.R / self.N
        self.L = delta((M, M))

        # model based parameters, child class should redefine these parameters
        self.y0 = np.array([self.s, self.i, self.r]).flatten()
        self.s_idx = np.arange(M * M)
        self.i_idx = np.arange(M * M, 2 * M * M)
        self.r_idx = np.arange(2 * M * M, 3 * M * M)
        self.idx_dic = {0: self.s_idx, 1: self.i_idx, 2: self.r_idx}
        self.labels = ['Susceptible', 'Infectious', 'Removed']
        self.colors = ['green', 'red', 'blue']

    def __call__(self, t):
        """
        return proportion of each group at particular time
        """
        # raise error if the solve method has not ran
        if self.sol is None:
            raise AssertionError('Covid class is callable only after solve has been run')

        return self.sol.sol(t)

    def __repr__(self):
        """
        model representation
        """
        return "Covid_SIR2D(b={}, k={}, p={})".format(self.b, self.k, self.p)

    def rhs(self, t, y):
        """
        Define sir model's differential equations
            ds(x, t) / dt = -b * s(x, t) * i(x, t) + p * L * s(x, t)
            di(x, t) / dt = b * s(x, t) * i(x, t) - k * i(x, t) + p * L * i(x, t)
            dr(x, t) / dt = k * i(x, t) + p * L * r(x, t)
        """
        s_ = -self.b * y[self.s_idx] * y[self.i_idx] + self.p * self.L @ y[self.s_idx]
        i_ = self.b * y[self.s_idx] * y[self.i_idx] - self.k * y[self.i_idx] + self.p * self.L @ y[self.i_idx]
        r_ = self.k * y[self.i_idx] + self.p * self.L @ y[self.r_idx]
        return np.array([s_, i_, r_]).flatten()

    def solve(self, t_bound=400, h=1, **kwargs):
        """
        solve this ODE system, with RHS=self.rhs, y0=self.y0,
        return
            self.sol - ode.OdeResult object, bunch object of ode results: t,y,sol, etc.

        parameters
            t_bound - compute this system from day(0) to day(t_bound)
            h - update step, default is 1 day.

        """
        t_span = (0, t_bound)
        t_eval = np.arange(0, t_bound, h)

        self.sol = solve_ivp(fun=self.rhs, t_span=t_span, y0=self.y0,
                             t_eval=t_eval, dense_output=True, **kwargs)

        return self.sol

    def tracePlot(self, save_path=None, decorators=True, show=True):
        """
        plot simulation result
        """
        n_com = len(self.idx_dic)
        # check the length of labels
        if len(self.labels) != n_com:
            self.labels = [str(i+1) for i in range(n_com)]
        # check the length of colors
        if len(self.colors) != n_com:
            self.colors = [None] * n_com
        # loop and plot
        for i in range(n_com):
            plt.plot(self.sol.t, np.mean(self.sol.y[self.idx_dic[i]], axis=0), label=self.labels[i], color=self.colors[i])

        if decorators:
            plt.title(self.__repr__())
            plt.ylabel("ratio")
            plt.xlabel("day")
            plt.legend()

        if save_path is not None:
            plt.savefig(save_path)

        if show:
            plt.show()

