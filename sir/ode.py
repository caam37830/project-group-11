from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

# todo 1. change plot color
# todo 2. extension: contact rate "b" -> "b(t)"
# todo 3. extension: vaccination related model modification


class covid():
    """
    A population in which Covid is present

    Methods for projecting and visualizing the trajectory given certain parameters are implemented
    """
    def __init__(self, b ,k, **kwargs):
        """
        initialize class with initial values for each compartment
        """
        # universal parameters, child class will have them
        self.S = kwargs.get('S', 999_995)
        self.I = kwargs.get('I', 5)

        assert b > 0, 'b must be a positive number'
        assert (k > 0 and k < 1), 'k must between 0 and 1'
        self.b = b
        self.k = k

        self.sol = None

        # model based parameters, child class should redefine these parameters
        self.y0 = np.array([999_995/1_000_000, 5/1_000_000])
        self.labels = ['Susceptible', 'Infectious']
        self.colors = ['green', 'red']

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
        return "Covid_constant(S={}, I{})".format(self.S, self.I)

    def rhs(self, t, y):
        """
        Right Hand Side (RHS) of Ordinary Differential Equations (ODE)
        null model rhs is zeros
        """
        return np.zeros(len(y))

    def solve(self, t_bound, h=1, **kwargs):
        """
        solve this ODE system, with RHS=self.rhs, y0=self.y0,
        return
            ode.OdeResult - bunch object of ode results: t,y,sol, etc.

        parameters
            t_bound - compute this system from day(0) to day(t_bound)
            h - update step, default is 1 day.

        """
        t_span = (0, t_bound)
        t_eval = np.arange(0, t_bound, h)

        self.sol = solve_ivp(fun=self.rhs, t_span=t_span, y0=self.y0,
                             t_eval=t_eval, dense_output=True, **kwargs)

        return self.sol

    def plot(self, save_path=None):
        """
        plot simulation result
        """
        n_com = len(self.y0)
        # check the length of labels
        if len(self.labels) != n_com:
            self.labels = [str(i+1) for i in range(n_com)]
        # check the length of colors
        if len(self.colors) != n_com:
            self.colors = [None] * n_com
        # loop and plot
        for i in range(n_com):
            plt.plot(self.sol.t, self.sol.y[i], label=self.labels[i], color=self.colors[i])
        plt.title(self.__repr__())
        plt.ylabel("ratio")
        plt.xlabel("day")
        plt.legend()
        if save_path is not None:
            plt.savefig(save_path)
        plt.show()

    def event_solve(self):
        """
        under some condition, the society take some actions, then solve the problem.
        """
        pass

def phase_plot(N, I, R, t, phase='I', bs=np.linspace(1, 10, 50), ks=np.linspace(0.01, .5, 50), save_path=None):
    """
    plot phase diagram
    :param N: total number of population
    :param I: infected number of poplation
    :param R: removed number of population
    :param t: time
    :param phase: plot which parameter's phase
    :param bs: discrete b
    :param ks: discrete k
    :param save_path:
    :return:
    """
    idx = 0
    if phase == 'I':
        idx = 1
    elif phase == 'R':
        idx = 2
    else:
        idx = 0

    cts = np.zeros((len(bs), len(ks)))
    for i, b in enumerate(bs):
        for j, k in enumerate(ks):
            ode_covid = covid(S=N-I-R, I=I, R=R, b=b, k=k)
            ode_covid.solve(t_bound=t)
            cts[i, j] = ode_covid.sol.y[idx, -1]

    fig, ax = plt.subplots()
    axcontour = ax.contour(ks, bs, cts)
    fig.colorbar(axcontour)
    ax.set_title('phase diagram [{}, t={}]'.format(phase, t))
    ax.set_xlabel('k')
    ax.set_ylabel('b')

    if save_path is not None:
        fig.savefig(save_path)
    plt.show()


class SIR(covid):
    def __init__(self, b, k, **kwargs):
        """
        init sir model parameters,
            ds / dt = -b * s * i
            di / dt = b * s * i - k * i
            dr / dt = k * i

        Parameter
            b - b is number of interactions per individual per day
            k - k is fraction of infectious period which recovers each day (0 < k < 1)

        Optional Parameter
            S - susceptible population
            I - infectious population
            R - recovered population
            N - total population
        """
        super().__init__(b, k, **kwargs)

        # init model related parameters
        self.R = kwargs.get('R', 0)
        self.N = kwargs.get('N', self.S + self.I + self.R)
        assert self.S + self.I + self.R == self.N, 'S+I+R should equal to N'
        self.s = self.S / self.N
        self.i = self.I / self.N
        self.r = self.R / self.N

        # redefine self.y0, self.labels, self.colors.
        self.y0 = np.array([self.s, self.i, self.r])
        self.labels = ['Susceptible', 'Infectious', 'Removed']
        self.colors = ['green', 'red', 'blue']

    def __repr__(self):
        return "Covid_SIR(s={}, i={}, r={})(b={}, k={})".format(self.s, self.i, self.r, self.b, self.k)

    def rhs(self, t, y):
        """
        Define sir model's differential equations
        """
        s_ = -self.b * y[0] * y[1]
        i_ = self.b * y[0] * y[1] - self.k * y[1]
        r_ = self.k * y[1]
        return np.array([s_, i_, r_])


class SIS(covid):
    def __init__(self, b, k, **kwargs):
        """
        init sis model parameters,
            ds / dt = - b * s * i + k * i
            di / dt = b * s * i - k * i

        Parameter
            b - b is number of interactions per individual per day
            k - k is fraction of infectious period which recovers each day (0 < k < 1)

        Optional Parameter
            S - susceptible population
            I - infectious population
            N - total population
        """
        super().__init__(b, k, **kwargs)

        # init sis model related parameters
        self.N = kwargs.get('N', self.S + self.I)
        assert self.S + self.I == self.N, 'S+I should equal to N'
        self.s = self.S / self.N
        self.i = self.I / self.N

        # redefine self.y0, self.labels, self.colors.
        self.y0 = np.array([self.s, self.i])
        self.labels = ['Susceptible', 'Infectious']
        self.colors = ['green', 'red']

    def __repr__(self):
        """
        redefine model representation
        """
        return "Covid_SIS(s={}, i={})(b={}, k={})".format(self.s, self.i, self.b, self.k)

    def rhs(self, t, y):
        """
        Define sir model's differential equations
        """
        s_ = -self.b * y[0] * y[1] + self.k * y[1]
        i_ = self.b * y[0] * y[1] - self.k * y[1]
        return np.array([s_, i_])


class SIRD(covid):
    def __init__(self, b, k, mu, **kwargs):
        """
        init SIRD model parameters,
            ds / dt = -b * s * i
            di / dt = b * s * i - k * i - mu * i
            dr / dt = k * i
            dd / dt = mu * i

        Parameter
            b - b is number of interactions per individual per day
            k - k is fraction of infectious period which recovers each day (0 < k < 1)
            mu - mu is the rate of mortality

        Optional Parameter
            S - susceptible population
            I - infectious population
            R - recovered population
            D - decreased population
            N - total population
        """
        super().__init__(b, k, **kwargs)

        # init model related parameters
        self.mu = mu
        self.R = kwargs.get('R', 0)
        self.D = kwargs.get('D', 0)
        self.N = kwargs.get('N', self.S + self.I + self.R + self.D)
        assert self.S + self.I + self.R + self.D == self.N, 'S+I+R+D should equal to N'
        self.s = self.S / self.N
        self.i = self.I / self.N
        self.r = self.I / self.N
        self.d = self.D / self.N

        # redefine self.y0, self.labels, self.colors.
        self.y0 = np.array([self.s, self.i, self.r, self.d])
        self.labels = ['Susceptible', 'Infectious', 'Removed', 'Decreased']
        self.colors = ['green', 'red', 'blue', 'black']

    def __repr__(self):
        """
        redefine the model representation
        """
        return f"Covid_SIRD(s={self.s}, i={self.i}, r={self.r}, d={self.d})(b={self.b}, k={self.k}, mu={self.mu})"

    def rhs(self, t, y):
        """
        Define SIRD model's differential equations
        """
        s_ = -self.b * y[0] * y[1]
        i_ = self.b * y[0] * y[1] - self.k * y[1] - self.mu * y[1]
        r_ = self.k * y[1]
        d_ = self.mu * y[1]
        return np.array([s_, i_, r_, d_])


class MSIR(covid):
    def __init__(self, lam, sigma, b, k, mu, **kwargs):
        """
        init msir model parameters,
            dm / dt = lam - sigma * m - mu * m
            ds / dt = sigma * m - b * s * i - mu * s
            di / dt = b * s * i - k * i - mu * i
            dr / dt = k * i - mu * r

        Parameter
            lam - the rate of new born
            sigma - the rate of maternally derived immunity change to susceptible population
            b - b is number of interactions per individual per day
            k - k is fraction of infectious period which recovers each day (0 < k < 1)
            mu - mu is the rate of mortality

        Optional Parameter
            M - maternally derived immunity population
            S - susceptible population
            I - infectious population
            R - recovered population
            N - total population
        """
        super().__init__(b, k, **kwargs)

        # init model related parameters
        self.lam = lam
        self.sigma = sigma
        self.mu = mu
        self.R = kwargs.get('R', 0)
        self.M = kwargs.get('M', 0)
        self.N = kwargs.get('N', self.S + self.I + self.R + self.M)
        assert self.S + self.I + self.R + self.M == self.N, 'M+S+I+R should equal to N'
        self.s = self.S / self.N
        self.i = self.I / self.N
        self.r = self.I / self.N
        self.m = self.M / self.N

        # redefine self.y0, self.labels, self.colors.
        self.y0 = np.array([self.m, self.s, self.i, self.r])
        self.labels = ['MDI', 'Susceptible', 'Infectious', 'Removed']
        self.colors = ['grey', 'green', 'red', 'blue']

    def __repr__(self):
        """
        redefine the model representation
        """
        return f"Covid_MSIR(m={self.m}, s={self.s}, i={self.i}, r={self.r})(lam={self.lam}, sigma={self.sigma}, b={self.b}, k={self.k}, mu={self.mu})"

    def rhs(self, t, y):
        """
        Define MSIR model's differential equations
        """
        m_ = self.lam - self.sigma * y[0] - self.mu * y[0]
        s_ = self.sigma * y[0] - self.b * y[1] * y[2] - self.mu * y[1]
        i_ = self.b * y[1] * y[2] - self.k * y[2] - self.mu * y[2]
        r_ = self.k * y[2] - self.mu * y[3]
        return np.array([m_, s_, i_, r_])


class SIRC(covid):
    def __init__(self, b, k, c1, c2, **kwargs):
        """
        init SIRC model parameters, [susceptible-infectious-recovered-carrier]
            ds / dt = - b * s * i
            di / dt = b * s * i - k * i - c1 * i + c2 * c
            dc / dt = c1 * i - c2 * c
            dr / dt = k * i

        Parameter
            b - b is number of interactions per individual per day
            k - k is fraction of infectious period which recovers each day (0 < k < 1)
            c1 - the rate of infectious to carrier
            c2 - the rate of carrier to infectious

        Optional Parameter
            S - susceptible population
            I - infectious population
            R - recovered population
            C - carrier population
            N - total population
        """
        super().__init__(b, k, **kwargs)

        # init model related parameters
        self.c1 = c1
        self.c2 = c2
        self.R = kwargs.get('R', 0)
        self.C = kwargs.get('C', 0)
        self.N = kwargs.get('N', self.S + self.I + self.R + self.C)
        assert self.S + self.I + self.R + self.C == self.N, 'S+I+R+C should equal to N'
        self.s = self.S / self.N
        self.i = self.I / self.N
        self.r = self.I / self.N
        self.c = self.C / self.N

        # redefine self.y0, self.labels, self.colors.
        self.y0 = np.array([self.s, self.i, self.c,self.r])
        self.labels = ['Susceptible', 'Infectious', 'Carrier', 'Removed']
        self.colors = ['green', 'red', 'grey', 'blue']

    def __repr__(self):
        """
        redefine the model representation
        """
        return f"Covid_SIRC(s={self.s}, i={self.i}, r={self.r}, c={self.c})(b={self.b}, k={self.k}, c1={self.c1}, c2={self.c2})"

    def rhs(self, t, y):
        """
        Define SIRC model's differential equations
        """
        s_ = - self.b * y[0] * y[1]
        i_ = self.b * y[0] * y[1] - self.k * y[1] - self.c1 * y[1] + self.c2 * y[2]
        c_ = self.c1 * y[1] - self.c2 * y[2]
        r_ = self.k * y[1]
        return np.array([s_, i_, c_, r_])


class SEIR(covid):
    def __init__(self, lam, b, k, a, mu, **kwargs):
        """
        init SEIR model parameters, [susceptible-exposed-infectious-recovered]
            ds / dt = lam - mu * s - b * s * i
            de / dt = b * s * i - (mu + a) * e
            di / dt = a * e - (k + mu) * i
            dr / dt = k * i - mu * r

        Parameter
            lam - birth rate of total population
            b - b is number of interactions per individual per day
            k - k is fraction of infectious period which recovers each day (0 < k < 1)
            a - 1/a is the mean incubation period of exponential distribution
            mu - population decrease rate

        Optional Parameter
            S - susceptible population
            E - exposed population
            I - infectious population
            R - recovered population
            N - total population
        """
        super().__init__(b, k, **kwargs)

        # init model related parameters
        self.lam = lam
        self.a = a
        self.mu = mu
        self.R = kwargs.get('R', 0)
        self.E = kwargs.get('E', 0)
        self.N = kwargs.get('N', self.S + self.I + self.R + self.E)
        assert self.S + self.I + self.R + self.E == self.N, 'S+E+I+R should equal to N'
        self.s = self.S / self.N
        self.i = self.I / self.N
        self.r = self.I / self.N
        self.e = self.E / self.N

        # redefine self.y0, self.labels, self.colors.
        self.y0 = np.array([self.s, self.e, self.i, self.r])
        self.labels = ['Susceptible', 'Exposed', 'Infectious', 'Removed']
        self.colors = ['green', 'yellow', 'red', 'blue']

    def __repr__(self):
        """
        redefine the model representation
        """
        return f"Covid_SEIR(s={self.s}, e={self.e}, i={self.i}, r={self.r})(lam={self.lam}, b={self.b}, k={self.k}, a={self.a}, mu={self.mu})"

    def rhs(self, t, y):
        """
        Define SEIR model's differential equations
        """
        s_ = self.lam - self.b * y[0] * y[2] - self.mu * y[0]
        e_ = self.b * y[0] * y[2] - (self.mu + self.a) * y[1]
        i_ = self.a * y[1] - (self.mu + self.k) * y[2]
        r_ = self.k * y[2] - self.mu * y[3]
        return np.array([s_, e_, i_, r_])


class SEIS(covid):
    def __init__(self, lam, b, k, a, mu, **kwargs):
        """
        init SEIS model parameters, [susceptible-exposed-infectious-susceptible]
            ds / dt = lam - mu * s - b * s * i + k * i
            de / dt = b * s * i - (mu + a) * e
            di / dt = a * e - (k + mu) * i

        Parameter
            lam - birth rate of total population
            b - b is number of interactions per individual per day
            k - k is fraction of infectious period which recovers each day (0 < k < 1)
            a - 1/a is the mean incubation period of exponential distribution
            mu - population decrease rate

        Optional Parameter
            S - susceptible population
            E - exposed population
            I - infectious population
            N - total population
        """
        super().__init__(b, k, **kwargs)

        # init model related parameters
        self.lam = lam
        self.a = a
        self.mu = mu
        self.E = kwargs.get('E', 0)
        self.N = kwargs.get('N', self.S + self.I + self.E)
        assert self.S + self.I + self.E == self.N, 'S+E+I should equal to N'
        self.s = self.S / self.N
        self.i = self.I / self.N
        self.e = self.E / self.N

        # redefine self.y0, self.labels, self.colors.
        self.y0 = np.array([self.s, self.e, self.i])
        self.labels = ['Susceptible', 'Exposed', 'Infectious']
        self.colors = ['green', 'yellow', 'red']

    def __repr__(self):
        """
        redefine the model representation
        """
        return f"Covid_SEIS(s={self.s}, e={self.e}, i={self.i})(lam={self.lam}, b={self.b}, k={self.k}, a={self.a}, mu={self.mu})"

    def rhs(self, t, y):
        """
        Define SEIR model's differential equations
        """
        s_ = self.lam - self.b * y[0] * y[2] - self.mu * y[0] + self.k * y[2]
        e_ = self.b * y[0] * y[2] - (self.mu + self.a) * y[1]
        i_ = self.a * y[1] - (self.mu + self.k) * y[2]
        return np.array([s_, e_, i_])


class MSEIR(covid):
    def __init__(self, lam, sigma, b, k, a, mu, **kwargs):
        """
        init MSEIR model parameters, [maternally_derived_immunity-susceptible-exposed-infectious-recovered]
            dm / dt = lam - sigma * m - mu * m
            ds / dt = sigma * m - mu * s - b * s * i
            de / dt = b * s * i - (mu + a) * e
            di / dt = a * e - (k + mu) * i
            dr / dt = k * i - mu * r

        Parameter
            lam - birth rate of total population
            sigma - the rate of changing from maternally_derived_immunity to susceptible
            b - b is number of interactions per individual per day
            k - k is the fraction of infectious to recovered each day (0 < k < 1)
            a - 1/a is the mean incubation period of exponential distribution
            mu - population decrease rate

        Optional Parameter
            M - maternally-derived-immunity population
            S - susceptible population
            E - exposed population
            I - infectious population
            R - recovered population
            N - total population
        """
        super().__init__(b, k, **kwargs)

        # init model related parameters
        self.lam = lam
        self.sigma = sigma
        self.a = a
        self.mu = mu
        self.M = kwargs.get('M', 0)
        self.R = kwargs.get('R', 0)
        self.E = kwargs.get('E', 0)
        self.N = kwargs.get('N', self.S + self.I + self.R + self.E + self.M)
        assert self.S + self.I + self.R + self.E + self.M == self.N, 'M+S+E+I+R should equal to N'
        self.m = self.M / self.N
        self.s = self.S / self.N
        self.i = self.I / self.N
        self.r = self.I / self.N
        self.e = self.E / self.N

        # redefine self.y0, self.labels, self.colors.
        self.y0 = np.array([self.m, self.s, self.e, self.i, self.r])
        self.labels = ['MDI', 'Susceptible', 'Exposed', 'Infectious', 'Removed']
        self.colors = ['grey', 'green', 'yellow', 'red', 'blue']

    def __repr__(self):
        """
        redefine the model representation
        """
        return f"Covid_MSEIR(m={self.m}, s={self.s}, e={self.e}, i={self.i}, r={self.r})(lam={self.lam}, sigma={self.sigma}, b={self.b}, k={self.k}, a={self.a}, mu={self.mu})"

    def rhs(self, t, y):
        """
        Define MSEIR model's differential equations
        """
        m_ = self.lam - self.sigma * y[0] - self.mu * y[0]
        s_ = self.sigma * y[0] - self.b * y[1] * y[3] - self.mu * y[1]
        e_ = self.b * y[1] * y[3] - (self.mu + self.a) * y[2]
        i_ = self.a * y[2] - (self.mu + self.k) * y[3]
        r_ = self.k * y[3] - self.mu * y[4]
        return np.array([m_, s_, e_, i_, r_])


class MSEIRS(covid):
    def __init__(self, lam, sigma, b, k, a, mu, l, **kwargs):
        """
        init MSEIRS model parameters, [maternally_derived_immunity-susceptible-exposed-infectious-recovered-susceptible]
            dm / dt = lam - sigma * m - mu * m
            ds / dt = sigma * m + l * r - mu * s - b * s * i
            de / dt = b * s * i - (a + mu) * e
            di / dt = a * e - (k + mu) * i
            dr / dt = k * i - (l + mu) * r

        Parameter
            lam - birth rate of total population
            sigma - the rate of changing from maternally_derived_immunity to susceptible
            b - b is number of interactions per individual per day
            k - k is the fraction of infectious to recovered each day (0 < k < 1)
            a - 1/a is the mean incubation period of exponential distribution
            mu - population decrease rate
            l - temporary immunity R would become S, 1/l is the mean immunity period of exponential distribution

        Optional Parameter
            M - maternally-derived-immunity population
            S - susceptible population
            E - exposed population
            I - infectious population
            R - recovered population
            N - total population
        """
        super().__init__(b, k, **kwargs)

        # init model related parameters
        self.lam = lam
        self.sigma = sigma
        self.a = a
        self.mu = mu
        self.l = l
        self.M = kwargs.get('M', 0)
        self.R = kwargs.get('R', 0)
        self.E = kwargs.get('E', 0)
        self.N = kwargs.get('N', self.S + self.I + self.R + self.E + self.M)
        assert self.S + self.I + self.R + self.E + self.M == self.N, 'M+S+E+I+R should equal to N'
        self.m = self.M / self.N
        self.s = self.S / self.N
        self.i = self.I / self.N
        self.r = self.I / self.N
        self.e = self.E / self.N

        # redefine self.y0, self.labels, self.colors.
        self.y0 = np.array([self.m, self.s, self.e, self.i, self.r])
        self.labels = ['MDI', 'Susceptible', 'Exposed', 'Infectious', 'Removed']
        self.colors = ['grey', 'green', 'yellow', 'red', 'blue']

    def __repr__(self):
        """
        redefine the model representation
        """
        return f"Covid_MSEIRS(m={self.m}, s={self.s}, e={self.e}, i={self.i}, r={self.r})(lam={self.lam}, sigma={self.sigma}, b={self.b}, k={self.k}, a={self.a}, mu={self.mu}, l={self.l})"

    def rhs(self, t, y):
        """
        Define MSEIR model's differential equations
        """
        m_ = self.lam - self.sigma * y[0] - self.mu * y[0]
        s_ = self.sigma * y[0] + self.l * y[4] - self.b * y[1] * y[3] - self.mu * y[1]
        e_ = self.b * y[1] * y[3] - (self.mu + self.a) * y[2]
        i_ = self.a * y[2] - (self.mu + self.k) * y[3]
        r_ = self.k * y[3] - (self.mu + self.l) * y[4]
        return np.array([m_, s_, e_, i_, r_])
