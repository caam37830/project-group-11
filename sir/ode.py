from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import OptimizeResult

# todo 1. change plot color
# todo 2. extension: contact rate "b" -> "b(t)"
# todo 3. extension: vaccination related model modification
# todo 4. phase plot need to be updated


class covid():
    """
    A population in which Covid is present

    Methods for projecting and visualizing the trajectory given certain parameters are implemented
    """
    def __init__(self, b, k, **kwargs):
        """
        initialize class with initial values for each compartment
        """
        # universal parameters, child class will have them
        self.S = kwargs.get('S', 999_995)
        self.I = kwargs.get('I', 5)

        assert b > 0, 'b must be a positive number'
        assert (k > 0 and k < 1), 'k must between 0 and 1'

        self.parameters = {}
        self.parameters['b'] = b
        self.parameters['k'] = k
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

    def solve(self, t_bound=400, h=1, y0=None, **kwargs):
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
        if y0 is None:
            y0 = self.y0
        self.sol = solve_ivp(fun=self.rhs, t_span=t_span, y0=y0,
                             t_eval=t_eval, dense_output=True, **kwargs)

        return self.sol

    def plot(self, save_path=None, decorators=True, show=True):
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

        if decorators:
            plt.title(self.__repr__())
            plt.ylabel("ratio")
            plt.xlabel("day")
            plt.legend()

        if save_path is not None:
            plt.savefig(save_path)

        if show:
            plt.show()


class OdeResult(OptimizeResult):
    pass


def intervention_solve(model, events, t_bound=400, h=1, **kwargs):
    """
    under some condition, the society take some actions, then solve the problem.
    """
    # check terminal
    assert isinstance(model, MSEIQRDS), 'should use most general model'
    for event in events:
        assert event.terminal == True, 'all intervention would cause model change, so event.terminal must be True'
        try:
            assert isinstance(event.intervention, dict), 'events intervention should be a dict'
        except:
            raise AttributeError('this event function object has not define a attribute named "intervention"')

    # if some events happen, solve it until time is up.
    def findidx(l):
        for i, item in enumerate(l):
            if len(item) > 0:
                return i
        return None

    sol = model.solve(t_bound, h=h, events=events)
    event_idx = findidx(sol.t_events)

    if event_idx is None:
        return OdeResult(t=sol.t, y=sol.y, t_events=None, interventions=None)

    # implement model.solve until t_solved >= t_bound
    t_events = []
    interventions = []

    t_solved = sol.t[-1]
    t_events.append(t_solved)
    new_intervention = events[event_idx].intervention
    interventions.append(new_intervention.copy())
    intervention_dic = new_intervention
    y = sol.y
    y0 = sol.y[:, -1]
    events.remove(events[event_idx])
    while t_solved < t_bound-1:
        model = MSEIQRDS(**intervention_dic)
        sol = model.solve(t_bound=t_bound-t_solved, y0=y0, h=h, events=events)

        # update state variables
        t_solved += sol.t[-1]
        y = np.hstack((y, sol.y[:, 1:]))
        y0 = sol.y[:, -1]

        # if solve to t_bound, then directly return results.
        if t_solved == t_bound-1:
            return OdeResult(t=np.arange(0, t_solved+1, h), y=y, t_events=t_events, interventions=interventions)

        # else find intervention and update intervention dict
        t_events.append(t_solved)
        event_idx = findidx(sol.t_events)
        new_intervention = events[event_idx].intervention
        interventions.append(new_intervention.copy())
        intervention_dic.update(new_intervention)
        events = events.remove(events[event_idx])

    return OdeResult(t=np.arange(0, t_solved+1, h), y=y, t_events=t_events, interventions=interventions)


def intervention_plot(sol, save_path=None, show=True):
    model = MSEIQRDS()
    for i in range(6):
        plt.plot(sol.t, sol.y[i], label=model.labels[i], color=model.colors[i])
    for i, t in enumerate(sol.t_events):
        plt.vlines(t, ymin=0, ymax=1, colors='black', linestyles='dashed', label='%dst intervention' % (i+1))
    plt.title('MSEIQRDS model with interventions\n'+str(sol.interventions))
    plt.ylabel("ratio")
    plt.xlabel("day")
    plt.legend()

    if save_path is not None:
        plt.savefig(save_path)
    if show:
        plt.show()



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
            ode_covid = SIR(b=b, k=k, S=N-I-R, I=I, R=R)
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
        return f"Covid_SIRD(s={self.s}, i={self.i}, r={self.r}, d={self.d})\n(b={self.b}, k={self.k}, mu={self.mu})"

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
        self.r = self.R / self.N
        self.m = self.M / self.N

        # redefine self.y0, self.labels, self.colors.
        self.y0 = np.array([self.m, self.s, self.i, self.r])
        self.labels = ['MDI', 'Susceptible', 'Infectious', 'Removed']
        self.colors = ['grey', 'green', 'red', 'blue']

    def __repr__(self):
        """
        redefine the model representation
        """
        return f"Covid_MSIR(m={self.m}, s={self.s}, i={self.i}, r={self.r})\n(lam={self.lam}, sigma={round(self.sigma,4)}, b={self.b}, k={self.k}, mu={self.mu})"

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
        self.r = self.R / self.N
        self.c = self.C / self.N

        # redefine self.y0, self.labels, self.colors.
        self.y0 = np.array([self.s, self.i, self.c,self.r])
        self.labels = ['Susceptible', 'Infectious', 'Carrier', 'Removed']
        self.colors = ['green', 'red', 'grey', 'blue']

    def __repr__(self):
        """
        redefine the model representation
        """
        return f"Covid_SIRC(s={self.s}, i={self.i}, r={self.r}, c={self.c})\n(b={self.b}, k={self.k}, c1={self.c1}, c2={self.c2})"

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
        self.r = self.R / self.N
        self.e = self.E / self.N

        # redefine self.y0, self.labels, self.colors.
        self.y0 = np.array([self.s, self.e, self.i, self.r])
        self.labels = ['Susceptible', 'Exposed', 'Infectious', 'Removed']
        self.colors = ['green', 'yellow', 'red', 'blue']

    def __repr__(self):
        """
        redefine the model representation
        """
        return f"Covid_SEIR(s={self.s}, e={self.e}, i={self.i}, r={self.r})\n(lam={self.lam}, b={self.b}, k={round(self.k, 4)}, a={round(self.a,4)}, mu={self.mu})"

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
        return f"Covid_SEIS(s={self.s}, e={self.e}, i={self.i})\n(lam={self.lam}, b={self.b}, k={self.k}, a={self.a}, mu={self.mu})"

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
        self.r = self.R / self.N
        self.e = self.E / self.N

        # redefine self.y0, self.labels, self.colors.
        self.y0 = np.array([self.m, self.s, self.e, self.i, self.r])
        self.labels = ['MDI', 'Susceptible', 'Exposed', 'Infectious', 'Removed']
        self.colors = ['grey', 'green', 'yellow', 'red', 'blue']

    def __repr__(self):
        """
        redefine the model representation
        """
        return f"Covid_MSEIR(m={self.m}, s={self.s}, e={self.e}, i={self.i}, r={self.r})\n(lam={self.lam}, sigma={round(self.sigma,4)}, b={self.b}, k={self.k}, a={round(self.a,4)}, mu={self.mu})"

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
        self.r = self.R / self.N
        self.e = self.E / self.N

        # redefine self.y0, self.labels, self.colors.
        self.y0 = np.array([self.m, self.s, self.e, self.i, self.r])
        self.labels = ['MDI', 'Susceptible', 'Exposed', 'Infectious', 'Removed']
        self.colors = ['grey', 'green', 'yellow', 'red', 'blue']

    def __repr__(self):
        """
        redefine the model representation
        """
        return f"Covid_MSEIRS(m={self.m}, s={self.s}, e={self.e}, i={self.i}, r={self.r})\n" \
               f"(lam={self.lam}, sigma={round(self.sigma,3)}, b={self.b}, k={self.k}, a={round(self.a,3)}, mu={self.mu}, l={round(self.l,3)})"

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


class MSEIQRDS(covid):
    def __init__(self, **kwargs):
        """
        init MSEIQRSD model parameters,
        [maternally_derived_immunity-susceptible-exposed-infectious-quarantine-recovered-decreased]
        q percent infectious will be become quarantine, they will not spread the virus
        so quarantine will not be included in the ode model,
            dm / dt = lam - sigma * m - mu * m
            ds / dt = sigma * m + re * r - mu * s - b * s * (1-q)i
            de / dt = b * s * (1-q)i - (mu + a) * e
            di / dt = a * e - (k + mu) * i - d * log(i/1-i)
            dd / dt = d * log(i/1-i)
            dr / dt = k * i - (mu +re) * r

        Parameter
            lam - birth rate of total population
            sigma - the rate of changing from maternally_derived_immunity to susceptible
            b - b is number of interactions per individual per day
            k - k is fraction of infectious period which recovers each day (0 < k < 1)
            q - quarantine rate
            a - 1/a is the mean incubation period of exponential distribution
            mu - population decrease rate
            dr - death/decrease rate
            re - re-susceotible rate

        Optional Parameter
            M -
            S - susceptible population
            E - exposed population
            I - infectious population
            R - recovered population
            D -
            N - total population
        """

        # init model related parameters
        # lam=3e-5, sigma=1/720, b=3, k=1/10, q=0, a=1/14, mu=3e-5, d=0.3, re=1/360,
        self.parameters = {}
        self.parameters['lam'] = kwargs.get('lam', 3e-5)
        self.parameters['sigma'] = kwargs.get('sigma', 1/720)
        self.parameters['b'] = kwargs.get('b', 3)
        self.parameters['k'] = kwargs.get('k', 0.1)
        self.parameters['a'] = kwargs.get('a', 1/14)
        self.parameters['mu'] = kwargs.get('mu', 3e-5)
        self.parameters['q'] = kwargs.get('q', 0)
        self.parameters['dr'] = kwargs.get('dr', 0.3)
        self.parameters['re'] = kwargs.get('re', 1/360)

        self.S = kwargs.get('S', 999_995)
        self.I = kwargs.get('I', 5)
        self.R = kwargs.get('R', 0)
        self.E = kwargs.get('E', 0)
        self.M = kwargs.get('M', 0)
        self.D = kwargs.get('D', 0)
        self.N = kwargs.get('N', self.S + self.I + self.R + self.E + self.M + self.D)
        assert self.S + self.I + self.R + self.E + self.M + self.D == self.N, 'M+S+E+I+R+D should equal to N'
        self.s = kwargs.get('s', self.S / self.N)
        self.i = kwargs.get('i', self.I / self.N)
        self.r = kwargs.get('r', self.R / self.N)
        self.e = kwargs.get('e', self.E / self.N)
        self.m = kwargs.get('m', self.M / self.N)
        self.d = kwargs.get('d', self.D / self.N)

        self.sol = None
        # redefine self.y0, self.labels, self.colors.
        self.y0 = np.array([self.m, self.s, self.e, self.i, self.r, self.d])
        self.labels = ['MDI', 'Susceptible', 'Exposed', 'Infectious', 'Death', 'Removed']
        self.colors = ['yellow', 'green', 'grey', 'red', 'black', 'blue']

    def __repr__(self):
        """
        redefine the model representation
        """
        return f"Covid_MSEIQRDS(m={self.m}, s={self.s}, e={self.e}, i={self.i}, r={self.r}, d={self.d})" \
               f"\n(lam={self.parameters['lam']}, sigma={round(self.parameters['sigma'],4)}, b={self.parameters['b']}, k={round(self.parameters['k'], 4)}, q={self.parameters['q']}, " \
               f"\na={round(self.parameters['a'],4)}, mu={self.parameters['mu']}, dr={self.parameters['dr']}, re={round(self.parameters['re'],4)})"

    def rhs(self, t, y):
        """
        Define SEIR model's differential equations
            dm / dt = lam - sigma * m - mu * m
            ds / dt = sigama * m + re * r - mu * s - b * s * (1-q)i
            de / dt = b * s * (1-q)i - (mu + a) * e
            di / dt = a * e - (k + mu) * i - d * i^2 / (1-i)
            dd / dt = d * i^2 / (1-i)
            dr / dt = k * i - (mu + re) * r
        """
        m_ = self.parameters['lam'] - (self.parameters['sigma'] + self.parameters['mu']) * y[0]
        s_ = self.parameters['sigma'] * y[0] + self.parameters['re'] * y[5] - self.parameters['b'] * y[1] * (1 - self.parameters['q']) * y[3] - self.parameters['mu'] * y[1]
        e_ = self.parameters['b'] * y[1] * (1 - self.parameters['q']) * y[3] - (self.parameters['mu'] + self.parameters['a']) * y[2]
        i_ = self.parameters['a'] * y[2] - (self.parameters['mu'] + self.parameters['k']) * y[3] - self.parameters['dr'] * y[3] / (1-y[3]) * y[3]
        d_ = self.parameters['dr'] * y[3] / (1-y[3]) * y[3]
        r_ = self.parameters['k'] * y[3] - (self.parameters['mu'] + self.parameters['re']) * y[5]
        return np.array([m_, s_, e_, i_, d_, r_])
