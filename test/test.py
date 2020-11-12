import unittest
from sir.ode import *
from sir.abm import *

class TestABM(unittest.TestCase):
    def setUp(self):
        pass

    def test_Remove(self):

        infected_Test = Person('I')
        infected_Test.remove()

        self.assertEqual(infected_Test.cur_state(), 'R')

    def test_Infect(self):

        suscept_Test = Person('S')
        suscept_Test.infect()

        self.assertEqual(suscept_Test.cur_state(), 'I')
    
    def test_closed_abm(self):
        
        N = 100_000
        pop = [Person() for i in range(N)]
        for i in range(100):
            pop[i].infect()

        T = 100

        pop, S, I, R = abm_pop_sim(pop, 4, 0.1, T)

        s = np.true_divide(S, N)
        i = np.true_divide(I, N)
        r = np.true_divide(R, N)
        
        n = s + i + r

        self.assertTrue(np.all(np.abs(n-1)) < 1e-8)
        

class TestODE(unittest.TestCase):
    def setUp(self):
        pass

    def test_ode(self):
        N = 1_000_000
        I = 5
        R = 0
        S = N - I - R
        k = 0.01
        b = 1
        ode_covid = covid(S, I, R, b, k)
        ode_covid.solve(t_bound=100)
        t = np.arange(5, 95, 5)
        s, i, r = ode_covid(t)

        lhs = (ode_covid(t + 0.1) - ode_covid(t)) / 0.1
        rhs = np.array([-b * s * i, b * s * i - k * i, k * i])
        self.assertTrue(np.all(np.abs(lhs - rhs) < 0.01))
    
    def test_closed_ode(self):
        N = 1_000_000
        I = 5
        R = 0
        S = N - I - R
        k = 0.01
        b = 1

        ode_covid = covid(S, I, R, b, k)
        ode_covid.solve(t_bound=100)
        t = np.arange(5, 95, 5)
        s, i, r = ode_covid(t)
        n = s + i + r

        self.assertTrue(np.all(np.abs(n - 1)) < 1e-8)


if __name__ == '__main__':
    unittest.main()
