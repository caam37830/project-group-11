import unittest
from sir.ode import *

class TestABM(unittest.TestCase):
    def setUp(self):
        pass

    def test_something(self):
        self.assertEqual(True, False)


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


if __name__ == '__main__':
    unittest.main()
