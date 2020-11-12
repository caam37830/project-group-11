import sys
sys.path.append("../")
from sir.ode import covid, phase_plot

# init population
N = 1_000_000
I = 5
R = 0
S = N - I - R

# init the covid class
ode_covid = covid(S, I, R, b=3, k=0.01)
# solve that class
ode_covid.solve(t_bound=365)
# plot the numerical solution
ode_covid.plot('../doc/checkpoint/figures/ode1.png')
# s, i, r of day 150
print(ode_covid(100))


# init the covid class
ode_covid = covid(S, I, R, b=3, k=0.1)
ode_covid.solve(t_bound=365)
ode_covid.plot('../doc/checkpoint/figures/ode2.png')

# init the covid class
ode_covid = covid(S, I, R, b=0.8, k=0.01)
ode_covid.solve(t_bound=365)
ode_covid.plot('../doc/checkpoint/figures/ode3.png')


phase_plot(N, I, R, t=5, phase='I', save_path='../doc/checkpoint/figures/phase_diagram1.png')
phase_plot(N, I, R, t=10, phase='I', save_path='../doc/checkpoint/figures/phase_diagram2.png')
phase_plot(N, I, R, t=50, phase='I', save_path='../doc/checkpoint/figures/phase_diagram3.png')