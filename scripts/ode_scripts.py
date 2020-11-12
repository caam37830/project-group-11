from sir.ode import covid, phase_plot

# init population
N = 100_000
I = 100
R = 0
S = N - I - R

# init the covid class
ode_covid = covid(S, I, R, b=3, k=0.1)
# solve that class
ode_covid.solve(t_bound=100)
# plot the numerical solution
ode_covid.plot(save_path=None)
# s, i, r of day 150
print(ode_covid(100))


# init the covid class
ode_covid = covid(S, I, R, b=3, k=0.1)
ode_covid.solve(t_bound=100)
ode_covid.plot(save_path=None)

# init the covid class
ode_covid = covid(S, I, R, b=0.8, k=0.1)
ode_covid.solve(t_bound=100)
ode_covid.plot(save_path=None)

# Phase diagram with b, k axes with portion of population that is infected at time t 
phase_plot(N, I, R, t=5, phase='I', save_path=None)
phase_plot(N, I, R, t=10, phase='I', save_path=None)
phase_plot(N, I, R, t=50, phase='I', save_path=None)