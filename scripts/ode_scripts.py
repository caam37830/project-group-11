from sir.ode import covid

# init population
N = 1_000_000
I = 5
R = 0
S = N - I - R

# init the covid class
ode_covid = covid(S, I, R, k=0.05, b=3)
# solve that class
ode_covid.solve(t_bound=365)
# plot the numerical solution
ode_covid.plot()
# s, i, r of day 150
print(ode_covid(100))
