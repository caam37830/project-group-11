import sys
sys.path.append("../")
from sir.ode import *

# init population
N = 1_000_000
I = 5
R = 0
S = N - I - R

# init the covid class
ode_SIR = SIR(b=3, k=0.01, S=S, I=I, R=R)
# solve that class
ode_SIR.solve(t_bound=365)
# plot the numerical solution
ode_SIR.plot('../doc/checkpoint/figures/ode1.png')
# s, i, r of day 150
print(ode_SIR(100))


# init the sir class with different parameters
ode_SIR = SIR(b=3, k=0.1, S=S, I=I, R=R)
ode_SIR.solve(t_bound=365)
ode_SIR.plot('../doc/checkpoint/figures/ode2.png')

# init the sir class with different parameters
ode_SIR = SIR(b=0.8, k=0.01, S=S, I=I, R=R)
ode_SIR.solve(t_bound=365)
ode_SIR.plot('../doc/checkpoint/figures/ode3.png')

# Phase diagram with b, k axes with portion of population that is infected at time t 
phase_plot(N, I, R, t=5, phase='I', save_path='../doc/checkpoint/figures/phase_diagram1.png')
phase_plot(N, I, R, t=10, phase='I', save_path='../doc/checkpoint/figures/phase_diagram2.png')
phase_plot(N, I, R, t=50, phase='I', save_path='../doc/checkpoint/figures/phase_diagram3.png')
