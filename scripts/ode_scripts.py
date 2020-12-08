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

# other models
from sir.ode import *
ode_SIS = SIS(b=3, k=0.1)
ode_SIS.solve(t_bound=365)
ode_SIS.plot('./doc/final/figures/ode_sis.png')

ode_SIRD = SIRD(b=3, k=0.1, mu=0.05)
ode_SIRD.solve(t_bound=365)
ode_SIRD.plot('./doc/final/figures/ode_sird.png')

ode_MSIR = MSIR(lam=0.0003, sigma=1/180, b=3, k=0.1, mu=0.0003)
ode_MSIR.solve(t_bound=365)
ode_MSIR.plot('./doc/final/figures/ode_msir.png')

ode_SIRC = SIRC(b=3, k=0.1, c1=0.1, c2=0.1)
ode_SIRC.solve(t_bound=365)
ode_SIRC.plot('./doc/final/figures/ode_sirc.png')

ode_SEIR = SEIR(lam=0.0003, b=3, k=0.1, a=1/14, mu=0.0003)
ode_SEIR.solve(t_bound=365)
ode_SEIR.plot('./doc/final/figures/ode_seir.png')

ode_SEIS = SEIS(lam=0.0003, b=3, k=0.1, a=1/14, mu=0.0003)
ode_SEIS.solve(t_bound=365)
ode_SEIS.plot('./doc/final/figures/ode_seis.png')

ode_MSEIR = MSEIR(lam=0.0003, sigma=1/180, b=3, k=0.1, a=1/14, mu=0.0003)
ode_MSEIR.solve(t_bound=365)
ode_MSEIR.plot('./doc/final/figures/ode_mseir.png')

ode_MSEIRS = MSEIRS(lam=0.0003, sigma=1/180, b=3, k=0.1, a=1/14, mu=0.0003, l=1/180)
ode_MSEIRS.solve(t_bound=365)
ode_MSEIRS.plot('./doc/final/figures/ode_mseirs.png')

ode_MSEIQRDS = MSEIQRDS()
ode_MSEIQRDS.solve(t_bound=365)
ode_MSEIQRDS.plot('./doc/final/figures/ode_mseiqrds.png')
