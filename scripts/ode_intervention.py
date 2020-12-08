from sir.ode import *
def iEq1(t, y):
    return y[3] - 0.05

# terminate if y=0 and the direction is from positive to negetive.
iEq1.terminal = True
iEq1.direction = 1
iEq1.intervention = {'q': 0.8}

def iEq2(t, y):
    return y[3] - 0.1

# terminate if y=0 and the direction is from positive to negetive.
iEq2.terminal = True
iEq2.direction = 1
iEq2.intervention = {'q': 0.99}

model = MSEIQRDS()
sol = intervention_solve(model, events=[iEq1, iEq2], t_bound=365, h=1)

intervention_plot(sol)
