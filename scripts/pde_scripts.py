from sir.pde import *

# test different p
p_list = [1, 10]
for p in p_list:
    # init model, covid source from middle
    model = covid2D(p=p, source_loc=np.array([[99], [99]]), M=200)
    model.solve(t_bound=365)
    model.tracePlot('./doc/final/figures/pde_p%f.png' % p)
    del model

# test different p
p_list = [0.1, 0.3, 0.5, 0.7, 0.9]
for p in p_list:
    model = covid2D(p=p, source_loc=np.array([[99], [99]]), M=200)
    model.solve(t_bound=365)
    model.tracePlot('./doc/final/figures/pde_p%.1f.png' % p)
    del model

# test different p
p_list = [0.02, 0.04, 0.06, 0.08]
for p in p_list:
    model = covid2D(p=p, source_loc=np.array([[99], [99]]), M=200)
    model.solve(t_bound=365)
    model.tracePlot('./doc/final/figures/pde_p%.2f.png' % p)
    del model

# choose p=0.1
p = 0.1

# init model with covid source from the corner
model = covid2D(p=p, source_loc=np.array([[0], [0]]), M=200)
model.solve(t_bound=365)
model.tracePlot('./doc/final/figures/pde_p%.1f_corner.png' % p)
del model

# init model with covid source from uniform
model = covid2D(p=p, sourcepoint=None, M=200)
model.solve(t_bound=365)
model.tracePlot('./doc/final/figures/pde_p%.1f_unif.png' % p)
del model

# init model with covid source from random points
model = covid2D(p=p, sourcepoint=5, M=200)
model.solve(t_bound=365)
model.tracePlot('./doc/final/figures/pde_p%.1f_random.png' % p)
del model
