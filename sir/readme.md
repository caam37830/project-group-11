This is the description file for this sir package. Main python models like SIR-ODE and SIR-ABM are in this package.

## ode models
In [ode.py](ode.py) file we have various ode model for desease spreading. We have,
- [SIR model](ode.py#L147-L190)
- [SIS model](ode.py#L193-L234)
- [SIRD model](ode.py#L237-L290)
- [MSIR model](ode.py)
- [SIRC model](ode.py)
- [SEIR model](ode.py)
- [SEIS model](ode.py)
- [MSEIR](ode.py)
- [MSEIRS](ode.py)



### 1. how to use ode model?
Here are several lines to illustrate it.

This line we model a SIS model with some parameter,
for more detailed parameter description, please check the source code sir/ode.py
    
    ode_SIS = SIS(b=0.8, k=0.01)
    
This line we give a t_bound to simulate the system from day0 to day t_bound
    
    ode_SIS.solve(t_bound=300)
    
This line we plot the simulation result, if you want to save the plot, please plugin a path
    
    ode_SIS.plot()

With only three lines of code, you can simulate a ode model for covid-19
we have following examples for more information of this ode system,
please check the source code to find more
    
    # after we solve this model, we can find any value between [0, t_bound]
    ode_SIS(t=20)  
    # this return ode.OdeResult object, including t, y, sol, etc.
    ode_SIS.sol  
    
### 2. all models
    
    from sir.ode import *

#### 2.1 SIR model
    
    ode_sir = SIR(b, k, **kwargs)
        
        ODE Model
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
    
#### 2.2 SIS model

    ode_sis = SIS(b, k, **kwargs)
        
        ODE Model
            ds / dt = - b * s * i + k * i
            di / dt = b * s * i - k * i
            
        Parameter
            b - b is number of interactions per individual per day
            k - k is fraction of infectious period which recovers each day (0 < k < 1)

        Optional Parameter
            S - susceptible population
            I - infectious population
            N - total population

#### 2.3 SIRD model

    ode_sird = SIRD(b, k, mu, **kwargs)
    
        ODE model
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

#### 2.4 MSIR model

    ode_msir(lam, sigma, b, k, mu, **kwargs)
    
        ODE model
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

#### 2.5 SIRC model

    ode_sirc(b, k, c1, c2, **kwargs)
    
        ODE model [susceptible-infectious-recovered-carrier]
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

#### 2.6 SEIR model

    ode_seir(lam, b, k, a, mu, **kwargs)
    
        ODE model [susceptible-exposed-infectious-recovered]
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

#### 2.7 SEIS model

    ode_seis(lam, b, k, a, mu, **kwargs)
    
        ODE model [susceptible-exposed-infectious-susceptible]
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

#### 2.8 MSEIR model

    ode_mseir(lam, sigma, b, k, a, mu, **kwargs)
    
        ODE model [maternally_derived_immunity-susceptible-exposed-infectious-recovered]
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

#### 2.9 MSEIRS model
    
    ode_mseirs(lam, sigma, b, k, a, mu, l, **kwargs)
    
        ODE model [maternally_derived_immunity-susceptible-exposed-infectious-recovered-susceptible]
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
