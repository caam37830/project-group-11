# Agent-Based and Ordinary Differential Equation Modelling of SIR Model -- With Extended Functionality

## Abstract

## 1. Introduction to the SIR Model
The SIR model is a compartmental epidemiological model which assumes everyone in a population falls into one of three compartments: Susceptible, Infectious, and Recovered (sometimes referred to instead as Removed). The Susceptible population consists of those who have not yet gotten the virus and thus are _susceptible_ to becoming infected. Infectious individuals are those who carry the virus and are capable of infecting others for some period of time before eventually transitioning into the Removed group. Lastly, Removed individuals are those who are no longer part of the system because they have either acquired immunity or died from the disease. 

The trajectory of the disease generally depends on the parameters `b`, the per capita number of interactions capable of spreading the disease that occur per day, and `k`, the rate at which Infected individuals move into the Removed compartment. Notably, the basic SIR model assumes that all agents in the population interact randomly so every Susceptible individual shares an equal chance of becoming infected in each period. I.e. the probability for any given Susceptible individual to become infected in a period is `b * I/N`. Additionally, `k` can be interpreted as the reciprocal of the mean duration of infection when duration of infections are exponentially distributed (1). For Covid-19, no conclusive figure has been established as the mean duration of infection. However, the US CDC estimates that the majority of cases become no longer infectious within 10 days, though some serious cases can remain infectious for up to 20 days (2). 

From the dicussion above, we can note that for any period, we can predict if the infection will continue to spread or being to decrease in the next period. Since the (expected) number of new infections can be given by `S * b * I/N` and the (expected) number of removals is given by `k * I`, if the ratio of the two, `(S * b) / (k*N) > 1` then the infection will continue to spread since the number of new infections outpaces the number of recoveries; whereas if `(S * b) / k*N) < 1`, the number of recoveries exceed te number of infections and so the infection will begin to die out on its own. Within epidemiolgy, this relationship is commonly known as the 'effective reproduction number' (3).

## 2. Spatial SIR Model

### The Spatial PDE Model
Firstly, let's look at the spatial SIR model implemented by updating our ODE model to a similar PDE model. The system of PDEs is:

    \partial s(x, t)/ \partial t = -b * s(x, t) * i(x, t) + p * L * s(x, t)
    \partial r(x, t)/ \partial t = k * i(x, t) + p * L * r(x, t)
    \partial i(x, t)/ \partial t = b * s(x, t) * i(x, t) - k * i(x, t) + p * L * i(x, t)

Here we have a new parameter `p` which is to weight the diffusion term. A larger `p` means each compartment diffuse quickly in space. 

Set `b=3`,`k=0.1` test different p values.

p=0.02|p=0.1|p=1|p=10
-|-|-|-
![](figures/pde_p0.02.png)|![](figures/pde_p0.1_center.png)|![](figures/pde_p1.png)|![](figures/pde_p10.png)

From above plots, we can find that recovered population of all plots approximate to 100% when t goes large. 
This means for different `p` all population will be infected at last. 
However, the ratio of total infected population (recovered and infectious compartments) grows up at different speed. 
When p is large, the trajectory of spatial SIR model is similar to SIR model. 
This is because large p means a person can interact with more people, 
so when p is large enough the trajectory of spatial SIR model would be very similar to that of a non-spatial SIR model.

Therefore choosing a relative small p will be interesting for spatial SIR model. I select `p=0.1`
By setting i(x,0) appropriately, we can investigate the location effect of initial infected individuals.

- set i((0,0), 0) to be nonzero but other i to be zeros, we can simulate diseases spreading started in a single corner of the square.
- set i((M/2,M/2), 0) to be nonzero but other i to be zeros, we can simulate diseases spreading started in the center of the square.
- randomly select several points of i to be nonzero but other i to be zeros, we can simulate random spread out situation.

corner|center|random
-|-|-
![](figures/pde_p0.1_corner.png)|![](figures/pde_p0.1_center.png)|![](figures/pde_p0.1_random.png)

From the figures above, we can find that

- if the disease start at a corner of the square, it will take a longer time for the disease to spread out. The ratio of illness population at one time will be small or the curve be "smoothed". This `square world` will a longer time to prepare.
- if the disease start at the center of the square, it will take a middle time for the disease to spread out and the the curve is higher than the corner one. 
- if the disease start at some random points of the square, this is fastest situation.

### The Spatial Agent-Based Model
In addition to our PDE model, we also implemented a spatial component into our existing agent-based simulations. Upon initialization, agents in our spatial simulations were given two dimensional coordinates to designate their starting position contained within the unit square.

The daily steps of the simulation remain largely the same, with a specified proportion (`k`) of those infected being removed with the remaining infected interacting, and possibly infecting, the rest of the populace. Rather than designating a set number of interactions for the infected, each agent in the spatial simulation will move once in a random direction each day and interact with whichever agents happen to be where they end up in the unit square. Infected agents will infect whoever they reach and susceptible agents will get infected if they end up near a diseased agent. Recovered agents will not move at all as they cannot influence any part of the simulation in that state.

Daily movements will be scaled by a new parameter (`p`) and their net of interactions will be determined by the circle with radius `q` which will effectively replace our previous parameter `b`. Any daily movements that would see the agent exceed the bounds of the unit square will have them instead remain in the same location for the day.

Based on the earlier SIR model, we've elected to choose `k` = 0.1 and `q` $\approx$ 0.02 to determine the effect of `p` on simulations. The latter was determined based on the below rough relationship with `b`

 $q = \sqrt(\frac{b}{(N * \pi)}$

The below simulations started with an infection in the center of the unit square.

p=0.02|p=0.10|p=0.40|p=0.75
-|-|-|-
![](figures/abmStepSize0.02.png)|![](figures/abmStepSize0.10.png)|![](figures/abmStepSize0.40.png)|![](figures/abmStepSize0.75.png)

Similar to our PDE model, we see that no value of `p` is completely preventative of forever delaying a completely infected populace. However, we can note a few differences in the agent-based simulation's relationship with `p`.

Our slowest growing diseases come from our lowest and highest values of `p`. This comes from the infected simply making very meagre progress across the unit square, or from making such large moves that they exceed the boundaries of the square and do not end up moving. A middle ground sees infected able to cover large parts of the unit square without often stepping out of bounds.

To identify the relationship between our spatial component and the duration of our epidemic, we then can look at different simulations where the infection starts in a Corner, the Center, or at random points across the unit square. For a more moderate virus expansion, we've selected a step size of `p` = 0.02.


corner|center|random
-|-|-
![](figures/abm2dCorner.png)|![](figures/abm2dCenter.png)|![](figures/abm2dCRndom.png)

The agent-based figures largely compare well to the PDE counterparts:

- Corner-based diseases are the slowest to perpetuate across the unit square. Even with a small step size, infected agents are very likely to overstep the two edges they start next to which could greatly impede disease progress.
- Center-based diseases are faster, simply as they are unlikely to overstep any edges of the square at such a low step size. However, having all of the infected start in one location delays the spread a little bit.
- Randomly initialized diseases have the easiest time spreading as the infection can start in different `communities` of the square from the beginning.

## 3. Extensions to the SIR Model

### Allowing for different rates of contact between certain individuals

### Incorporate Interventions
For ODE model, use the most general model.

![](./figures/MSEIQRDS.png)

Q is not a new compartment, it only means `q` percent of `I` will not spread the virus, because they quarantine themselves.

    [maternally_derived_immunity-susceptible-exposed-infectious-quarantine-recovered-decreased]
    dm / dt = lam - sigma * m - mu * m
    ds / dt = sigama * m + re * r - mu * s - b * s * (1-q)i
    de / dt = b * s * (1-q)i - (mu + a) * e
    di / dt = a * e - (k + mu) * i - d * log(i/1-i)
    dd / dt = d * log(i/1-i)
    dr / dt = k * i - (mu +re) * r

If the society does not take any action, the trajectory will look like this. Parameters are set to be

(lam=3e-05, sigma=1/720, b=3, k=1/10, q=0, a=1/14, mu=3e-05, dr=0.3, re=1/360)

![](./figures/no_intervention.png)

By define some events functions, we can implement ODE with interventions. We will compare with the figure above with the trajectory of interventions.
    
#### mask
![](/figures/mask1.png)


For ODE model, wearing a mask will decrease `b` of the model. Therefore, we define a event when `i=0.05`, people are asked to wear masks and fortunately most people follow this advice `b` change from 3 to 0.1 .

wearing a mask|no intervention
-|-
![](figures/intervention_b.png)| ![](figures/no_intervention.png)

the dashed line shows when people follow this advice, before and after the dashed line shows the effect of the intervention of wearing a mask.
Compare two plots, we can find that wearing a mask significantly reduce the ratio of disease caused death population. Most people will not suffer the pain of illness and remain to be susceptible. 

#### vaccine

#### quarantine


For ODE model, quarantine will limit `I`'s movement, reduce the interaction of `S` and `I`. Here we define two events, when `i=0.05`, `q=0.8`, when `i=0.1`, `q=0.99`.

suggest quarantine | no intervention
-|-
![](figures/intervention_q.png) | ![](figures/no_intervention.png)

From above plots, we can find that by taking interventions total death population reduced significantly. 

### Additional Compartments by ODE method
In this part we extend SIR model to various other ODE models and have a look at the plot of those models.
For detailed model description please check in appendices. Here we only show the result.

![](figures/ode_sir.png)|![](figures/ode_sis.png)
![](figures/ode_sird.png)|![](figures/ode_msir.png)
![](figures/ode_sirc.png)|![](figures/ode_seir.png)
![](figures/ode_seis.png)|![](figures/ode_mseir.png)
![](figures/ode_mseirs.png)|![](figures/ode_mseiqrds.png)

## 3. Discussion and Conclusion


## 4. Bibliography
