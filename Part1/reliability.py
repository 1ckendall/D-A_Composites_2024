import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
from pprint import pprint
from Class import Lamina, Laminate

# plot controls (Booleans)
plot_1a = True # plot 1a: engineering constants for laminae (default = True)
plot_1b = False # plot 1b: stress-strain (default = True)

# savefig controls
savefig = False # default = False

# functions
def generate_cdf(mean, std_dev, num_points=int(1e4), is_sorted = False):  
    # num points: use 1e4 during developement, use atleast 1e5 for final simulation
    # Generate random samples from a normal distribution
    samples = np.random.normal(mean, std_dev, num_points)
    
    if is_sorted:
        samples = np.sort(samples)
    cdf = np.arange(1, len(samples) + 1) / len(samples)
    
    return samples, cdf

# UD Lamina Mean Properties 
E1_mean = 145.3E9 # [Pa]
E2_mean = 8.5E9 # [Pa]
v12_mean = 0.31 # [--]
G12_mean = 4.58E9 # [Pa]
Xt_mean = 1932E6 # [Pa]
Yt_mean = 108E6 # [Pa]
Xc_mean = 1480E6 # [Pa]
Yc_mean = 220E6 # [Pa]
S_mean = 132.8E6 # [Pa]
t_mean = 0.125E-3 # [m] # free variable, can be changed 

UD_mean = np.array((E1_mean, E2_mean, v12_mean, G12_mean, Xt_mean, Yt_mean, Xc_mean, Yc_mean, S_mean, t_mean))


# UD Lamina Std. Dev Properties 
E1_std = 3.28E9 # [Pa]
E2_std = 1.28E9 # [Pa]
v12_std = 0.018 # [--]
G12_std = 0.83E9 # [Pa]
Xt_std = 128.3E6 # [Pa]
Yt_std = 8.2E6 # [Pa]
Xc_std = 98.28E6 # [Pa] # assumed, same fraction as: Xc_std = Xc_mean * Xt_std/Xt_mean
Yc_std = 16.70E6 # [Pa]  # assumed, same fraction as: Yc_std = Yc_mean * Yt_std/Yt_mean
S_std = 6.21E6 # [Pa] 
t_std = 0.002E-3 # [m] # assumed, from certificate of analysis of a pre-preg material (from Solvay)

UD_std = np.array((E1_std, E2_std, v12_std, G12_std, Xt_std, Yt_std, Xc_std, Yc_std, S_std, t_std))


# Question 3: Reliability Analysis

layup =  [0, 90, +45, -45, - 45, + 45, 90, 0, 0, 90, +45, -45, - 45, + 45, 90, 0] # laminate [0/90/±45]_2s
N_load = 850 # 1250 [N], applied load
theta = 30 # [deg], inclination of the load vector w.r.t. x-axis
Nx = N_load * np.cos(np.radians(theta)) # [N]
Ny = N_load * np.cos(np.radians(theta)) # [N]

n_vars = len(UD_mean) # number of independent, Gaussian random variables 
iterations = 10000 # 1E8 Monte Carlo: number of rounds of simulations (R)
Pf_arr = np.zeros((n_vars, iterations)) # array for probabability of failure, registered for each round of simulations and for every random variable


# simulation loop
# loop through all random variables 
for i in range(n_vars):
    samples, cdf = generate_cdf(mean = UD_mean[i], std_dev = UD_std[i])
    # loop through rounds of simulations (R)
    for j in range(iterations): 
        failed = False
        N = 0 
        # determine number of simulations (N) needed for failure 
        while not failed:
            pass # Puck failure criterion
            N += 1
        Pf_arr[i,j] = 1/N
        
        
    




