import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
from pprint import pprint
from Class import Lamina, Laminate

import time

start_time = time.time()

# plot controls (Booleans)

# savefig controls
savefig = False # default = False

# functions
def generate_cdf(mean, std_dev, num_points=int(1e5), is_sorted = False):  
    """num points: 
        use 1e4 during developement:  
        use atleast 1e5 for final simulation: 
        
        COMPUTATIONAL TIME:
            1e4: ~15s
            1e5: ~30s
            1e6: ~170s
    """
    distance = norm(loc=mean, scale=std_dev)
    cdf = np.linspace(mean - 3*std_dev, mean + 3*std_dev, num_points)
    samples = distance.cdf(cdf)
    return samples, cdf

def inverse_transform_sampling(samples, cdf, num_points=1):
    # Generate random numbers uniformly distributed between 0 and 1
    u = np.random.rand(num_points)
    
    #sort the samples for interpolation
    # sorted_samples = np.sort(samples)
    sorted_samples = samples

    # Use inverse transform sampling to find corresponding points from the CDF
    interpolated_samples = np.interp(u, cdf, sorted_samples)
    
    return interpolated_samples

def gaussian_distribution(mu=0, sigma=1, num_points = 1000):
    x_arr = np.linspace(mu - 3*sigma, mu + 3*sigma, num_points)  # Generate x values
    gaussian_arr = (1/(sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x_arr - mu) / sigma)**2)
    return x_arr, gaussian_arr


# UD Lamina Material Properties
UD_names = np.array(['E1', 'E2', 'v12', 'G12', 'Xt', 'Xc', 'Yt', 'Yc', 'S', 't'])
n_vars_active = np.array([True, True, True, True, True, False, True, False, True, False])
count_n_vars_active = np.count_nonzero(n_vars_active)

'''
Index of each property variable in UD array
Property | Index
------------------
E1       | 0
E2       | 1
v12      | 2
G12      | 3
Xt       | 4
Xc       | 5
Yt       | 6
Yc       | 7
S        | 8
t        | 9       
'''

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
# t_std = 0.002E-3 # [m] # TOO HIGH (OVERSENSITIVE RESULTS) assumed, from certificate of analysis of a pre-preg material (from Solvay)
t_std = 0.002E-6 # [m] # assumed and adjusted
 
UD_std = np.array((E1_std, E2_std, v12_std, G12_std, Xt_std, Yt_std, Xc_std, Yc_std, S_std, t_std))

# Question 3: Reliability Analysis

layup =  [0, 90, +45, -45, - 45, + 45, 90, 0, 0, 90, +45, -45, - 45, + 45, 90, 0] # laminate [0/90/±45]_2s

N_load = 0.1e6 # (pf = 0.7142857142857143)


theta = 30 # [deg], inclination of the load vector w.r.t. x-axis
Nx = N_load * np.cos(np.radians(theta)) # [N/m]
Ny = N_load * np.sin(np.radians(theta)) # [N/m]

n_vars = len(UD_mean) # number of independent, Gaussian random variables 
iterations = 10 # 3 # 1E8 Monte Carlo: number of rounds of simulations (R)
N_max = 200 # maximum number of simulations per round (N)
Pf_arr = np.zeros((n_vars, iterations)) # array for probabability of failure, registered for each round of simulations and for every random variable

# simulation loop
UD = UD_mean
loop_counter = 0 
# loop through all random variables 
#for i in range(n_vars): # full simulation (without thickness)
for i in range(n_vars): # full simulation
    samples, cdf = generate_cdf(mean = UD_mean[i], std_dev = UD_std[i])
    for j in range(iterations): 
        firstplyfailureoccurence = False
        damageoccurence = False
        N = 0
        # determine number of simulations (N) needed for failure 
        while firstplyfailureoccurence == False and N < (N_max):
            N += 1
            loop_counter += 1
            # sample random variable from CDF
            sample = inverse_transform_sampling(samples, cdf, num_points=1) #TODO: optimise this (it is in the third layer of the loop) eg: use scipy, or static inline
            #print(f'i = {i}, Sampled value: {sample}) #TOCHECK (bias in sampling)
            UD[i] = sample 
            # instantiate a Laminate object
            plylist_ijk = [] 
            for angle in layup:
                plylist_ijk.append(Lamina(angle, *UD))
            Laminate_ijk = Laminate(plylist_ijk, Nx=Nx, Ny=Ny, Ns=0, Mx=0, My=0, Ms=0)
            
            Laminate_ijk.getStressStrain()
            
            # failuretracking = 0 
            # checking for failure per lamina. Output: update boolean 'firstplyfailureoccurence'
            for idx, ply in enumerate(plylist_ijk):
                failuretracking = 0 
                # assign computed stresses as lamina attributes
                ply.sigma_1 = Laminate_ijk.sigma11[idx]
                ply.sigma_2 = Laminate_ijk.sigma22[idx]
                ply.tau_21 = Laminate_ijk.sigma12[idx]
                
                # Puck failure criterion
                ply.PuckFibreFail(sigma_1 = ply.sigma_1, sigma_2 = ply.sigma_2, sigma_3 = ply.tau_21, R_para_t = UD[4], R_para_c = UD[5], v_perppara = UD[2], E_para = UD[0])  #Puck fiber failure #### currently using Minor poisson, to make it consistent with Minor poisson at failure
                ply.PuckIFF(sigma_22 = ply.sigma_2, sigma_21 = ply.tau_21 , sigma_22T_u = UD[6], sigma_22C_u=UD[7], sigma_12_u=UD[8])  #Puck fiber failure #### currently using Minor poisson, to make it consistent with Minor poisson at failure          
                
                if ply.failuremode == 'FFT' or ply.failuremode == "FFC" or ply.failuremode == 'IFF A' or ply.failuremode == "IFF B" or ply.failuremode == "IFF C":
                    failuretracking = 2
                    firstplyfailureoccurence = True
                    Pf_arr[i,j] = 1/N
            print(f'loop count: {loop_counter}, i = {i}, j = {j}, k = {N}, isFPF?: {firstplyfailureoccurence} ')

Pf_vars_mean = np.zeros(n_vars)
Pf_vars_std = np.zeros(n_vars)
Pf_vars_mean_norm = np.zeros(n_vars)
Pf_vars_std_norm = np.zeros(n_vars)

for n in range(n_vars):
    Pf_vars_mean[n] = np.mean(Pf_arr[n,:])
    Pf_vars_std[n] = np.std(Pf_arr[n,:])

Pf_vars_mean_norm = Pf_vars_mean-Pf_vars_mean
Pf_vars_std_norm = Pf_vars_std/np.sqrt(iterations)
    
Pf_mean = np.mean(Pf_arr)
Pf_std = np.std(Pf_arr)
print(f'\nLoad: {N_load}[N/m] => Probability of Failure: {Pf_mean}')
print(f'\nIterations (R): {iterations}, Nmax: {N_max} =>  Std. Dev: {Pf_std}')
print(f'\nLoad: {N_load}[N/m] => Probability of Failure: {Pf_mean}, Std. dev: {Pf_std}')


# Record end time
end_time = time.time()

# Compute elapsed time
elapsed_time = end_time - start_time

print("Time taken:", elapsed_time, "seconds")


# convergence: generate Gaussian distribution
exact_x_arr, exact_gaussian_arr = gaussian_distribution()

# normalise (mu, sigma) -> (0,1)
Pf_mean_norm = Pf_mean-Pf_mean  # 0 
Pf_std_norm = Pf_std/(np.sqrt(iterations*n_vars))
# print(f'Normalised, Gaussian Monte-Carlo Variable (Probability of Failure): mean = {Pf_mean_norm}, std. = {Pf_std_norm}')

# error bar values w/ different -/+ errors that
# also vary with the x-position
lower_error =  np.zeros(n_vars)
upper_error =  np.zeros(n_vars)
for i in range(n_vars):
    if Pf_mean>Pf_vars_mean[i]:
        lower_error[i] = Pf_mean-Pf_vars_mean[i]
    elif Pf_mean<=Pf_vars_mean[i]:
        upper_error[i] = Pf_vars_mean[i]-Pf_mean
asymmetric_error = np.array(list(zip(lower_error, upper_error))).T

#filename = f'Monte_Carlo_Puck_Load={N_load/10**3}_Nslmm_R={iterations}_Nmax={N_max}'
figname = f'Monte Carlo (Puck Criterion): Load = {N_load/10**3}N/mm, (R = {iterations}, Nmax = {N_max})'
plt.figure('1')
plt.clf()
plt.title(figname)
plt.errorbar(np.arange(n_vars), Pf_mean*np.ones(n_vars), asymmetric_error, color = 'blue', fmt='o', ecolor='red', capsize=6)
plt.axhline(Pf_mean, color='blue', linestyle='--', label=f'Probability of Failure: {Pf_mean}')
plt.xticks(np.arange(n_vars), UD_names, fontsize = 16)
plt.yticks(fontsize = 16)
plt.legend(fontsize = 16)
plt.grid(True,alpha=0.5)
plt.ylabel(r"Probability of Failure", fontsize = 16)


# figname = f'Monte Carlo Convergence (Central Limit Theorem): Load = {N_load/10**3}N/mm, (R = {iterations}, Nmax = {N_max})'
# plt.figure('2')
# plt.clf()
# plt.title(figname)
# # plt.plot(angleinputvector, damage_tol_Puck, linewidth=1, color='blue', linestyle='--', label='Puck')
# plt.axhline(Pf_mean, color='blue', linestyle='--', label=f'Probability of Failure: {Pf_mean}')
# plt.xticks(np.arange(n_vars), UD_names, fontsize = 16)
# plt.yticks(fontsize = 16)
# plt.legend(fontsize = 16)
# plt.grid(True,alpha=0.5)
# plt.ylabel(r"Probability of Failure", fontsize = 16)

plt.show()
