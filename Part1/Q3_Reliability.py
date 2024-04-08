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
def generate_cdf(mean, std_dev, num_points=int(1e4), is_sorted = False):  
    # num points: use 1e4 during developement, use atleast 1e5 for final simulation
    # Generate random samples from a normal distribution
    samples = np.random.normal(mean, std_dev, num_points)
    
    if is_sorted:
        samples = np.sort(samples)
    cdf = np.arange(1, len(samples) + 1) / len(samples)
    
    return samples, cdf

def inverse_transform_sampling(samples, cdf, num_points=1):
    # Generate random numbers uniformly distributed between 0 and 1
    u = np.random.rand(num_points)
    
    #sort the samples for interpolation
    sorted_samples = np.sort(samples)

    # Use inverse transform sampling to find corresponding points from the CDF
    interpolated_samples = np.interp(u, cdf, sorted_samples)
    
    return interpolated_samples

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
t_std = 0.002E-3 # [m] # assumed, from certificate of analysis of a pre-preg material (from Solvay)

UD_std = np.array((E1_std, E2_std, v12_std, G12_std, Xt_std, Yt_std, Xc_std, Yc_std, S_std, t_std))


# Question 3: Reliability Analysis

layup =  [0, 90, +45, -45, - 45, + 45, 90, 0, 0, 90, +45, -45, - 45, + 45, 90, 0] # laminate [0/90/±45]_2s
# benchmark N_load for MaxStress
# N_load = 1.20e6 # 1.5e6 (pf = 0.0) # 2e6 (pf = 0.0) # 3e6 (pf = 0.31027094717668485) # 3.4e6 # (pf = 0.3610857) # 3.4375e6 (pf = 0.9715166666666667)  # 1.2e6 [N/m] (pf = 0.0003) # 0.85e6 [N/m] (pf = 0), applied load # INCORRECT AS ALL FAILURE MODES (INCLD. DAMAGE) WERE CONSIDERED FAILURE

# benchmark N_load for Puck
N_load = 0.6e6  # 0.6e6 (pf =  0.35666666666666663) # 0.7e6 (pf =  1.0) # 0.5e6 (pf =  0.0) # 0.1e6 (pf =  0.0) # 1e6 (pf =  1.0) # 2e6 (pf =  1.0) # 3e6 (pf =  1.0) # 1.2e6 [N/m] (pf = 0.0003) # 0.85e6 [N/m] (pf = 0), applied load

theta = 30 # [deg], inclination of the load vector w.r.t. x-axis
Nx = N_load * np.cos(np.radians(theta)) # [N]
Ny = N_load * np.cos(np.radians(theta)) # [N]

n_vars = len(UD_mean) # number of independent, Gaussian random variables 
iterations = 3 # 1E8 Monte Carlo: number of rounds of simulations (R)
N_max = 100 # maximum number of simulations per round (N)
Pf_arr = np.zeros((n_vars, iterations)) # array for probabability of failure, registered for each round of simulations and for every random variable

# simulation loop
UD = UD_mean
loop_counter = 0 
# loop through all random variables 
for i in range(n_vars): # full simulation
#for i in range(0, 1): # testing (only treat E1 as random variable)
    samples, cdf = generate_cdf(mean = UD_mean[i], std_dev = UD_std[i])#, num_points = iterations)
    # loop through rounds of simulations (R)
    for j in range(iterations): 
        firstplyfailureoccurence = False
        N = 0
        # determine number of simulations (N) needed for failure 
        while firstplyfailureoccurence == False and N < (N_max):
            N += 1
            loop_counter += 1
            # sample random variable from CDF
            sample = inverse_transform_sampling(samples, cdf, num_points=1) #TODO: optimise this (it is in the third layer of the loop) eg: use scipy, or static inline
            UD[i] = sample 
            # instantiate a Laminate object
            plylist_ijk = [] 
            for angle in layup:
                # plylist_ijk.append(Lamina(angle,E1,E2,G12,v12,Xt,Xc,Yt,Yc,S,t))
                plylist_ijk.append(Lamina(angle, *UD))
            Laminate_ijk = Laminate(plylist_ijk, Nx=Nx, Ny=Ny, Ns=0, Mx=0, My=0, Ms=0)
            
            Laminate_ijk.getStressStrain()
            
            failuretracking = 0 
            # checking for failure per lamina. Output: update boolean 'firstplyfailureoccurence'
            for idx, ply in enumerate(plylist_ijk):
                # assign computed stresses as lamina attributes
                ply.sigma_1 = Laminate_ijk.sigma11[idx]
                ply.sigma_2 = Laminate_ijk.sigma22[idx]
                ply.tau_21 = Laminate_ijk.sigma12[idx]
                
                # print(f'Laminate_ijk.sigma11: {Laminate_ijk.sigma11[idx]}, Laminate_ijk.sigma22: {Laminate_ijk.sigma22[idx]}, Laminate_ijk.sigma12: {Laminate_ijk.sigma12[idx]}')
                # print(f'ply sigma_1: {ply.sigma_1}, ply.sigma_2: {ply.sigma_2}, ply.tau_21: {ply.tau_21}')
                
                # # Maximum stress failure criterion (Placeholder)
                # ply.maxStressFibreFail() #fiber failure
                # ply.maxStressInterFibreFail() #InterFiberfailure
                # ply.maxStressShearFail() #Shearfailure
                
                # if ply.failuremode == "Tensile Fibre" or ply.failuremode == "Compressive":
                #   failuretracking = 2
                # elif ply.failuremode == "Tensile Inter-Fibre" or ply.failuremode=="Compressive Inter-Fibre" or ply.failuremode ==  "Shear Parallel to Fibres":
                #   failuretracking += 1
                
                
                
                # Puck failure criterion
                
                ply.PuckFibreFail(sigma_1 = ply.sigma_1, sigma_2 = ply.sigma_2, sigma_3 = ply.tau_21, R_para_t = UD[4], R_para_c = UD[5], v_perppara = UD[2], E_para = UD[0])  #Puck fiber failure #### currently using Minor poisson, to make it consistent with Minor poisson at failure
                ply.PuckIFF(sigma_22 = ply.sigma_2, sigma_21 = ply.tau_21 , sigma_22T_u = UD[6], sigma_22C_u=UD[7], sigma_12_u=UD[8])  #Puck fiber failure #### currently using Minor poisson, to make it consistent with Minor poisson at failure          
                
                if  ply.failuremode == "FFT" or ply.failuremode == "FFC":
                    failuretracking = 2
                elif ply.failuremode == 'IFF A' or ply.failuremode == "IFF B" or ply.failuremode == "IFF C":
                    failuretracking += 1
                
                
                  # Puck failure criterion
            #     ply.PuckFibreFail(sigma_2, gamma_21, epsilon1T, epsilon1C, epsilon1)
            #     RAperpperp = ply.get_RAperpperp(S21, p_perppara_minus, Yc)
            #     p_perpperp_minus = ply.get_p_perpperp_minus(p_perppara_minus, RAperpperp, S21)
            #     tau_21c = ply.get_tau_21c(S21, p_perpperp_minus)
            #     ply.PuckIFF(sigma_2, sigma_1, sigma_1D, tau_21, tau_21c, S21, RAperpperp, p_perppara_plus, p_perppara_minus, p_perpperp_minus, Y_T, Y_C)
           
                # FPF Condition
                # if ply.failuremode is not None:    
                # Max Stress: (any failure is considered FPF) # INCORRECT
                # if ply.failuremode == "Tensile Fibre" or ply.failuremode == "Compressive" or ply.failuremode == "Tensile Inter-Fibre" or ply.failuremode=="Compressive Inter-Fibre" or ply.failuremode ==  "Shear Parallel to Fibres":
                    
                if failuretracking == 2:
                    firstplyfailureoccurence = True
                    Pf_arr[i,j] = 1/N
            print(f'loop count: {loop_counter}, i = {i}, j = {j}, k = {N}, isFPF?: {firstplyfailureoccurence} ')

Pf = np.mean(Pf_arr)
print(f'\nLoad: {N_load}[N/m] => Probability of Failure: {Pf}')


# Record end time
end_time = time.time()

# Compute elapsed time
elapsed_time = end_time - start_time

print("Time taken:", elapsed_time, "seconds")