import matplotlib.pyplot as plt
import numpy as np
from pprint import pprint
from Class import Lamina, Laminate

# UD Lamina Mean Properties 
E1 = 145.3E3 # [MPa]
E2 = 8.5E3 # [MPa]
v12 = 0.31 # [--]
G12 = 4.58E3 # [MPa]
Xt = 1932 # [MPa]
Yt = 108 # [MPa]
Xc = 1480 # [MPa]
Yc = 220 # [MPa]
S = 132.8 # [MPa]
t= 0.125 # [mm] # free variable, can be changed 

# UD Lamina Std. Dev Properties 
pass

# Question 1: ABD Matrices and Stress Analysis
# Question 1a: Engineering Constants
# Consider symmtric laminate [15, +theta, - theta, 75, 75]_ns

n_arr = np.array((10,)) #,2,5,10)) # number of repetitions of symmetric laminate
theta_arr = np.arange(0,90,5) # theta variation in laminate

# initialise engineering constants
Ex_arr = np.zeros((len(n_arr), len(theta_arr)))
Ey_arr = np.zeros((len(n_arr), len(theta_arr)))
Gxy_arr = np.zeros((len(n_arr), len(theta_arr)))
vxy_arr = np.zeros((len(n_arr), len(theta_arr)))
vyx_arr = np.zeros((len(n_arr), len(theta_arr)))

# simulation loop
for i in range(len(n_arr)):
    for j in range(len(theta_arr)):
        theta = theta_arr[j]
        # layup = [theta]
        layup = [15, theta, -theta, 75, 75, 75, 75, -theta, theta, 15]
        layup *= (i+1)
        
        plylist= [] 
        
        for angle in layup:
            plylist.append(Lamina(angle,E1,E2,G12,v12,Xt,Xc,Yt,Yc,S,t))
        
        Laminate_ij = Laminate(plylist)
        
        # Obtain engineering constants
        Ex_arr[i, j] = Laminate_ij.Ex
        Ey_arr[i, j] = Laminate_ij.Ey
        Gxy_arr[i, j] = Laminate_ij.Gxy
        vxy_arr[i, j] = Laminate_ij.vxy
        vyx_arr[i, j] = Laminate_ij.vyx
        
        
# Create a figure and subplots
plt.figure('1')
plt.clf()
for i in range(len(n_arr)):
    plt.plot(theta_arr, Ex_arr[i, :] , color='blue', label = 'Ex')
    plt.plot(theta_arr, Ey_arr[i, :] , color='red', label = 'Ey')
    plt.plot(theta_arr, Gxy_arr[i, :] , color='green', label = 'Gxy')
plt.grid()
plt.legend()
plt.title('Ex, Ey, Gxy')
plt.xlabel(r'$\theta$')
plt.ylabel('[MPa]')

plt.figure('2')
plt.clf()
for i in range(len(n_arr)):
    plt.plot(theta_arr, vxy_arr[i, :] ,color='red', label = 'vxy')
    plt.plot(theta_arr, vyx_arr[i, :] ,color='blue', label = 'vyx')
plt.grid()
plt.legend()
plt.title('vxy, vyx')
plt.xlabel(r'$\theta$')
plt.ylabel('[-]')

plt.show()

# if __name__ == "__main__":
#     print("Running Main")
#     E1 = 140
#     E2 = 10
#     G12 = 5
#     v12 = 0.3
#     t = 0.125
#     Xt = 1500
#     Xc = 1200
#     Yt = 50
#     Yc = 250
#     S = 70
#     anglelist = [0, 45, -45, 90, 90, -45, 45, 0]
#     plylist = []
#     for angle in anglelist:
#         plylist.append(Lamina(angle, E1, E2, G12, v12, Xt, Xc, Yt, Yc, S, t))
#     L = Laminate(plylist)
#     pprint(L.ABD)

#Question  2 damage progression 
anglelist = [0,90,45,-45,-45,45,90,0,0,90,45,-45,-45,45,90,0]
stressinputvector = np.linspace(0,10000,10) 
angleinputvector = np.radians(np.linspace(0,360,10))
plylist = []
failuretracking = np.full(len(anglelist), True, dtype=bool)
firstfailuremaxstress = []
firstfailurePUCK =[]
lastplyfailuremaxstress= []
lastplypuck=[]

for i in angleinputvector: 
  for j in stressinputvector: 
      m = np.cos(i)
      n = np.sin(i)
      stressloading = np.array([0,j,0])
      print(i)
      print(j)
      stresstransformmatrix =np.array([[m**2, n**2,-2*m*n],
                              [n**2,m**2,2*m*n],
                                [-m*n,m*n,m**2-n**2]])
      stressused = stresstransformmatrix @ stressloading
      print(stressused)
      