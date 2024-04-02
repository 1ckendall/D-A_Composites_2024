import matplotlib.pyplot as plt
import numpy as np
from pprint import pprint
from Class import Lamina, Laminate

# UD Lamina Mean Properties 
E1 = 145.3E9 # [Pa]
E2 = 8.5E9 # [Pa]
v12 = 0.31 # [--]
G12 = 4.58E9 # [Pa]
Xt = 1932E6 # [Pa]
Yt = 108E6 # [Pa]
Xc = 1480E6 # [Pa]
Yc = 220E6 # [Pa]
S = 132.8E6 # [Pa]
t= 0.125E-3 # [m] # free variable, can be changed 

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

Exb_arr = np.zeros((len(n_arr), len(theta_arr)))
Eyb_arr = np.zeros((len(n_arr), len(theta_arr)))
Gxyb_arr = np.zeros((len(n_arr), len(theta_arr)))
vxyb_arr = np.zeros((len(n_arr), len(theta_arr)))
vyxb_arr = np.zeros((len(n_arr), len(theta_arr)))

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
        
        Exb_arr[i, j] = Laminate_ij.Exb
        Eyb_arr[i, j] = Laminate_ij.Eyb
        Gxyb_arr[i, j] = Laminate_ij.Gxyb
        vxyb_arr[i, j] = Laminate_ij.vxyb
        vyxb_arr[i, j] = Laminate_ij.vyxb

        
        
fig, axes = plt.subplots(2, 2, figsize=(10, 8))
# Plot 1: Inplane Engineering Stiffness
for i in range(len(n_arr)):
    if i == 0:
        axes[0, 0].plot(theta_arr, Ex_arr[i, :], color='blue', label=r'$E_{x}$')
        axes[0, 0].plot(theta_arr, Ey_arr[i, :], color='red', label=r'$E_{y}$')
        axes[0, 0].plot(theta_arr, Gxy_arr[i, :], color='green', label=r'$G_{xy}$')
        if np.allclose(Ex_arr[1:,:], Ex_arr[:-1,:]):
            text2plot = 'n = ' + ','.join(map(str, n_arr))
            axes[0, 0].text(60, 0.4*np.max(Ex_arr), text2plot)
        if np.allclose(Ey_arr[1:,:], Ey_arr[:-1,:]):
            text2plot = 'n = ' + ','.join(map(str, n_arr))
            axes[0, 0].text(60, 0.9*np.max(Ey_arr), text2plot)
        if np.allclose(Gxy_arr[1:,:], Gxy_arr[:-1,:]):
            text2plot = 'n = ' + ','.join(map(str, n_arr))
            axes[0, 0].text(50, np.max(Gxy_arr), text2plot)
    else:
        axes[0, 0].plot(theta_arr, Ex_arr[i, :], color='blue')
        axes[0, 0].plot(theta_arr, Ey_arr[i, :], color='red')
        axes[0, 0].plot(theta_arr, Gxy_arr[i, :], color='green')  
axes[0, 0].grid()
axes[0, 0].legend()
axes[0, 0].set_title('Inplane Stiffness')
axes[0, 0].set_xlabel(r'$\theta$ [deg]')
axes[0, 0].set_ylabel('[Pa]')

# Plot 2: Inplane Engineering Poisson Ratio
for i in range(len(n_arr)):
    if i == 0:
        axes[0, 1].plot(theta_arr, vxy_arr[i, :], color='red', label=r'$\nu_{xy}$')
        axes[0, 1].plot(theta_arr, vyx_arr[i, :], color='blue', label=r'$\nu_{yx}$')
        if np.allclose(vxy_arr[1:,:], vxy_arr[:-1,:]):
            text2plot = 'n = ' + ','.join(map(str, n_arr))
            axes[0, 1].text(40, np.max(vxy_arr), text2plot)
        if np.allclose(vyx_arr[1:,:], vyx_arr[:-1,:]):
            text2plot = 'n = ' + ','.join(map(str, n_arr))
            axes[0, 1].text(30, 0.9*np.max(vyx_arr), text2plot)
    else:
        axes[0, 1].plot(theta_arr, vxy_arr[i, :], color='red')
        axes[0, 1].plot(theta_arr, vyx_arr[i, :], color='blue')
axes[0, 1].grid()
axes[0, 1].legend()
axes[0, 1].set_title('Inplane Poisson Ratio')
axes[0, 1].set_xlabel(r'$\theta [deg]$')
axes[0, 1].set_ylabel('[-]')

# Plot 3: Flexural Engineering Stiffness
for i in range(len(n_arr)):
    if i == 0:
        axes[1, 0].plot(theta_arr, Exb_arr[i, :], color='blue', label=r'$E_{xb}$')
        axes[1, 0].plot(theta_arr, Eyb_arr[i, :], color='red', label=r'$E_{yb}$')
        axes[1, 0].plot(theta_arr, Gxyb_arr[i, :], color='green', label=r'$G_{xyb}$')
    else:
        axes[1, 0].plot(theta_arr, Exb_arr[i, :], color='blue')
        axes[1, 0].plot(theta_arr, Eyb_arr[i, :], color='red')
        axes[1, 0].plot(theta_arr, Gxyb_arr[i, :], color='green')
    text2plot = f'n = {n_arr[i]}'
    axes[1, 0].text(theta_arr[np.argmax(Exb_arr[i, :])], np.max(Exb_arr[i, :]), text2plot)
    axes[1, 0].text(theta_arr[np.argmax(Eyb_arr[i, :])], np.max(Eyb_arr[i, :]), text2plot)
    axes[1, 0].text(theta_arr[np.argmax(Gxyb_arr[i, :])], np.max(Gxyb_arr[i, :]), text2plot)
axes[1, 0].grid()
axes[1, 0].legend()
axes[1, 0].set_title('Flexural Stiffness')
axes[1, 0].set_xlabel(r'$\theta [deg]$')
axes[1, 0].set_ylabel('[Pa]')

# Plot 4: Flexural Engineering Poisson Ratio
for i in range(len(n_arr)):
    if i == 0:
        axes[1, 1].plot(theta_arr, vxyb_arr[i, :], color='red', label=r'$\nu_{xyb}$')
        axes[1, 1].plot(theta_arr, vyxb_arr[i, :], color='blue', label=r'$\nu_{yxb}$')
    else:
        axes[1, 1].plot(theta_arr, vxyb_arr[i, :], color='red')
        axes[1, 1].plot(theta_arr, vyxb_arr[i, :], color='blue')
    text2plot = f'n = {n_arr[i]}'
    axes[1, 1].text(theta_arr[np.argmax(vxyb_arr[i, :])], np.max(vxyb_arr[i, :]), text2plot)
    axes[1, 1].text(theta_arr[np.argmax(vyxb_arr[i, :])], np.max(vyxb_arr[i, :]), text2plot)
axes[1, 1].grid()
axes[1, 1].legend()
axes[1, 1].set_title('Flexural Poisson Ratio')
axes[1, 1].set_xlabel(r'$\theta [deg]$')
axes[1, 1].set_ylabel('[-]')

plt.tight_layout()
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
      