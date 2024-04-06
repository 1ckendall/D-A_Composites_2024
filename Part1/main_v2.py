import matplotlib.pyplot as plt
import numpy as np
from pprint import pprint
from Class import Lamina, Laminate

# plot controls (Booleans)
plot_1a = True # plot 1a: engineering constants for laminae (default = True)
plot_1b = False # plot 1b: stress-strain (default = True)

# savefig controls
savefig = False # default = False


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

n_arr = np.array((1, 5, 10,)) #,2,5,10)) # number of repetitions of symmetric laminate
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
        
if plot_1a:      
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    # Plot 1: Inplane Engineering Stiffness
    for i in range(len(n_arr)):
        if i == 0:
            axes[0, 0].plot(theta_arr, Ex_arr[i, :], color='blue', label=r'$E_{x}$')
            axes[0, 0].plot(theta_arr, Ey_arr[i, :], color='red', label=r'$E_{y}$')
            axes[0, 0].plot(theta_arr, Gxy_arr[i, :], color='green', label=r'$G_{xy}$')
            if np.allclose(Ex_arr[1:,:], Ex_arr[:-1,:]):
                text2plot = 'n = ' + ','.join(map(str, n_arr))
                axes[0, 0].text(60, 0.5*np.max(Ex_arr), text2plot, verticalalignment = 'top', color=(0, 0, 1))
            if np.allclose(Ey_arr[1:,:], Ey_arr[:-1,:]):
                text2plot = 'n = ' + ','.join(map(str, n_arr))
                axes[0, 0].text(60, 0.9*np.max(Ey_arr), text2plot,  verticalalignment = 'top', color=(1, 0, 0))
            if np.allclose(Gxy_arr[1:,:], Gxy_arr[:-1,:]):
                text2plot = 'n = ' + ','.join(map(str, n_arr))
                axes[0, 0].text(50, np.max(Gxy_arr), text2plot,  verticalalignment = 'bottom', color=(0, 0.5, 0))
        else:
            axes[0, 0].plot(theta_arr, Ex_arr[i, :], color='blue')
            axes[0, 0].plot(theta_arr, Ey_arr[i, :], color='red')
            axes[0, 0].plot(theta_arr, Gxy_arr[i, :], color='green')  
    #axes[0, 0].grid()
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
                axes[0, 1].text(40, 1.1*np.max(vxy_arr), text2plot, color=(1, 0, 0))
            if np.allclose(vyx_arr[1:,:], vyx_arr[:-1,:]):
                text2plot = 'n = ' + ','.join(map(str, n_arr))
                axes[0, 1].text(30, np.max(vyx_arr), text2plot, color=(0, 0, 1))
        else:
            axes[0, 1].plot(theta_arr, vxy_arr[i, :], color='red')
            axes[0, 1].plot(theta_arr, vyx_arr[i, :], color='blue')
    #axes[0, 1].grid()
    axes[0, 1].legend()
    axes[0, 1].set_title('Inplane Poisson Ratio')
    axes[0, 1].set_xlabel(r'$\theta$ [deg]')
    axes[0, 1].set_ylabel('[-]')
    
    # Plot 3: Flexural Engineering Stiffness
    opacity_control = 0.3
    valign_control = 0.05
    for i in range(len(n_arr)):
        alpha = 1 - opacity_control*i
        if i == 0:
            axes[1, 0].plot(theta_arr, Exb_arr[i, :], alpha = alpha,color='blue', label=r'$E_{xb}$')
            axes[1, 0].plot(theta_arr, Eyb_arr[i, :], alpha = alpha,color='red', label=r'$E_{yb}$')
            axes[1, 0].plot(theta_arr, Gxyb_arr[i, :], alpha = alpha,color='green', label=r'$G_{xyb}$')
        else:
            axes[1, 0].plot(theta_arr, Exb_arr[i, :], alpha = alpha,color='blue')
            axes[1, 0].plot(theta_arr, Eyb_arr[i, :], alpha = alpha, color='red')
            axes[1, 0].plot(theta_arr, Gxyb_arr[i, :],alpha = alpha, color='green')
        text2plot = f'n = {n_arr[i]}'
        if i == 0: # hacky formatting fix
            axes[1, 0].text(theta_arr[np.argmax(Exb_arr[i, :])], 1.02*np.max(Exb_arr[i, :]), text2plot, verticalalignment = 'top', color=(0, 0, 1, alpha))
        elif i == 1: 
            axes[1, 0].text(theta_arr[np.argmax(Exb_arr[i, :])], 1.05*np.max(Exb_arr[i, :]), text2plot,verticalalignment = 'top', color=(0, 0, 1, alpha))
        else: 
            axes[1, 0].text(theta_arr[np.argmax(Exb_arr[i, :])], 0.9*np.max(Exb_arr[i, :]), text2plot,verticalalignment = 'bottom', color=(0, 0, 1, alpha))
        axes[1, 0].text(theta_arr[np.argmax(Eyb_arr[i, :])], (0.9+valign_control*i)*np.max(Eyb_arr[i, :]), text2plot, color=(1, 0, 0, alpha))
        axes[1, 0].text(theta_arr[np.argmax(Gxyb_arr[i, :])], (0.8+2*valign_control*i)*np.max(Gxyb_arr[i, :]), text2plot, color=(0, 0.5, 0, alpha), fontsize = 9)
    #axes[1, 0].grid()
    axes[1, 0].legend()
    axes[1, 0].set_title('Flexural Stiffness')
    axes[1, 0].set_xlabel(r'$\theta$ [deg]')
    axes[1, 0].set_ylabel('[Pa]')
    
    # Plot 4: Flexural Engineering Poisson Ratio
    for i in range(len(n_arr)):
        alpha = 1 - opacity_control*i
        if i == 0:
            axes[1, 1].plot(theta_arr, vxyb_arr[i, :], alpha = alpha, color='red', label=r'$\nu_{xyb}$')
            axes[1, 1].plot(theta_arr, vyxb_arr[i, :], alpha = alpha, color='blue', label=r'$\nu_{yxb}$')
        else:
            axes[1, 1].plot(theta_arr, vxyb_arr[i, :],  alpha = alpha,color='red')
            axes[1, 1].plot(theta_arr, vyxb_arr[i, :],  alpha = alpha,color='blue')
        text2plot = f'n = {n_arr[i]}'
        if i != 2: # hacky formatting fix
            axes[1, 1].text(theta_arr[np.argmax(vxyb_arr[i, :])], np.max(vxyb_arr[i, :]), text2plot, color=(1, 0, 0, alpha), fontsize = 8)
        else:
            axes[1, 1].text(theta_arr[np.argmin(vxyb_arr[i, :])], 0.85*np.min(vxyb_arr[i, :]), text2plot, color=(1, 0, 0, alpha), verticalalignment = 'bottom', fontsize = 8)

        axes[1, 1].text(theta_arr[np.argmax(vyxb_arr[i, :])], np.max(vyxb_arr[i, :]), text2plot, color=(0, 0, 1, alpha), fontsize = 8)
    #axes[1, 1].grid()
    axes[1, 1].legend()
    axes[1, 1].set_title('Flexural Poisson Ratio')
    axes[1, 1].set_xlabel(r'$\theta$ [deg]')
    axes[1, 1].set_ylabel('[-]')
    plt.tight_layout()
    plt.show()
    
    if savefig:
        plt.savefig('1a_Engineering_Constants_v1.png', dpi=500)

# Question 1b: Lamina Stress and Strain
Nx = 0.2E2 # [N/m]
Ny = 1.8E4 # [N/m]
Mx = 18E3 # [N]

layup =  [0, 0, 90, 30, 90]  #verification (abd): [0, 30, -30, 90] ##transverse symmetric: [90, 30, 90, 0, 0, 0, 0, 90, 30, 90] # symmetric: [0, 0, 90, 30, 90, 90, 30, 90, 0, 0]  # truncated: [0, 0, 90, 30, 90] 
plylist = []

for angle in layup:
    plylist.append(Lamina(angle,E1,E2,G12,v12,Xt,Xc,Yt,Yc,S,t))

Laminate1 = Laminate(plylist, Nx=Nx, Ny=Ny, Ns=0, Mx=Mx, My=0, Ms=0)

Laminate1.getStressStrain()


if plot_1b:
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    # fig.subplots_adjust(hspace=0.5)  # Adjust vertical spacing between subplots
    
    y2plot = Laminate1.z
    
    # Plotting Strain along e11
    x2plot = Laminate1.e11
    x2plot = np.append(x2plot, Laminate1.e11[-1])
    ax = axes[0, 0]
    ax.step(x2plot, y2plot, where='pre', color='red', linewidth=2)  
    ax.plot((x2plot[0], 0), (y2plot[0], y2plot[0]), color='red', linewidth=2)  
    ax.plot((x2plot[-1], 0), (y2plot[-1], y2plot[-1]), color='red', linewidth=2)  
    ax.plot((0, 0), (y2plot[0], y2plot[-1]), color='black')  
    ax.plot((0, 0), (y2plot[0], y2plot[-1]), color='black')  # Adjusted for midplane
    ax.plot((np.min(x2plot), np.max(x2plot)), (0, 0), color='black')
    setscale = 0.2
    textxpos = np.min(x2plot) + (np.max(x2plot)-np.min(x2plot))*setscale
    for j in range(len(Laminate1.plys)):
        text2plot = f'k={j + 1} '
        textypos = Laminate1.z_lamina_midplane[j]
        ax.text(textxpos, textypos, text2plot, verticalalignment='bottom', horizontalalignment='right')
    for z in Laminate1.z[1:-1]:
        ax.plot((np.min(x2plot), np.max(x2plot)), (z, z), color='black', linestyle='--', linewidth = 1)
    ax.set_xlabel(r'$\epsilon_{11}$ [-]')
    ax.set_ylabel('z [m]')
    
    
    # Plotting Strain along e22
    x2plot = Laminate1.e22
    x2plot = np.append(x2plot, Laminate1.e22[-1])
    ax = axes[0, 1]
    ax.step(x2plot, y2plot, where='pre', color='red', linewidth=2)  
    ax.plot((x2plot[0], 0), (y2plot[0], y2plot[0]), color='red', linewidth=2)  
    ax.plot((x2plot[-1], 0), (y2plot[-1], y2plot[-1]), color='red', linewidth=2)  
    ax.plot((0, 0), (y2plot[0], y2plot[-1]), color='black')  # Adjusted for midplane
    ax.plot((np.min(x2plot), np.max(x2plot)), (0, 0), color='black')
    setscale = 0.2
    textxpos = np.min(x2plot) + (np.max(x2plot)-np.min(x2plot))*setscale
    for j in range(len(Laminate1.plys)):
        text2plot = f'k={j + 1} '
        textypos = Laminate1.z_lamina_midplane[j]
        ax.text(textxpos, textypos, text2plot, verticalalignment='bottom', horizontalalignment='right')
    for z in Laminate1.z[1:-1]:
        ax.plot((np.min(x2plot), np.max(x2plot)), (z, z), color='black', linestyle='--', linewidth = 1)
    ax.set_xlabel(r'$\epsilon_{22}$ [-]')
    # ax.set_ylabel('z [m]')
    ax.set_yticklabels([])
    
    
    # Plotting Strain along e12
    x2plot = Laminate1.e12
    x2plot = np.append(x2plot, Laminate1.e12[-1])
    ax = axes[0, 2]
    ax.step(x2plot, y2plot, where='pre', color='red', linewidth=2)  
    ax.plot((x2plot[0], 0), (y2plot[0], y2plot[0]), color='red', linewidth=2)  
    ax.plot((x2plot[-1], 0), (y2plot[-1], y2plot[-1]), color='red', linewidth=2)  
    ax.plot((0, 0), (y2plot[0], y2plot[-1]), color='black')  # Adjusted for midplane
    ax.plot((np.min(x2plot), np.max(x2plot)), (0, 0), color='black')
    setscale = 0.2
    textxpos = np.min(x2plot) + (np.max(x2plot)-np.min(x2plot))*setscale
    for j in range(len(Laminate1.plys)):
        text2plot = f'k={j + 1} '
        textypos = Laminate1.z_lamina_midplane[j]
        ax.text(textxpos, textypos, text2plot, verticalalignment='bottom', horizontalalignment='right')
    for z in Laminate1.z[1:-1]:
        ax.plot((np.min(x2plot), np.max(x2plot)), (z, z), color='black', linestyle='--', linewidth = 1)
    ax.set_xlabel(r'$\epsilon_{12}$ []')
    # ax.set_ylabel('z [m]')
    ax.set_yticklabels([])
    
    # Plotting Stress along sigma11
    x2plot = Laminate1.sigma11
    x2plot = np.append(x2plot, Laminate1.sigma11[-1])
    ax = axes[1, 0]
    ax.step(x2plot, y2plot, where='pre', color='blue', linewidth=2)  
    ax.plot((x2plot[0], 0), (y2plot[0], y2plot[0]), color='blue', linewidth=2)  
    ax.plot((x2plot[-1], 0), (y2plot[-1], y2plot[-1]), color='blue', linewidth=2)  
    ax.plot((0, 0), (y2plot[0], y2plot[-1]), color='black')  # Adjusted for midplane
    ax.plot((np.min(x2plot), np.max(x2plot)), (0, 0), color='black')
    setscale = 0.2
    textxpos = np.min(x2plot) + (np.max(x2plot)-np.min(x2plot))*setscale
    for j in range(len(Laminate1.plys)):
        text2plot = f'k={j + 1} '
        textypos = Laminate1.z_lamina_midplane[j]
        ax.text(textxpos, textypos, text2plot, verticalalignment='bottom', horizontalalignment='right')
    for z in Laminate1.z[1:-1]:
        ax.plot((np.min(x2plot), np.max(x2plot)), (z, z), color='black', linestyle='--', linewidth = 1)
    ax.set_xlabel(r'$\sigma_{11}$ [Pa]')
    ax.set_ylabel('z [m]')
    
    
    # Plotting Stress along sigma22
    x2plot = Laminate1.sigma22
    x2plot = np.append(x2plot, Laminate1.sigma22[-1])
    ax = axes[1, 1]
    ax.step(x2plot, y2plot, where='pre', color='blue', linewidth=2)  
    ax.plot((x2plot[0], 0), (y2plot[0], y2plot[0]), color='blue', linewidth=2)  
    ax.plot((x2plot[-1], 0), (y2plot[-1], y2plot[-1]), color='blue', linewidth=2)  
    ax.plot((0, 0), (y2plot[0], y2plot[-1]), color='black')  # Adjusted for midplane
    ax.plot((np.min(x2plot), np.max(x2plot)), (0, 0), color='black')
    setscale = 0.2
    textxpos = np.min(x2plot) + (np.max(x2plot)-np.min(x2plot))*setscale
    for j in range(len(Laminate1.plys)):
        text2plot = f'k={j + 1} '
        textypos = Laminate1.z_lamina_midplane[j]
        ax.text(textxpos, textypos, text2plot, verticalalignment='bottom', horizontalalignment='right')
    for z in Laminate1.z[1:-1]:
        ax.plot((np.min(x2plot), np.max(x2plot)), (z, z), color='black', linestyle='--', linewidth = 1)
    ax.set_xlabel(r'$\sigma_{22}$ [Pa]')
    # ax.set_ylabel('z [m]')
    ax.set_yticklabels([])
    
    # Plotting Stress along sigma12
    x2plot = Laminate1.sigma12
    x2plot = np.append(x2plot, Laminate1.sigma12[-1])
    ax = axes[1, 2]
    ax.step(x2plot, y2plot, where='pre', color='blue', linewidth=2)  
    ax.plot((x2plot[0], 0), (y2plot[0], y2plot[0]), color='blue', linewidth=2)  
    ax.plot((x2plot[-1], 0), (y2plot[-1], y2plot[-1]), color='blue', linewidth=2)  
    ax.plot((0, 0), (y2plot[0], y2plot[-1]), color='black')  # Adjusted for midplane
    ax.plot((np.min(x2plot), np.max(x2plot)), (0, 0), color='black')
    setscale = 0.2
    textxpos = np.min(x2plot) + (np.max(x2plot)-np.min(x2plot))*setscale
    for j in range(len(Laminate1.plys)):
        text2plot = f'k={j + 1} '
        textypos = Laminate1.z_lamina_midplane[j]
        ax.text(textxpos, textypos, text2plot, verticalalignment='bottom', horizontalalignment='right')
    for z in Laminate1.z[1:-1]:
        ax.plot((np.min(x2plot), np.max(x2plot)), (z, z), color='black', linestyle='--', linewidth = 1)
    ax.set_xlabel(r'$\sigma_{12}$ [Pa]')
    # ax.set_ylabel('z [m]')
    ax.set_yticklabels([])
    
    plt.tight_layout()
    plt.show()
    

    if savefig:
        plt.savefig('1b_stress_strain_v2.png', dpi=500)


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
        
anglelist = [0,90,45,-45,-45,45,90,0,0,90,45,-45,-45,45,90,0]
stressinputvector = np.arange(1E3,1E8,500)
# angleinputvector = np.arange(45,46,1)
angleinputvector = np.arange(0,1,1)



firstfailuremaxstress = np.zeros((len(angleinputvector),2))
firstfaileglobalstrain =np.zeros((len(angleinputvector),2))
firstfailurePUCK =np.empty((0,2),dtype=int)
lastplyfailuremaxstress=np.zeros((len(angleinputvector),2))
lastplygolbalstrain = np.zeros((len(angleinputvector),2))
lastplypuck=np.empty((0,2),dtype=int)

for i in range(len(angleinputvector)):
  #initialiazing fresh lamina properties for new angle
  
 
  E1 = np.full(len(anglelist),E1)
  E2 = np.full(len(anglelist),E2)
  v12 = np.full(len(anglelist),v12)
  G12 = np.full(len(anglelist),G12)
  firstplyfailureoccurence = False
  lastplyfailureoccurence = False
  
  failuretracking = np.zeros(len(anglelist))
  failedplys = [False] * len(anglelist)

  
  print(f'Load Angle [deg]: {angleinputvector[i]}')
  for j in stressinputvector: 
     if lastplyfailureoccurence==False:
            plylist = []
            
            stressloading = np.array([0,j,j])
            

            stressused = stressloading
            
            nx=stressused[0]
            ny=stressused[1]
            ns=stressused[2]
            for k in range(len(anglelist)): 
                plylist.append(Lamina(anglelist[k],E1[k],E2[k],G12[k],v12[k],Xt,Xc,Yt,Yc,S,t))
            LaminateQ3 = Laminate(plylist,Nx=nx,Ny=ny,Ns=ns,Mx=0,My=0,Ms=0,Loadangle=angleinputvector[i])
            LaminateQ3.getstressstrainEnvelope()




            
            
            
            for k in range(len( anglelist)): 
                    laminaangle = anglelist[k]
                    Sigma1 = LaminateQ3.sigma11[k]
                    Sigma2 = LaminateQ3.sigma22[k]
                    tau21 =  LaminateQ3.sigma12[k]
                    failurechecking1  = Lamina(laminaangle,E1[k],E2[k],G12[k],v12[k],Xt,Xc,Yt,Yc,S,t,sigma_1=Sigma1,sigma_2=Sigma2,tau_21=tau21)
                    failurechecking2  = Lamina(laminaangle,E1[k],E2[k],G12[k],v12[k],Xt,Xc,Yt,Yc,S,t,sigma_1=Sigma1,sigma_2=Sigma2,tau_21=tau21)
                    failurechecking3  = Lamina(laminaangle,E1[k],E2[k],G12[k],v12[k],Xt,Xc,Yt,Yc,S,t,sigma_1=Sigma1,sigma_2=Sigma2,tau_21=tau21)
                    failurechecking1.maxStressFibreFail() #fiber failure
                    failurechecking2.maxStressInterFibreFail() #InterFiberfailure
                    failurechecking3.maxStressShearFail() #Shearfailure
                    if failurechecking1.failuremode == 'Tensile Fibre' or failurechecking1.failuremode == 'Compressive Fibre' :
                        failuretracking[k] = 2
                        E1[k] = 1E-25
                        E2[k]=1E-25
                        v12[k] =1E-25
                        G12[k]=1E-25
                        failedplys[k] =True
                    if failurechecking2.failuremode == 'Tensile Inter-Fibre' or failurechecking2.failuremode=="Compressive Inter-Fibre":
                            failuretracking[k] +=1
                            E2[k] = 0.1*E2[k] 
                            v12[k] =0.1*v12[k]
                            G12[k]=0.1*G12[k]
                    if failurechecking3.failuremode ==  "Shear Parallel to Fibres":
                            failuretracking[k] +=1
                            E2[k] = 0.1*E2[k] 
                            v12[k] =0.1*v12[k]
                            G12[k]=0.1*G12[k]
                    for firsplyfail in  failuretracking : 
                        if firsplyfail >= 2 and  not firstplyfailureoccurence: 
                            firstplyfailureoccurence = True
                            firstfailureloadY =LaminateQ3.sigmayprime
                            # print(failurechecking1.failuremode)
                            # print('firstply',k)
                            print(f'First-Ply Failure [FPF]: Ply index: {k}, Fibre angle: {laminaangle}')
                            firstfailureloadS = LaminateQ3.sigmaxyprime
                            
                            firstfailuremaxstress[i,0] = firstfailureloadY
                            firstfailuremaxstress[i,1] = firstfailureloadS
                            firstfaileglobalstrain[i,0] = LaminateQ3.eyglobal
                            firstfaileglobalstrain[i,1] = LaminateQ3.esglobal
                            print('loads',firstfailureloadY,firstfailureloadS)            
                    if 2<=  failuretracking[k] :
                        
                        E1[k] = 1E-25
                        E2[k]=1E-25
                        v12[k] =1E-25
                        G12[k]=1E-25
                    if lastplyfailureoccurence == False: 
                        for entry in failuretracking:
                            if all(entry >=2 for entry in failuretracking):
                                lastplyfailureoccurence = True 
                                print(failurechecking1.failuremode)
                                print('lastply',k)
                                lastfailureloadY =LaminateQ3.sigmayprime
                                lastfailureloadS = LaminateQ3.sigmaxyprime
                                print('lastY',lastfailureloadY)
                                print('lastS',lastfailureloadS)
                                lastplyfailuremaxstress[i,0] = lastfailureloadY
                                lastplyfailuremaxstress[i,1] =lastfailureloadS 
                                lastplygolbalstrain[i,0] = LaminateQ3.eyglobal
                                lastplygolbalstrain[i,1] = LaminateQ3.esglobal
                                print(failuretracking)
                    
        
                 
print(failuretracking)
        
#making the scatter plot maxstress and puck 
fig,(ax1,ax2) =plt.subplots(1,2)
#failure envelope for max stress criteria 

ax1.scatter(firstfailuremaxstress[:,0]/1E6,firstfailuremaxstress[:,1]/1E6,marker = '^',color='red') #first ply failure 
ax1.scatter(lastplyfailuremaxstress[:,0]/1E6,lastplyfailuremaxstress[:,1]/1E6,marker='o',color='blue') #last ply failure 
ax1.set_xlabel('sigma_y')
ax1.set_ylabel('sigma_s')
ax1.set_title('Max stress Criteria failure envelope')
plt.show()

ax2.scatter(firstfaileglobalstrain[:,0],firstfaileglobalstrain[:,1],marker = '^',color='red') #first ply failure puck
ax2.scatter(lastplygolbalstrain[:,0],lastplygolbalstrain[:,1],marker='o',color='blue') #last ply failure puck
ax2.set_xlabel('e_y')
ax2.set_ylabel('e_s') 
ax2.set_title('Puck Criteria failure envelope')

plt.tight_layout()


      