import matplotlib.pyplot as plt
import numpy as np
from pprint import pprint
from Class import Lamina, Laminate
import time

start_time = time.time()


# plot controls (Booleans)
plot_2a= True # plot 2a: failure envelopes (default = True)
plot_2b= True # plot 2b: damage tolerance (default = True)


# run controls (Booleans)
run_maxStress = True # run failure envelope for max stress (default = True)
run_Puck = True # run failure envelope for Puck (default = True)


# savefig controls (Booleans)
savefig = False # default = False


# UD Lamina Mean Properties 
E1_float = 145.3E9 # [Pa]
E2_float = 8.5E9 # [Pa]
v12_float = 0.31 # [--]
G12_float = 4.58E9 # [Pa]
Xt = 1932E6 # [Pa]
Yt = 108E6 # [Pa]
Xc = 1480E6 # [Pa]
Yc = 220E6 # [Pa]
S = 132.8E6 # [Pa]
t= 0.125E-3 # [m] # free variable, can be changed 


# UD Lamina Std. Dev Properties 
pass

#Question 2: Progressive Damage Analysis 

## 2a: Failure Envelope
anglelist = [0,90,45,-45,-45,45,90,0,0,90,45,-45,-45,45,90,0]

stressinputvector = np.arange(1E5,1E7,10000) ## UNCOMMENT FOR SHORT RUN
angleinputvector = np.arange(0,450,90)  ## UNCOMMENT FOR SHORT RUN

# stressinputvector = np.arange(1E3,1E8,500) ## UNCOMMENT FOR FINAL RUN
# angleinputvector = np.arange(0,390,30) ## UNCOMMENT FOR FINAL RUN

## 2ai: Failure Theory of Maximum Stress 
firstfailure_MaxStress = np.zeros((len(angleinputvector),3))
firstfailure_globalstrain_MaxStress =np.zeros((len(angleinputvector),3))

lastfailure_MaxStress=np.zeros((len(angleinputvector),3))
lastfailure_globalstrain_MaxStress = np.zeros((len(angleinputvector),3))

if run_maxStress:
    print('Computing Failure Envelope based on Max. Stress criterion...')
    for i in range(len(angleinputvector)):
      #initialiazing fresh lamina properties for new angle
      E1 = np.full(len(anglelist),E1_float)
      E2 = np.full(len(anglelist),E2_float)
      v12 = np.full(len(anglelist),v12_float)
      G12 = np.full(len(anglelist),G12_float)
      firstplyfailureoccurence = False
      lastplyfailureoccurence = False
      
      failuretracking = np.zeros(len(anglelist))
      failedplys = [False] * len(anglelist)
    
      
      print(f'Load Angle [deg]: {angleinputvector[i]}\n')
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
                # print(f'Midplane strain (global): ey:{LaminateQ3.eyglobal}, es:{LaminateQ3.esglobal}')
                
                for k in range(len( anglelist)): 
                        laminaangle = anglelist[k]
                        Sigma1 = LaminateQ3.sigma11[k]
                        Sigma2 = LaminateQ3.sigma22[k]
                        tau21 =  LaminateQ3.sigma12[k]
                        failurechecking  = Lamina(laminaangle,E1[k],E2[k],G12[k],v12[k],Xt,Xc,Yt,Yc,S,t,sigma_1=Sigma1,sigma_2=Sigma2,tau_21=tau21)
                        failurechecking.maxStressFibreFail() #fiber failure
                        failurechecking.maxStressInterFibreFail() #InterFiberfailure
                        failurechecking.maxStressShearFail() #Shearfailure
                        if failurechecking.failuremode == 'Tensile Fibre' or failurechecking.failuremode == 'Compressive Fibre' :
                            failuretracking[k] = 2
                            E1[k] = 1E-25
                            E2[k]=1E-25
                            v12[k] =1E-25
                            G12[k]=1E-25
                            failedplys[k] =True
                        if failurechecking.failuremode == 'Tensile Inter-Fibre' or failurechecking.failuremode=="Compressive Inter-Fibre":
                                failuretracking[k] +=1
                                E2[k] = 0.1*E2[k] 
                                v12[k] =0.1*v12[k]
                                G12[k]=0.1*G12[k]
                        if failurechecking.failuremode ==  "Shear Parallel to Fibres":
                                failuretracking[k] +=1
                                E2[k] = 0.1*E2[k] 
                                v12[k] =0.1*v12[k]
                                G12[k]=0.1*G12[k]
                                
                        # first ply failure (FPF): store failure mode, load, stress, strain ply
                        for firsplyfail in  failuretracking : 
                            if firsplyfail >= 2 and  not firstplyfailureoccurence: 
                                firstplyfailureoccurence = True
                                
                                firstfailure_loadX =LaminateQ3.sigmaxprime
                                firstfailure_loadY =LaminateQ3.sigmayprime
                                firstfailure_loadS = LaminateQ3.sigmaxyprime
                                firstfailure_eX = LaminateQ3.exglobal
                                firstfailure_eY = LaminateQ3.eyglobal
                                firstfailure_eS = LaminateQ3.esglobal
                                
                                firstfailure_MaxStress[i,0] = firstfailure_loadX
                                firstfailure_MaxStress[i,1] = firstfailure_loadY
                                firstfailure_MaxStress[i,2] = firstfailure_loadS

                                firstfailure_globalstrain_MaxStress[i,0] = firstfailure_eX
                                firstfailure_globalstrain_MaxStress[i,1] = firstfailure_eY
                                firstfailure_globalstrain_MaxStress[i,2] = firstfailure_eS
                                print(f'First-Ply Failure [FPF]: Ply index: {k}, Fibre angle: {laminaangle}, \nFailure Mode: {failurechecking.failuremode}, \nLoadY: {firstfailure_loadY}, LoadS: {firstfailure_loadS}\ney: {firstfailure_eY}, es: {firstfailure_eS}\n')
                        
                        # ply removal: zero the properties of that lamina
                        if failuretracking[k] >=2:
                            E1[k] = 1E-25
                            E2[k]=1E-25
                            v12[k] =1E-25
                            G12[k]=1E-25                        
                            for entry in failuretracking:
                                if lastplyfailureoccurence == False: 
                                    if all(entry >=2 for entry in failuretracking):
                                        lastplyfailureoccurence = True 
                                        
                                        lastfailure_loadX =LaminateQ3.sigmaxprime
                                        lastfailure_loadY =LaminateQ3.sigmayprime
                                        lastfailure_loadS = LaminateQ3.sigmaxyprime
                                        lastfailure_eX = LaminateQ3.exglobal
                                        lastfailure_eY = LaminateQ3.eyglobal
                                        lastfailure_eS = LaminateQ3.esglobal                           
                                        
                                        lastfailure_MaxStress[i,0] = lastfailure_loadX
                                        lastfailure_MaxStress[i,1] = lastfailure_loadY
                                        lastfailure_MaxStress[i,2] =lastfailure_loadS 
                                        lastfailure_globalstrain_MaxStress[i,0] = lastfailure_eX
                                        lastfailure_globalstrain_MaxStress[i,1] = lastfailure_eY
                                        lastfailure_globalstrain_MaxStress[i,2] = lastfailure_eS
                                        print(f'Last-Ply Failure [LPF]: Ply index: {k}, Fibre angle: {laminaangle}, \nFailure Mode: {failurechecking.failuremode}, LoadY: {lastfailure_loadY}, LoadS: {lastfailure_loadS}\ney: {lastfailure_eY}, es: {lastfailure_eS}\n\n')
                                    
                                        assert np.allclose(failuretracking, 2)
                                        break
            

## 2ai: Failure Theory of Puck
firstfailure_globalstrain_Puck =np.zeros((len(angleinputvector),3))
firstfailure_Puck = np.zeros((len(angleinputvector),3))

lastfailure_globalstrain_Puck = np.zeros((len(angleinputvector),3))
lastfailure_Puck = np.zeros((len(angleinputvector),3))
if run_Puck:
    print('Computing Failure Envelope based on Puck criterion...')
    for i in range(len(angleinputvector)):
      #initialiazing fresh lamina properties for new angle
      E1 = np.full(len(anglelist),E1_float)
      E2 = np.full(len(anglelist),E2_float)
      v12 = np.full(len(anglelist),v12_float)
      G12 = np.full(len(anglelist),G12_float)
      firstplyfailureoccurence = False
      lastplyfailureoccurence = False
      
      failuretracking = np.zeros(len(anglelist))
      failedplys = [False] * len(anglelist)
    
      
      print(f'Load Angle [deg]: {angleinputvector[i]}\n')
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
                        
                        epsilon1 = LaminateQ3.e11[k]
                        epsilon2 = LaminateQ3.e22[k]
                        Gamma21 = LaminateQ3.e12[k]
                        
                        failurechecking  = Lamina(laminaangle,E1[k],E2[k],G12[k],v12[k],Xt,Xc,Yt,Yc,S,t,sigma_1=Sigma1,sigma_2=Sigma2,tau_21=tau21)
                        
                        # Puck 
                        failurechecking.PuckFibreFail(sigma_1 = Sigma1, sigma_2 = Sigma2, sigma_3 = tau21, R_para_t = Xt, R_para_c = Xc, v_perppara = v12[k], E_para = E1[k])  #Puck fiber failure #### currently using Minor poisson, to make it consistent with Minor poisson at failure
                        failurechecking.PuckIFF(sigma_22 = Sigma2, sigma_21 = tau21, sigma_22T_u = Yt, sigma_22C_u=Yc, sigma_12_u=S)  #Puck fiber failure #### currently using Minor poisson, to make it consistent with Minor poisson at failure          
    
                        if failurechecking.failuremode == "FFT" or failurechecking.failuremode == "FFC" :
                            failuretracking[k] = 2
                            E1[k] = 1E-25
                            E2[k]=1E-25
                            v12[k] =1E-25
                            G12[k]=1E-25
                            failedplys[k] =True
                        if failurechecking.failuremode == "IFF A" or failurechecking.failuremode=="IFF B" or failurechecking.failuremode=="IFF C":
                                failuretracking[k] +=1
                                E2[k] = 0.1*E2[k] 
                                v12[k] =0.1*v12[k]
                                G12[k]= 0.1*G12[k]
                                
                        # first ply failure (FPF): store failure mode, load, stress, strain ply
                        for firsplyfail in  failuretracking : 
                            if firsplyfail >= 2 and  not firstplyfailureoccurence: 
                                firstplyfailureoccurence = True
                                
                                firstfailure_loadX = LaminateQ3.sigmaxprime
                                firstfailure_loadY =LaminateQ3.sigmayprime
                                firstfailure_loadS = LaminateQ3.sigmaxyprime
                                firstfailure_eX = LaminateQ3.exglobal
                                firstfailure_eY = LaminateQ3.eyglobal
                                firstfailure_eS = LaminateQ3.esglobal

                                firstfailure_Puck[i,0] = firstfailure_loadX
                                firstfailure_Puck[i,1] = firstfailure_loadY
                                firstfailure_Puck[i,2] = firstfailure_loadS
                                firstfailure_globalstrain_Puck[i,0] = firstfailure_eX
                                firstfailure_globalstrain_Puck[i,1] = firstfailure_eY
                                firstfailure_globalstrain_Puck[i,2] = firstfailure_eS
                                print(f'First-Ply Failure [FPF]: Ply index: {k}, Fibre angle: {laminaangle}, \nFailure Mode: {failurechecking.failuremode}, \nLoadY: {firstfailure_loadY}, LoadS: {firstfailure_loadS}\ney: {firstfailure_eY}, es: {firstfailure_eS}\n')
                        
                        # ply removal: zero the properties of that lamina
                        if failuretracking[k] >=2:
                            E1[k] = 1E-25
                            E2[k]=1E-25
                            v12[k] =1E-25
                            G12[k]=1E-25                        
                            for entry in failuretracking:
                                if lastplyfailureoccurence == False: 
                                    if all(entry >=2 for entry in failuretracking):
                                        lastplyfailureoccurence = True 
                                        
                                        lastfailure_loadX =LaminateQ3.sigmaxprime
                                        lastfailure_loadY =LaminateQ3.sigmayprime
                                        lastfailure_loadS = LaminateQ3.sigmaxyprime
                                        lastfailure_eX = LaminateQ3.exglobal
                                        lastfailure_eY = LaminateQ3.eyglobal
                                        lastfailure_eS = LaminateQ3.esglobal
                                        
                                        lastfailure_Puck[i,0] = lastfailure_loadX
                                        lastfailure_Puck[i,1] = lastfailure_loadY
                                        lastfailure_Puck[i,2] =lastfailure_loadS 
                                        lastfailure_globalstrain_Puck[i,0] = lastfailure_eX
                                        lastfailure_globalstrain_Puck[i,1] = lastfailure_eY
                                        lastfailure_globalstrain_Puck[i,2] = lastfailure_eS
                                        print(f'Last-Ply Failure [LPF]: Ply index: {k}, Fibre angle: {laminaangle}, \nFailure Mode: {failurechecking.failuremode}, LoadY: {lastfailure_loadY}, LoadS: {lastfailure_loadS}\ney: {lastfailure_eY}, es: {lastfailure_eS}\n\n')
                                        
                                        # print(f'Failure index (all plies): {failuretracking}')
                                        assert np.allclose(failuretracking, 2)
                                        break
                                    
# Record end time
end_time = time.time()

# Compute elapsed time
elapsed_time = end_time - start_time

print("Time taken:", elapsed_time, "seconds")
    
if plot_2a:
    # making failure envelope plot (maxstress) sigmay-sigmas  
    plt.figure(1)
    plt.clf()
    plt.plot(firstfailure_MaxStress[:,1]/1E6,firstfailure_MaxStress[:,2]/1E6, linewidth = 1,  color='red', marker = 'x', markersize = 5, label = 'Max Stress: FPF')
    plt.plot(lastfailure_MaxStress[:,1]/1E6,lastfailure_MaxStress[:,2]/1E6, linewidth = 1, color='blue', marker = 'x', markersize = 5, label ='Max Stress: LPF')
    plt.plot(firstfailure_Puck[:,1]/1E6,firstfailure_Puck[:,2]/1E6, linewidth = 1,  linestyle = '--', color='green', label = 'Puck: FPF')
    plt.plot(lastfailure_Puck[:,1]/1E6,lastfailure_Puck[:,2]/1E6, linewidth = 1, linestyle = '--', color='orange', label ='Puck: LPF')
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    plt.legend(fontsize = 16)
    plt.grid(True,alpha=0.5)
    plt.xlabel(r"$\sigma_y$' [MPa]", fontsize = 16)
    plt.ylabel(r"$\sigma_s$' [MPa]", fontsize = 16)
            
    if savefig:
        loadincrement = stressinputvector[1]-stressinputvector[0]
        angleincrement = angleinputvector[1]-angleinputvector[0]
        filename = f'2a_Failure_Envelope_Stress_loadincrement={loadincrement}_angleincrement={angleincrement}.png'
        plt.savefig(filename, dpi=500)
        
    # plt.figure(2)
    # plt.clf()
    # plt.plot(angleinputvector,firstfailure_MaxStress[:,0]/1E6, linewidth = 1,  color='red', marker = 'x', markersize = 5, label ='Max Stress')
    # plt.plot(angleinputvector,firstfailure_Puck[:,0]/1E6, linewidth = 1,  color='blue', linestyle = '--', label = 'Puck')
    # plt.legend()
    # plt.xticks(fontsize = 16)
    # plt.yticks(fontsize = 16)
    # plt.grid(True,alpha=0.5)
    # plt.xlabel(r'Load-Angle $\phi = \frac{N_s}{N_y}$ [deg]', fontsize = 16)
    # plt.ylabel(r"$\sigma_x$' [MPa]", fontsize = 16)
    # plt.xlim(np.min(angleinputvector), np.max(angleinputvector))
    
    
    ###### PLOT FPF STRESS STRAIN ######
    # fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig = plt.figure(2, figsize=(15, 10))
    plt.clf()
    axes = fig.subplots(2, 3)
    # sigma_x'
    ax = axes[0, 0]
    ax.plot(angleinputvector, firstfailure_MaxStress[:,0]/1E6, linewidth=1, color='red', marker='x', markersize=5, label='Max Stress')
    ax.plot(angleinputvector, firstfailure_Puck[:,0]/1E6, linewidth=1, color='blue', linestyle='--', label='Puck')
    ax.legend()
    #ax.set_xlabel(r'Load-Angle $\phi = \frac{N_s}{N_y}$ [deg]', fontsize=16)
    ax.set_ylabel(r"$\sigma_x'$ [MPa]")
    ax.set_xlim(np.min(angleinputvector), np.max(angleinputvector))
    
    # sigma_y'
    ax = axes[0, 1]
    ax.plot(angleinputvector, firstfailure_MaxStress[:,1]/1E6, linewidth=1, color='red', marker='x', markersize=5, label='Max Stress')
    ax.plot(angleinputvector, firstfailure_Puck[:,1]/1E6, linewidth=1, color='blue', linestyle='--', label='Puck')
    ax.legend()
    #ax.set_xlabel(r'Load-Angle $\phi = \frac{N_s}{N_y}$ [deg]', fontsize=16)
    ax.set_ylabel(r"$\sigma_y'$ [MPa]")
    ax.set_xlim(np.min(angleinputvector), np.max(angleinputvector))
    
    # sigma_s'
    ax = axes[0, 2]
    ax.plot(angleinputvector, firstfailure_MaxStress[:,2]/1E6, linewidth=1, color='red', marker='x', markersize=5, label='Max Stress')
    ax.plot(angleinputvector, firstfailure_Puck[:,2]/1E6, linewidth=1, color='blue', linestyle='--', label='Puck')
    ax.legend()
    #ax.set_xlabel(r'Load-Angle $\phi = \frac{N_s}{N_y}$ [deg]', fontsize=16)
    ax.set_ylabel(r"$\sigma_s'$ [MPa]")
    ax.set_xlim(np.min(angleinputvector), np.max(angleinputvector))
    
    # e_x
    ax = axes[1, 0]
    ax.plot(angleinputvector, firstfailure_globalstrain_MaxStress[:,0], linewidth=1, color='red', marker='x', markersize=5, label='Max Stress')
    ax.plot(angleinputvector, firstfailure_globalstrain_Puck[:,0], linewidth=1, color='blue', linestyle='--', label='Puck')
    ax.legend()
    ax.set_xlabel(r'Load-Angle $\phi = \frac{N_s}{N_y}$ [deg]')
    ax.set_ylabel(r"$\epsilon_x$ [-]")
    ax.set_xlim(np.min(angleinputvector), np.max(angleinputvector))
    
    # e_y
    ax = axes[1, 1]
    ax.plot(angleinputvector, firstfailure_globalstrain_MaxStress[:,1], linewidth=1, color='red', marker='x', markersize=5, label='Max Stress')
    ax.plot(angleinputvector, firstfailure_globalstrain_Puck[:,1], linewidth=1, color='blue', linestyle='--', label='Puck')
    ax.legend()
    ax.set_xlabel(r'Load-Angle $\phi = \frac{N_s}{N_y}$ [deg]')
    ax.set_ylabel(r"$\epsilon_y$ [-]")
    ax.set_xlim(np.min(angleinputvector), np.max(angleinputvector))

    # e_s
    ax = axes[1, 2]
    ax.plot(angleinputvector, firstfailure_globalstrain_MaxStress[:,2], linewidth=1, color='red', marker='x', markersize=5, label='Max Stress')
    ax.plot(angleinputvector, firstfailure_globalstrain_Puck[:,2], linewidth=1, color='blue', linestyle='--', label='Puck')
    ax.legend()
    ax.set_xlabel(r'Load-Angle $\phi = \frac{N_s}{N_y}$ [deg]')
    ax.set_ylabel(r"$\epsilon_s$ [-]")
    ax.set_xlim(np.min(angleinputvector), np.max(angleinputvector))
    
    plt.tight_layout()
    
    if savefig:
        filename = f'2a_Failure_FPF_Stress_Strain_loadincrement={loadincrement}_angleincrement={angleincrement}.png'
        plt.savefig(filename, dpi=500)
    
        
    ###### PLOT LPF STRESS STRAIN ######
    fig = plt.figure(3, figsize=(15, 10))
    plt.clf()
    axes = fig.subplots(2, 3)
    
    # sigma_x'
    ax = axes[0, 0]
    ax.plot(angleinputvector, lastfailure_MaxStress[:,0]/1E6, linewidth=1, color='red', marker='x', markersize=5, label='Max Stress')
    ax.plot(angleinputvector, lastfailure_Puck[:,0]/1E6, linewidth=1, color='blue', linestyle='--', label='Puck')
    ax.legend()
    ax.set_ylabel(r"$\sigma_x'$ [MPa]")
    ax.set_xlim(np.min(angleinputvector), np.max(angleinputvector))
    
    # sigma_y'
    ax = axes[0, 1]
    ax.plot(angleinputvector, lastfailure_MaxStress[:,1]/1E6, linewidth=1, color='red', marker='x', markersize=5, label='Max Stress')
    ax.plot(angleinputvector, lastfailure_Puck[:,1]/1E6, linewidth=1, color='blue', linestyle='--', label='Puck')
    ax.legend()
    ax.set_ylabel(r"$\sigma_y'$ [MPa]")
    ax.set_xlim(np.min(angleinputvector), np.max(angleinputvector))
    
    # sigma_s'
    ax = axes[0, 2]
    ax.plot(angleinputvector, lastfailure_MaxStress[:,2]/1E6, linewidth=1, color='red', marker='x', markersize=5, label='Max Stress')
    ax.plot(angleinputvector, lastfailure_Puck[:,2]/1E6, linewidth=1, color='blue', linestyle='--', label='Puck')
    ax.legend()
    ax.set_ylabel(r"$\sigma_s'$ [MPa]")
    ax.set_xlim(np.min(angleinputvector), np.max(angleinputvector))
    
    # e_x
    ax = axes[1, 0]
    ax.plot(angleinputvector, lastfailure_globalstrain_MaxStress[:,0], linewidth=1, color='red', marker='x', markersize=5, label='Max Stress')
    ax.plot(angleinputvector, lastfailure_globalstrain_Puck[:,0], linewidth=1, color='blue', linestyle='--', label='Puck')
    ax.legend()
    ax.set_xlabel(r'Load-Angle $\phi = \frac{N_s}{N_y}$ [deg]')
    ax.set_ylabel(r"$\epsilon_x$ [-]")
    ax.set_xlim(np.min(angleinputvector), np.max(angleinputvector))
    
    # e_y
    ax = axes[1, 1]
    ax.plot(angleinputvector, lastfailure_globalstrain_MaxStress[:,1], linewidth=1, color='red', marker='x', markersize=5, label='Max Stress')
    ax.plot(angleinputvector, lastfailure_globalstrain_Puck[:,1], linewidth=1, color='blue', linestyle='--', label='Puck')
    ax.legend()
    ax.set_xlabel(r'Load-Angle $\phi = \frac{N_s}{N_y}$ [deg]')
    ax.set_ylabel(r"$\epsilon_y$ [-]")
    ax.set_xlim(np.min(angleinputvector), np.max(angleinputvector))
    
    # e_s
    ax = axes[1, 2]
    ax.plot(angleinputvector, lastfailure_globalstrain_MaxStress[:,2], linewidth=1, color='red', marker='x', markersize=5, label='Max Stress')
    ax.plot(angleinputvector, lastfailure_globalstrain_Puck[:,2], linewidth=1, color='blue', linestyle='--', label='Puck')
    ax.legend()
    ax.set_xlabel(r'Load-Angle $\phi = \frac{N_s}{N_y}$ [deg]')
    ax.set_ylabel(r"$\epsilon_s$ [-]")
    ax.set_xlim(np.min(angleinputvector), np.max(angleinputvector))
    
    plt.tight_layout()

    if savefig:
        filename = f'2a_Failure_LPF_Stress_Strain_loadincrement={loadincrement}_angleincrement={angleincrement}.png'
        plt.savefig(filename, dpi=500)
    

        
## 2b: Post-Processing (incld. Damage Tolerance) 

###### PLOT Damage Tolerance ######

damage_tol_MaxStress = np.sqrt((lastfailure_MaxStress[:,1]/1E6- firstfailure_MaxStress[:,1]/1E6)**2 +  (lastfailure_MaxStress[:,2]/1E6- firstfailure_MaxStress[:,2]/1E6)**2)
damage_tol_Puck = np.sqrt((lastfailure_Puck[:,1]/1E6- firstfailure_Puck[:,1]/1E6)**2 +  (lastfailure_Puck[:,2]/1E6- firstfailure_Puck[:,2]/1E6)**2)

plt.figure(4)
plt.clf()
plt.plot(angleinputvector, damage_tol_MaxStress, linewidth=1, color='red', marker='x', markersize=5, label='Max Stress')
plt.plot(angleinputvector, damage_tol_Puck, linewidth=1, color='blue', linestyle='--', label='Puck')
plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)
plt.legend(fontsize = 16)
plt.grid(True,alpha=0.5)
plt.xlabel(r'Load-Angle $\phi = \frac{N_s}{N_y}$ [deg]',  fontsize = 16)
plt.ylabel(r"Damage Tolerance [MPa]", fontsize = 16)

if savefig:
    filename = f'2b_DamageTol_loadincrement={loadincrement}_angleincrement={angleincrement}.png'
    plt.savefig(filename, dpi=500)


plt.show()