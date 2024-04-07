import matplotlib.pyplot as plt
import numpy as np
from pprint import pprint
from Class import Lamina, Laminate

# plot controls (Booleans)
plot_2a= True # plot 2a: failure envelopes (default = True)
plot_2b= True # plot 2b: damage tolerance (default = True)


# run controls (Booleans)
run_maxStress = True # run failure envelope for max stress (default = True)
run_Puck = True # run failure envelope for Puck (default = True)


# savefig controls (Booleans)
savefig = True # default = False


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

stressinputvector = np.arange(1E3,1E8,5000) ## UNCOMMENT FOR SHORT RUN
angleinputvector = np.arange(0,375,15)  ## UNCOMMENT FOR SHORT RUN

# stressinputvector = np.arange(1E3,1E8,500) ## UNCOMMENT FOR FINAL RUN
# angleinputvector = np.arange(0,361,1) ## UNCOMMENT FOR FINAL RUN

## 2ai: Failure Theory of Maximum Stress 
firstfailuremaxstress = np.zeros((len(angleinputvector),2))
firstfaileglobalstrain =np.zeros((len(angleinputvector),2))

lastplyfailuremaxstress=np.zeros((len(angleinputvector),2))
lastplygolbalstrain = np.zeros((len(angleinputvector),2))

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
                                firstfailureloadY =LaminateQ3.sigmayprime
                                firstfailureloadS = LaminateQ3.sigmaxyprime
                                
                                firstfailuremaxstress[i,0] = firstfailureloadY
                                firstfailuremaxstress[i,1] = firstfailureloadS
                                firstfaileglobalstrain[i,0] = LaminateQ3.eyglobal
                                firstfaileglobalstrain[i,1] = LaminateQ3.esglobal
                                print(f'First-Ply Failure [FPF]: Ply index: {k}, Fibre angle: {laminaangle}, \nFailure Mode: {failurechecking.failuremode}, \nLoadY: {firstfailureloadY}, LoadS: {firstfailureloadS}\ney: {firstfaileglobalstrain[i,0]}, es: {firstfaileglobalstrain[i,1]}\n')
                        
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
                                        lastfailureloadY =LaminateQ3.sigmayprime
                                        lastfailureloadS = LaminateQ3.sigmaxyprime
                                        lastplyfailuremaxstress[i,0] = lastfailureloadY
                                        lastplyfailuremaxstress[i,1] =lastfailureloadS 
                                        lastplygolbalstrain[i,0] = LaminateQ3.eyglobal
                                        lastplygolbalstrain[i,1] = LaminateQ3.esglobal
                                        print(f'Last-Ply Failure [LPF]: Ply index: {k}, Fibre angle: {laminaangle}, \nFailure Mode: {failurechecking.failuremode}, LoadY: {lastfailureloadY}, LoadS: {lastfailureloadS}\ney: {lastplygolbalstrain[i,0]}, es: {lastplygolbalstrain[i,1]}\n\n')
                                    
                                        assert np.allclose(failuretracking, 2)
            
    if plot_2a:
        # making failure envelope plot (maxstress) sigmay-sigmas  
        plt.figure(1)
        plt.clf()
        plt.plot(firstfailuremaxstress[:,0]/1E6,firstfailuremaxstress[:,1]/1E6, linewidth = 1,  color='red', label = 'FPF')
        plt.plot(lastplyfailuremaxstress[:,0]/1E6,lastplyfailuremaxstress[:,1]/1E6, linewidth = 1, color='blue', label ='LPF')
        plt.xticks(fontsize = 16)
        plt.yticks(fontsize = 16)
        plt.legend(fontsize = 16)
        plt.grid(True,alpha=0.5)
        plt.xlabel(r'$\sigma_y$ [MPa]', fontsize = 16)
        plt.ylabel(r'$\sigma_s$ [MPa]', fontsize = 16)
        
        if savefig:
            loadincrement = stressinputvector[1]-stressinputvector[0]
            angleincrement = angleinputvector[1]-angleinputvector[0]
            filename = f'2a_Failure_Envelope_Stress_Max_Stressloadincrement={loadincrement}_angleincrement={angleincrement}.png'
            plt.savefig(filename, dpi=500)
            
        # making failure envelope plot (maxstress) ey-es
        plt.figure(2)
        plt.clf()
        plt.plot(firstfaileglobalstrain[:,0],firstfaileglobalstrain[:,1], linewidth = 1,  color='red', label = 'FPF')
        plt.plot(lastplygolbalstrain[:,0],lastplygolbalstrain[:,1], linewidth = 1, color='blue', label ='LPF')
        plt.xticks(fontsize = 16)
        plt.yticks(fontsize = 16)
        plt.legend(fontsize = 16)
        plt.grid(True,alpha=0.5)
        plt.xlabel(r'$\epsilon_y$ [-]', fontsize = 16)
        plt.ylabel(r'$\epsilon_s$ [-]', fontsize = 16)
        plt.show()
        
        if savefig:
            loadincrement = stressinputvector[1]-stressinputvector[0]
            angleincrement = angleinputvector[1]-angleinputvector[0]
            filename = f'2a_Failure_Envelope_Strain_Max_Stress_loadincrement={loadincrement}_angleincrement={angleincrement}.png'
            plt.savefig(filename, dpi=500)


## 2ai: Failure Theory of Puck
firstfaileglobalstrain =np.zeros((len(angleinputvector),2))
firstfailure_Puck = np.zeros((len(angleinputvector),2))

lastplygolbalstrain = np.zeros((len(angleinputvector),2))
lastfailure_Puck = np.zeros((len(angleinputvector),2))
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
                                firstfailureloadY =LaminateQ3.sigmayprime
                                firstfailureloadS = LaminateQ3.sigmaxyprime
                                
                                firstfailure_Puck[i,0] = firstfailureloadY
                                firstfailure_Puck[i,1] = firstfailureloadS
                                firstfaileglobalstrain[i,0] = LaminateQ3.eyglobal
                                firstfaileglobalstrain[i,1] = LaminateQ3.esglobal
                                print(f'First-Ply Failure [FPF]: Ply index: {k}, Fibre angle: {laminaangle}, \nFailure Mode: {failurechecking.failuremode}, LoadY: {firstfailureloadY}, LoadS: {firstfailureloadS}\ney: {firstfaileglobalstrain[i,0]}, es: {firstfaileglobalstrain[i,1]}\n')
                        
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
                                        lastfailureloadY =LaminateQ3.sigmayprime
                                        lastfailureloadS = LaminateQ3.sigmaxyprime
                                        lastfailure_Puck[i,0] = lastfailureloadY
                                        lastfailure_Puck[i,1] =lastfailureloadS 
                                        lastplygolbalstrain[i,0] = LaminateQ3.eyglobal
                                        lastplygolbalstrain[i,1] = LaminateQ3.esglobal
                                        print(f'Last-Ply Failure [LPF]: Ply index: {k}, Fibre angle: {laminaangle}, \nFailure Mode: {failurechecking.failuremode}, LoadY: {lastfailureloadY}, LoadS: {lastfailureloadS}\ney: {lastplygolbalstrain[i,0]}, es: {lastplygolbalstrain[i,1]}\n\n')
                                        
                                        # print(f'Failure index (all plies): {failuretracking}')
                                        assert np.allclose(failuretracking, 2)
    
    if plot_2a:
        plt.figure(3)
        plt.clf()
        plt.plot(firstfailure_Puck[:,0]/1E6,firstfailure_Puck[:,1]/1E6, linewidth = 1,  color='red', label = 'FPF')
        plt.plot(lastfailure_Puck[:,0]/1E6,lastfailure_Puck[:,1]/1E6, linewidth = 1, color='blue', label ='LPF')
        plt.xticks(fontsize = 16)
        plt.yticks(fontsize = 16)
        plt.legend(fontsize = 16)
        plt.grid(True,alpha=0.5)
        plt.xlabel(r'$\sigma_y$ [MPa]', fontsize = 16)
        plt.ylabel(r'$\sigma_s$ [MPa]', fontsize = 16)
        
        if savefig:
            loadincrement = stressinputvector[1]-stressinputvector[0]
            angleincrement = angleinputvector[1]-angleinputvector[0]
            filename = f'2a_Failure_Envelope_Stress_Puck_loadincrement={loadincrement}_angleincrement={angleincrement}.png'
            plt.savefig(filename, dpi=500)
          
        
        plt.figure(4)
        plt.clf()
        plt.plot(firstfaileglobalstrain[:,0],firstfaileglobalstrain[:,1], linewidth = 1,  color='red', label = 'FPF')
        plt.plot(lastplygolbalstrain[:,0],lastplygolbalstrain[:,1], linewidth = 1, color='blue', label ='LPF')
        plt.xticks(fontsize = 16)
        plt.yticks(fontsize = 16)
        plt.legend(fontsize = 16)
        plt.grid(True,alpha=0.5)
        plt.xlabel(r'$\epsilon_y$ [-]', fontsize = 16)
        plt.ylabel(r'$\epsilon_s$ [-]', fontsize = 16)
        plt.show()
        
        if savefig:
            loadincrement = stressinputvector[1]-stressinputvector[0]
            angleincrement = angleinputvector[1]-angleinputvector[0]
            filename = f'2a_Failure_Envelope_Strain_Puck_loadincrement={loadincrement}_angleincrement={angleincrement}.png'
            plt.savefig(filename, dpi=500)


## 2b: Post-Processing (incld. Damage Tolerance) 