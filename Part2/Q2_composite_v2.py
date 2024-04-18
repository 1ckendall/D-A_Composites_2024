import numpy as np
import matplotlib.pyplot as plt
from Class import Lamina, Laminate
# from Q1_aluminum_v2 import *


# Composite (UD tape) material properties 
Ex = 142*10**9 # [Pa]
Ey = 11.2*10**9 # [Pa]
Gxy = 5*10**9 # [Pa]
vxy = 0.3
Xt = 2200*10**6 # [Pa]
Xc = 1800*10**6 # [Pa]
Yt = 70*10**6 # [Pa]
Yc = 300*10**6 # [Pa]
S =	100*10**6 # [Pa]
rho  = 1610 # [kg/m^3]
t_ply = 0.135*10**-3 # [m]

# Coordinate system 
# x: pitch (+ve pointing outboard of left wing), 
# y: yaw (+ve pointing upwards), 
# z: roll (+ve pointing towards nose)

# Geometry
D_outer = 6 # [m]
t = 0.002 # [m] # TODO: placeholde
D_inner = D_outer - 2*t # [m]
R_outer = D_outer/2 # [m]
R_inner = D_inner/2 # [m]

# Knockdown/safety factors 
# TODO: find values from slides/literature
pass

# Load case
M_x = 15000000 # [Nm]
V_y = 1500000 # [N]

# define mass
def calc_mass_per_unit_length(t_arr):
    
    theta = np.linspace(0, np.pi/2, len(t_arr)+1)
    y = R_outer*np.sin(theta)

    
    for i in range(len(y)):
        if i > 0:
            theta[i] = np.arcsin(-(theta[0] - 1/R_outer*(y[i]-y[0])))

    arc_length = theta*D_outer/2 # [m]
    mass_panel_segment = (arc_length[1:]-arc_length[:-1])*t_arr*rho # [kg/m]
    mass_unit_length_quarter_fuselage = np.sum(mass_panel_segment)
    mass_unit_length = 4*mass_unit_length_quarter_fuselage
    return mass_unit_length


# compute Ixx. Assume double symmetric section
def calc_second_moment_of_area(t_arr):
    
    theta = np.linspace(0, np.pi/2, len(t_arr)+1)    
    Ixx_arr = np.zeros(len(t_arr))
    
    for i in range(len(t_arr)):
        Ixx_arr[i] = t_arr[i]*R_outer**3*(0.5*(theta[i+1]-theta[i])-0.25*(np.sin(2*theta[i+1]) - np.sin(2*theta[i])))

    Ixx = 4*np.sum(Ixx_arr)
    return Ixx



# Loading
#def calc_loadvector(V, M, t_arr, n_points=5):
def calc_loadvector(V, M, t_arr):

    theta = np.linspace(0, np.pi/2, len(t_arr)+1)
    y = R_outer*np.sin(theta)
    
    
    
    # sigma_z_times_t_arr = np.zeros(len(t_arr)+1)
    # qs_arr = np.zeros(len(t_arr)+1)
    sigma_z_times_t_arr = np.zeros(len(t_arr))
    qs_arr = np.zeros(len(t_arr))
    
    Ixx = calc_second_moment_of_area(t_arr)
    
    # computed at extemetities, conservative
    for i in range(len(t_arr)):
        # bending stress taken at the top of panel: y_i = max[y_{i+1}, y_{i}]
        sigma_z_times_t_arr[i] = M/Ixx*y[i+1]*t_arr[i]
        # shear flow taken at the bottom of panel: y_i = min[y_{i+1}, y_{i}]
        qs_arr[i] = V/Ixx*R_outer**2*np.cos(theta[i])*t_arr[i]
#        qs_arr[i] = V/(np.pi*R_outer)*np.sqrt(1-(y[i]/R_outer)**2)
        # to verify aranshu's code:
        test_sigma_z = M/Ixx*y[i+1]
        print(f'test_sigma_z:  {test_sigma_z*10**-6} MPa\nTest shear flow {qs_arr[i]} [N/m]')

    return sigma_z_times_t_arr, qs_arr

# check for buckling #TODO:which type of buckling (local? skin? conservative or not?)
def Buckling_check(Nx,Ny,Ns,Dmatrix,a = 1,b = 1,m = 1):
     AR = a/b
     R_buckling = 0
     #for compressive loading 
     N_0 = ((np.pi**2) *(Dmatrix[0,0]*m**4 + 2*(Dmatrix[0,1]+2*Dmatrix[2,2])*(m**2)*(AR**2)+Dmatrix[1,1]*(AR**4) ))/((a**2)*(m**2))
     #shear buckling: 
     beta = (Dmatrix[0,0]/Dmatrix[1,1])**(1/4)
     A = -0.27 + 0.185 *((Dmatrix[0,1]+2*Dmatrix[2,2])/(np.sqrt(Dmatrix[0,0]*Dmatrix[1,1])))
     
     
     B =0.82 + 0.46*((Dmatrix[0,1]+2*Dmatrix[2,2])/(np.sqrt(Dmatrix[0,0]*Dmatrix[1,1]))) -0.2*((Dmatrix[0,1]+2*Dmatrix[2,2])/(np.sqrt(Dmatrix[0,0]*Dmatrix[1,1])))**2
     K = 8.2 + 5 * ((Dmatrix[0,1]+2*Dmatrix[2,2])/(Dmatrix[0,0]*Dmatrix[1,1]))*(1/(10**(A/beta + B*beta)))
     Nxy = 4/(b**2) *( (Dmatrix[0,0]*Dmatrix[1,1]**3)**(1/4)) * K 
     print('buckling',N_0,Nxy)
     #check for compressive loading: 
     R_c = 0
     if  Nx < 0: 
         R_c = Nx / N_0 
     elif Ny < 0 : 
         R_c = Ny / N_0
     #for shear loading: 
     R_s = np.abs(Ns) / Nxy 
     R_buckling = R_c + R_s**2 
     if R_buckling >=1 : 
         print('buckling has occured')
         print('shear ratio',R_s**2)
         print('Compressive ration',R_c)
     elif R_buckling <1: 
         print('no buckling')
         

def is_MaxStress_FPF(Laminate):
    isFPF = False
    ply_fail_mode = None
    if not isFPF:
        for idx, ply in enumerate(Laminate.plys):
            
            ply.sigma_1 = Laminate.sigma11[idx]
            ply.sigma_2 = Laminate.sigma22[idx]
            ply.tau_21 = Laminate.sigma12[idx]
            
            ply.maxStressFibreFail() #fiber failure
            ply.maxStressInterFibreFail() #InterFiberfailure
            ply.maxStressShearFail() #Shearfailure
            if ply.failuremode:
                # ply_fail_sigma = sigma
                ply_fail_mode = ply.failuremode
                plyidx_fail = idx                
                isFPF = True
    return isFPF, ply_fail_mode 


# Black Aluminum design
# quasi isotropic laminate

# (1) initial guess: -8plys- total thickness = 1.08 [mm], too thin 
quasi_isotropic_layup = [90, +45, -45, 0, 0, -45, +45, 90] 

# (2) guess after 'black aluminum' design: -16plys- total thickness = 2.16 [mm], inadequate thickness
# quasi_isotropic_layup = [90, +45, -45, 0, 0, -45, +45, 90, 
#                           90, +45, -45, 0, 0, -45, +45, 90] 


# (3) slightly tailored quasi-isotropic, 0 degree bias: -18plys- total thickness = 2.43 [mm], excessive thickness (as putting more 0 deg increases Ex_planar)
# quasi_isotropic_layup = [90, +45, -45, 0, 0,  0, -45, +45, 90, 
#                          90, +45, -45, 0, 0, 0, -45, +45, 90] 

# (4) slightly more tailored quasi-isotropic removed two 90 degree plies, 0 degree bias: -16plys- total thickness = 2.16 [mm], slightly excessive thickness. Computation states that only 14.3 plies => rounded up to 15 plies => rounded up to 16 plies for balanced, symmetric 
# quasi_isotropic_layup = [90, +45, -45, 0, 0, -45, +45, 0, 
#                          0, +45, -45, 0, 0, -45, +45, 90] 

# (5) ever-so-slightly more tailored quasi-isotropic removed one more 90 degree ply, i put the one 90 degree ply in the middle: 0 degree bias: -15plys- total thickness = 2.16 [mm], slightly excessive thickness. Computation states that only 14.3 plies => rounded up to 15 plies => symmetric 
# disadvantage: +- 45 at the top/bottom instead of 90 
# quasi_isotropic_layup = [+45, -45, 0, 0, 0, -45, +45,
#                          90,
#                          +45, -45, 0, 0, 0, -45, +45] 

plylist= [] 
for angle in quasi_isotropic_layup:
    plylist.append(Lamina(angle,E1 = Ex,E2 = Ey,G12 =Gxy,v12=vxy,Xt=Xt,Xc=Xc,Yt=Yt,Yc=Yt,S=S,t=t_ply))
black_aluminum_laminate = Laminate(plylist)

# Engineering Coefficients
# In-plane
Ex_planar = black_aluminum_laminate.Ex
Ey_planar = black_aluminum_laminate.Ey
Gxy_planar = black_aluminum_laminate.Gxy
vxy_planar = black_aluminum_laminate.vxy
vyx_planar = black_aluminum_laminate.vyx

# Flexural
Exb_flexural = black_aluminum_laminate.Exb
Eyb_flexural = black_aluminum_laminate.Eyb
Gxyb_flexural = black_aluminum_laminate.Gxyb
vxyb_flexural = black_aluminum_laminate.vxyb
vyxb_flexural = black_aluminum_laminate.vyxb

print(f'--- START: Black Aluminum Design --- \nQuasi-Isotropic Layup: {quasi_isotropic_layup}')
print("Planar Properties:")
print(f"Ex_planar: {Ex_planar * 10**-9} GPa")
print(f"Ey_planar: {Ey_planar * 10**-9} GPa")
print(f"Gxy_planar: {Gxy_planar * 10**-6} MPa")
print(f"vxy_planar: {vxy_planar} MPa")
print(f"vyx_planar: {vyx_planar} MPa")

print("Flexural Properties:")
print(f"Exb_flexural: {Exb_flexural * 10**-9} GPa")
print(f"Eyb_flexural: {Eyb_flexural * 10**-9} GPa")
print(f"Gxyb_flexural: {Gxyb_flexural * 10**-9} GPa")
print(f"vxyb_flexural: {vxyb_flexural} MPa")
print(f"vyxb_flexural: {vyxb_flexural} MPa")


# stiffness matching, first thickness estimate
# (EA)_composite = (EA)_alu and (EI)_composite = (EI)_alu => (Et)_composite = (Et)_alu 
t_alu_arr = np.array([0.0012407, 0.00189225, 0.00189225, 0.00207727])
mass_unit_length_alu =  calc_mass_per_unit_length(t_alu_arr)


E_alu = 69.*10**9 #[Pa]
t_composite_arr = E_alu/ Ex_planar * t_alu_arr
n_plies = np.ceil(t_composite_arr/t_ply)
print('Stiffness matching of Alu/quasi isotropic composite')
print(f'Alu: {t_alu_arr*10**3}[mm]')
print(f'CFRP: {t_composite_arr*10**3}[mm]')
print(f'#plys: {n_plies}[mm]')
print(f'--- END: Black Aluminum Design --- \n')


# strength based sizing
# baseline laminate
baseline_layup = [+45, -45, 0, 0, 0, -45, +45,
                         90,
                 +45, -45, 0, 0, 0, -45, +45] # conservative, 15 ply design 
t_baseline = t_ply*len(baseline_layup) #t_composite_arr[-1] # exact



# hard-coded: three different thicknesses (varied at the same location as aluminum)
top_layup = baseline_layup
t_top = t_baseline


diagonal_layup = [+45, -45, 0, 0, -45, +45,
                         90,
                 +45, -45, 0, 0, -45, +45]
t_diagonal =  t_ply*len(diagonal_layup) # close enough to # t_composite_arr[-2] (dropped 2 plies- for symmetry instead of 3- optimal)


middle_layup = [+45, -45, 0, -45, +45,
                         90,
                 +45, -45, 0, -45, +45]
t_middle =  t_ply*len(middle_layup) # close enough to # t_composite_arr[-2] (dropped 4 plies- for symmetry instead of 5- optimal)


# update t_composite_arr based on above laminates
t_composites_arr = np.array([t_middle, t_diagonal, t_diagonal, t_top])


# initialise loadvector for ABD
n_x_arr, n_s_arr = calc_loadvector(V = V_y, M = M_x, t_arr = t_composite_arr)
# n_x_arr = np.zeros(4)
# n_s_arr = np.zeros(4)
n_y = 0
m_x = 0
m_y = 0
m_s = 0



# create required laminates
plylist= [] 
for angle in baseline_layup:
    plylist.append(Lamina(angle,E1 = Ex,E2 = Ey,G12 =Gxy,v12=vxy,Xt=Xt,Xc=Xc,Yt=Yt,Yc=Yt,S=S,t=t_ply))
baseline_laminate = Laminate(plylist,  Nx=n_x_arr[3], Ny=0, Ns=n_s_arr[0], Mx=0, My=0, Ms=0) # worst-case

plylist= [] 
for angle in top_layup:
    plylist.append(Lamina(angle,E1 = Ex,E2 = Ey,G12 =Gxy,v12=vxy,Xt=Xt,Xc=Xc,Yt=Yt,Yc=Yt,S=S,t=t_ply))
# top_laminate = Laminate(plylist,  Nx=n_x_arr[3], Ny=0, Ns=n_s_arr[2], Mx=0, My=0, Ms=0)
top_laminate = Laminate(plylist,  Nx=n_x_arr[3], Ny=0, Ns=n_s_arr[3], Mx=0, My=0, Ms=0)

plylist= [] 
for angle in diagonal_layup:
    plylist.append(Lamina(angle,E1 = Ex,E2 = Ey,G12 =Gxy,v12=vxy,Xt=Xt,Xc=Xc,Yt=Yt,Yc=Yt,S=S,t=t_ply))
#diagonal_laminate = Laminate(plylist,  Nx=n_x_arr[2], Ny=0, Ns=n_s_arr[1], Mx=0, My=0, Ms=0)
diagonal_laminate = Laminate(plylist,  Nx=n_x_arr[2], Ny=0, Ns=n_s_arr[1], Mx=0, My=0, Ms=0) # conservative

plylist= [] 
for angle in middle_layup:
    plylist.append(Lamina(angle,E1 = Ex,E2 = Ey,G12 =Gxy,v12=vxy,Xt=Xt,Xc=Xc,Yt=Yt,Yc=Yt,S=S,t=t_ply))
# middle_laminate = Laminate(plylist,  Nx=n_x_arr[1], Ny=0, Ns=n_s_arr[0], Mx=0, My=0, Ms=0)
middle_laminate = Laminate(plylist,  Nx=n_x_arr[0], Ny=0, Ns=n_s_arr[0], Mx=0, My=0, Ms=0)

# CLPT stress-strain
baseline_laminate.getStressStrain()
top_laminate.getStressStrain()
diagonal_laminate.getStressStrain()
middle_laminate.getStressStrain()

# get stiffness
Ex_planar_top = top_laminate.Ex
Ex_planar_diagonal = diagonal_laminate.Ex
Ex_planar_middle = middle_laminate.Ex


# Retrieve stress vectors in principal coordinate system
baseline_stress_vector = baseline_laminate.localstressVector

top_stress_vector = top_laminate.localstressVector
top_stress_sigma11 = top_laminate.sigma11
top_stress_sigma22 = top_laminate.sigma22
top_stress_sigma12 = top_laminate.sigma12



diagonal_stress_vector = diagonal_laminate.localstressVector
diagonal_stress_sigma11 = diagonal_laminate.sigma11
diagonal_stress_sigma22 = diagonal_laminate.sigma22
diagonal_stress_sigma12 = diagonal_laminate.sigma12

middle_stress_vector = middle_laminate.localstressVector
middle_stress_sigma11 = middle_laminate.sigma11
middle_stress_sigma22 = middle_laminate.sigma22
middle_stress_sigma12 = middle_laminate.sigma12

# Print the stress vectors
# print("\nBaseline Laminate Stress Vector:", baseline_stress_vector*10**-6, 'MPa')
# print("\nTop Laminate Stress Vector:", top_stress_vector*10**-6, 'MPa')
# print("\nDiagonal Laminate Stress Vector:", diagonal_stress_vector*10**-6, 'MPa')
# print("\nMiddle Laminate Stress Vector:", middle_stress_vector*10**-6, 'MPa')


# Print stress components in MPa
print("Top Laminate Stress Components (MPa):")
print(f'Ply angles: {baseline_layup} [deg]')
print(f'Ex Planar: {Ex_planar_top * 10**-9} GPa')
print(f"sigma11: {top_stress_sigma11 * 10**-6} MPa")
print(f"sigma22: {top_stress_sigma22 * 10**-6} MPa")
print(f"sigma12: {top_stress_sigma12 * 10**-6} MPa")
print()

print("Diagonal Laminate Stress Components (MPa):")
print(f'Ply angles: {diagonal_layup} [deg]')
print(f'Ex Planar: {Ex_planar_diagonal * 10**-9} GPa')
print(f"sigma11: {diagonal_stress_sigma11 * 10**-6} MPa")
print(f"sigma22: {diagonal_stress_sigma22 * 10**-6} MPa")
print(f"sigma12: {diagonal_stress_sigma12 * 10**-6} MPa")
print()

print("Middle Laminate Stress Components (MPa):")
print(f'Ply angles: {middle_layup} [deg]')
print(f'Ex Planar: {Ex_planar_middle * 10**-9} GPa')
print(f"sigma11: {middle_stress_sigma11 * 10**-6} MPa")
print(f"sigma22: {middle_stress_sigma22 * 10**-6} MPa")
print(f"sigma12: {middle_stress_sigma12 * 10**-6} MPa")


# part A: check max stress condition on given laminates. Plot reserve factors.
isFPF_top = is_MaxStress_FPF(top_laminate)[0]
isFPF_diagonal = is_MaxStress_FPF(diagonal_laminate)[0]
isFPF_middle = is_MaxStress_FPF(middle_laminate)[0]

print('FPF Check')
print(f'Top Laminate FPF?: {isFPF_top}')
print(f'Diagonal Laminate FPF?: {isFPF_diagonal}')
print(f'Middle Laminate FPF?: {isFPF_middle}')

# part B: check for buckling

print('\nTop panel buckling check')
Buckling_check(Nx = top_laminate.Nx, Ny = top_laminate.Ny, Ns = top_laminate.Ns, Dmatrix = top_laminate.D)

print('\nDiagonal panel buckling check')
Buckling_check(Nx = diagonal_laminate.Nx, Ny = diagonal_laminate.Ny, Ns = diagonal_laminate.Ns, Dmatrix = diagonal_laminate.D)

print('\nMiddle panel buckling check')
Buckling_check(Nx = middle_laminate.Nx, Ny = middle_laminate.Ny, Ns = middle_laminate.Ns, Dmatrix = middle_laminate.D)




# t_var_composite_arr = t_ply*np.array([len(middle_layup), len(diagonal_layup), len(top_layup)])
var_cfrp_mass = calc_mass_per_unit_length(t_composite_arr)
print(f'\nvar_cfrp_mass: {var_cfrp_mass} [kg/m]')