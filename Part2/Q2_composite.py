import numpy as np
import matplotlib.pyplot as plt
from Class import Lamina, Laminate
from Q1_aluminum import *


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
t = 0.002 # [m] TODO: placeholde
D_inner = D_outer - 2*t # [m]
R_outer = D_outer/2 # [m]
R_inner = D_inner/2 # [m]

# Knockdown/safety factors 
# TODO: find values from slides/literature
pass

# Load case
M_x = 15000000 # [Nm]
V_y = 1500000 # [N]


# Loading
def calc_loadvector(V, M, n_points=3):  
    y = np.linspace(0, R_outer, n_points)    
    sigma_z_times_t_arr = np.zeros(len(y))
    qs_arr = np.zeros(len(y))
    for i in range(len(y)):
        sigma_z_times_t_arr[i] = M*y[i]/(np.pi*R_outer**3)
        qs_arr[i] = V/(np.pi*R_outer)*np.sqrt(1-(y[i]/R_outer)**2)
    return sigma_z_times_t_arr, qs_arr


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

print(f'Black Aluminum Design --- Quasi-Isotropic Layup: {quasi_isotropic_layup}')
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
t_alu_arr = calc_variable_thickness(V = V_y, M = M_x, n_points=4)[1]
E_alu = 69.*10**9 #[Pa]
t_composite_arr = E_alu/ Ex_planar * t_alu_arr
n_plies = np.ceil(t_composite_arr/t_ply)
print('Stiffness matching of Alu/quasi isotropic composite')
print(f'Alu: {t_alu_arr*10**3}[mm]')
print(f'CFRP: {t_composite_arr*10**3}[mm]')
print(f'#plys: {n_plies}[mm]')


# strength based sizing
# initialise loadvector for ABD
n_x_arr, n_s_arr = calc_loadvector(V = V_y, M = M_x, n_points=4)
n_y = 0
m_x = 0
m_y = 0
m_s = 0

# baseline laminate
baseline_layup = [+45, -45, 0, 0, 0, -45, +45,
                         90,
                 +45, -45, 0, 0, 0, -45, +45] # conservative, 15 ply design 
t_baseline = t_composite_arr[-1] # exact



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


# create required laminates
plylist= [] 
for angle in baseline_layup:
    plylist.append(Lamina(angle,E1 = Ex,E2 = Ey,G12 =Gxy,v12=vxy,Xt=Xt,Xc=Xc,Yt=Yt,Yc=Yt,S=S,t=t_ply))
baseline_laminate = Laminate(plylist,  Nx=n_x_arr[3], Ny=0, Ns=n_s_arr[0], Mx=0, My=0, Ms=0) # worst-case

plylist= [] 
for angle in top_layup:
    plylist.append(Lamina(angle,E1 = Ex,E2 = Ey,G12 =Gxy,v12=vxy,Xt=Xt,Xc=Xc,Yt=Yt,Yc=Yt,S=S,t=t_ply))
top_laminate = Laminate(plylist,  Nx=n_x_arr[3], Ny=0, Ns=n_s_arr[2], Mx=0, My=0, Ms=0)

plylist= [] 
for angle in diagonal_layup:
    plylist.append(Lamina(angle,E1 = Ex,E2 = Ey,G12 =Gxy,v12=vxy,Xt=Xt,Xc=Xc,Yt=Yt,Yc=Yt,S=S,t=t_ply))
diagonal_laminate = Laminate(plylist,  Nx=n_x_arr[2], Ny=0, Ns=n_s_arr[1], Mx=0, My=0, Ms=0)

plylist= [] 
for angle in middle_layup:
    plylist.append(Lamina(angle,E1 = Ex,E2 = Ey,G12 =Gxy,v12=vxy,Xt=Xt,Xc=Xc,Yt=Yt,Yc=Yt,S=S,t=t_ply))
middle_laminate = Laminate(plylist,  Nx=n_x_arr[1], Ny=0, Ns=n_s_arr[0], Mx=0, My=0, Ms=0)

# CLPT stress-strain
baseline_laminate.getStressStrain()
top_laminate.getStressStrain()
diagonal_laminate.getStressStrain()
middle_laminate.getStressStrain()


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
print(f"sigma11: {top_stress_sigma11 * 10**-6} MPa")
print(f"sigma22: {top_stress_sigma22 * 10**-6} MPa")
print(f"sigma12: {top_stress_sigma12 * 10**-6} MPa")
print()

print("Diagonal Laminate Stress Components (MPa):")
print(f"sigma11: {diagonal_stress_sigma11 * 10**-6} MPa")
print(f"sigma22: {diagonal_stress_sigma22 * 10**-6} MPa")
print(f"sigma12: {diagonal_stress_sigma12 * 10**-6} MPa")
print()

print("Middle Laminate Stress Components (MPa):")
print(f"sigma11: {middle_stress_sigma11 * 10**-6} MPa")
print(f"sigma22: {middle_stress_sigma22 * 10**-6} MPa")
print(f"sigma12: {middle_stress_sigma12 * 10**-6} MPa")





def calc_mass_per_unit_length(t_arr):
    print('\nComputing mass...\n')
    # Ixx = np.pi*R_outer**3*t
    # Iyy = Ixx
    
    y = np.linspace(0, R_outer, len(t_arr)+1)
    theta = np.zeros(len(y))

    
    for i in range(len(y)):
        if i > 0:
            theta[i] = np.arcsin(-(theta[0] - 1/R_outer*(y[i]-y[0])))

    arc_length = theta*D_outer/2 # [m]
    mass_panel_segment = (arc_length[1:]-arc_length[:-1])*t_arr*rho # [kg/m]
    print(f'theta: {np.degrees(theta)} arc len: {arc_length}, mass_penl: {mass_panel_segment}')
    mass_unit_length_quarter_fuselage = np.sum(mass_panel_segment)
    mass_unit_length = 4*mass_unit_length_quarter_fuselage
    return y, t_arr, mass_unit_length

t_var_composite_arr = t_ply*np.array([len(middle_layup), len(diagonal_layup), len(top_layup)])
var_cfrp_mass = calc_mass_per_unit_length(t_var_composite_arr)
print(f'var_cfrp_mass: {var_cfrp_mass} [kg/m]')

