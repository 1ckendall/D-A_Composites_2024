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
def calc_loadvector(V, M, n_points=4):  
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
# quasi_isotropic_layup = [90, +45, -45, 0, 0, -45, +45, 90] 

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
quasi_isotropic_layup = [+45, -45, 0, 0, 0, -45, +45,
                         90,
                         +45, -45, 0, 0, 0, -45, +45] 

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
                 +45, -45, 0, 0, -45, +45,]
t_diagonal =  t_ply*len(diagonal_layup) # close enough to # t_composite_arr[-2] (dropped 2 plies- for symmetry instead of 3- optimal)


middle_layup = [+45, -45, 0, -45, +45,
                         90,
                 +45, -45, 0, -45, +45,]
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
baseline_laminate.getABD()
top_laminate.getABD()
diagonal_laminate.getABD()
middle_laminate.getABD()


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
    
    # Ixx = np.pi*R_outer**3*t
    # Iyy = Ixx
    
    y = np.linspace(0, R_outer, len(t_arr)+1)
    theta = np.zeros(len(y))

    
    for i in range(len(y)):
        if i > 0:
            theta[i] = np.arcsin(-(theta[0] - 1/R_outer*(y[i]-y[0])))

    arc_length = theta*D_outer/2 # [m]
    mass_panel_segment = (arc_length[1:]-arc_length[:-1])*t_arr*rho # [kg/m]
    print(f'arc len: {arc_length}, mass_penl: {mass_panel_segment}')
    mass_unit_length_quarter_fuselage = np.sum(mass_panel_segment)
    mass_unit_length = 4*mass_unit_length_quarter_fuselage
    return y, t_arr, mass_unit_length

t_var_composite_arr = t_ply*np.array([len(middle_layup), len(diagonal_layup), len(top_layup)])
var_cfrp_mass = calc_mass_per_unit_length(t_var_composite_arr)
print(f'var_cfrp_mass: {var_cfrp_mass} [kg/m]')



#checkling of buckling

def Buckling_check(Nx,Ny,Ns,Dmatrix,a,b,m):
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


# baseline laminate 
baselinelaminate_buckling = Buckling_check(Nx=n_x_arr[3],Ny=0,Ns=n_s_arr[0],Dmatrix=baseline_laminate.D,a=1,b=1,m=1)
toplam = Buckling_check(Nx=n_x_arr[3],Ny=0,Ns=n_s_arr[2],Dmatrix=top_laminate.D,a=1,b=1,m=1)
diagonal = Buckling_check(Nx=n_x_arr[2],Ny=0,Ns=n_s_arr[1],Dmatrix=diagonal_laminate.D,a=1,b=1,m=1)
middle = Buckling_check(Nx=n_x_arr[1],Ny=0,Ns=n_s_arr[0],Dmatrix=middle_laminate.D,a=1,b=1,m=1)


#panel checking function of buckling: 
#the skin is designed against buckling and the stiffeners stiffness is higher then the ultimate load of 1.5*skin buckling load 
def stiffened_panel_buckling(Ftot,shear_panel,ds,a,b,EA,EI,Amatrix_skin,Dmatrix_skin,t):
    #inputs 
    #Ftot the total force on the panel 
    #spacing between stiffeners
    #a  panel length 
    #b panel width 
    #EA stiffness of stiffener E*A KEEP IN MIND
    #EI of stiffener E*I KEEP IN MINd
    #A_matrix  of the panel 
    #D_matrix of the panel
    #t skin thickness
    #only check that skin or stiffener buckle seperately
    buckling = False

    #load distribution to the skin in between stiffeners 
    F_skin = ((Amatrix_skin[0,0])/(Amatrix_skin[0,0]+(EA)/ds) )* Ftot
    # Nskin load on the panel
    N_skin = F_skin / b
    N_panel = Ftot /b
    #shear load per skin 
    tau = shear_panel/ (a*t)
    #k* used to calculate the number of halfwaces 
    k_star = (Dmatrix_skin[1,1]/Dmatrix_skin[0,0])**(0.25) * (a/ds)
    k = int(k_star) +1  
    if k<1: 
        k = 1
    AR_bar = a/ds
    AR = a/b
    # m* used for panel- to check the skin buckles first and the panel doesnt 
    m_star = (Dmatrix_skin[1,1]/(Dmatrix_skin[0,0]+EI/ds))**(0.25) * (a/b)
    m = int(m_star)
    if m < 1 : 
        m =1
    Nx_panel = (((np.pi)**2)/(a**2)) * (Dmatrix_skin[0,0] * (m**2) + 2*(Dmatrix_skin[0,1]+2*Dmatrix_skin[2,2])*(AR_bar**2)+ Dmatrix_skin[1,1]*((AR_bar**4))/(m**2))

    Nx_skin = (((np.pi)**2)/(a**2)) * (Dmatrix_skin[0,0] * (m**2) + 2*(Dmatrix_skin[0,1]+2*Dmatrix_skin[2,2])*(AR_bar**2)+ Dmatrix_skin[1,1]*((AR_bar**4))/(m**2))
    
    #checking for shear on the skins
    beta = (Dmatrix_skin[0,0]/Dmatrix_skin[1,1])**(1/4)
    A = -0.27 + 0.185 *((Dmatrix_skin[0,1]+2*Dmatrix_skin[2,2])/(np.sqrt(Dmatrix_skin[0,0]*Dmatrix_skin[1,1])))
     
     
    B =0.82 + 0.46*((Dmatrix_skin[0,1]+2*Dmatrix_skin[2,2])/(np.sqrt(Dmatrix_skin[0,0]*Dmatrix_skin[1,1]))) -0.2*((Dmatrix_skin[0,1]+2*Dmatrix_skin[2,2])/(np.sqrt(Dmatrix_skin[0,0]*Dmatrix_skin[1,1])))**2
    K = 8.2 + 5 * ((Dmatrix_skin[0,1]+2*Dmatrix_skin[2,2])/(Dmatrix_skin[0,0]*Dmatrix_skin[1,1]))*(1/(10**(A/beta + B*beta)))
    Nxy = 4/(ds**2) *( (Dmatrix_skin[0,0]*Dmatrix_skin[1,1]**3)**(1/4)) * K 
    R_xy = np.abs(shear_panel) / Nxy
    



    lambda_buckling = (Amatrix_skin[0,0]+(EA)/ds) / Amatrix_skin[0,0] 
    Nx_bay = N_skin / lambda_buckling
    Nx_baycritical = Nx_skin / lambda_buckling
    
    #R_x_skin = N_skin / Nx_skin 
    R_x_skin = Nx_bay / Nx_baycritical


    EIadvised = Dmatrix_skin[0,0] * ds *((np.sqrt(Dmatrix_skin[1,1]/Dmatrix_skin[0,0]))*(2*lambda_buckling*(AR_bar**2)- (np.sqrt(Dmatrix_skin[1,1]/Dmatrix_skin[0,0]))*(AR**4))+((2*(Dmatrix_skin[0,1]+2*Dmatrix_skin[2,2]))/Dmatrix_skin[0,0])*(lambda_buckling*(AR_bar**2)-(AR**2))-1)
    EIneeded = 1.5 * EIadvised # this calculated what EI for stringer we would need such that we make sure the skin buckles first 
    #checking for shear buckling of the skin ?
    R_buckling = R_x_skin + R_xy**2 
    
    print('EIneeded',EIneeded)
    
    R_buckling = R_x_skin + R_xy**2 
    if R_buckling >=1 : 
         print('buckling has occured combined','Rx_skin',R_x_skin,'Rshear',R_xy**2,'combined',R_buckling)
        
         buckling = True
         
    elif R_buckling <1: 
         
         print(' NO buckling ','Rx_skin',R_x_skin,'Rshear',R_xy**2,'Combined',R_buckling)
    return Nx_skin,buckling,EIneeded





