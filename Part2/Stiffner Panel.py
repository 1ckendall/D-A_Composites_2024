from Class import Lamina, Laminate
from bucklingforpanels import stiffened_panel_buckling
from Test_discritization import compute_forces
import numpy as np

# Composite (UD tape) material properties
Ex = 142 * 10 ** 9  # [Pa]
Ey = 11.2 * 10 ** 9  # [Pa]
Gxy = 5 * 10 ** 9  # [Pa]
vxy = 0.3
Xt = 2200 * 10 ** 6  # [Pa]
Xc = 1800 * 10 ** 6  # [Pa]
Yt = 70 * 10 ** 6  # [Pa]
Yc = 300 * 10 ** 6  # [Pa]
S = 100 * 10 ** 6  # [Pa]
rho = 1610  # [kg/m^3]
t_ply = 0.135 * (10 ** -3)  # [m]

# Initializing 
#Laminate Design
Panel_3= top_layup =[+45, -45, +45, -45, 0, 0, -45, +45, -45, +45,
                        90,
                    +45, -45, +45, -45, 0, 0, -45, +45, -45, +45,
                        90
                    +45, -45, +45, -45, 0, 0, -45, +45, -45, +45,
                        90,
                        0, 0, 0, +45, -45, 0, 0, 
                        90, 90, 0, 90, 0, 90, 90,
                        0, 0, -45, 45, 0, 0, 0,
                    +45, -45, +45, -45, 0, 0, -45, +45, -45, +45,
                                90,
                    +45, -45, +45, -45, 0, 0, -45, +45, -45, +45,
                                90,
                    +45, -45, +45, -45, 0, 0, -45, +45, -45, +45]

Panel_3_lamina = [Lamina(i, Ex, Ey, Gxy, vxy, Xt, Xc, Yc, Yt, S, t_ply) for i in Panel_3]
T_p3=t_ply*len(Panel_3)
print('Thickness= ',T_p3)

# Creating Laminate
laminate_3=Laminate(Panel_3_lamina)

# calling A and D matrices
A_mat=laminate_3.A
#print(A_mat)
D_mat=laminate_3.D
Ex = laminate_3.Ex

print('stiffnest value for panel in MPa=', Ex/1e6)

# testing
# skin parameters
t_ply=0.135e-3
Eskin_val=[(45.23*10**9),(60.22*10**9),(71.21*10**9)]
T_skin=[T_p3,T_p3-0.0006,T_p3]

#stiffner parameters
A_stiffner=[0,0.000119,0.000149,0.000179]
Estiffner_val=[0, 48.216e9, 71.265e9, 74.07e9]


force_val,shear_val,weight,sig_stiffner,sig_boom=compute_forces(T_skin,Eskin_val,A_stiffner,Estiffner_val)
#print(shear_val)
print('weight=',weight)
force_panel_1_R=(np.max(force_val[0:9])*9)+(np.sum(force_val[135:144])*9)
force_panel_2_TR=np.max(np.abs(force_val[9:27]))*18
force_panel_3_T=np.max(np.abs(force_val[27:45]))*18
force_panel_2_TL=np.max(np.abs(force_val[45:63]))*18
force_panel_1_L=np.max(np.abs(force_val[63:81]))*18
force_panel_2_CL=np.max(np.abs(force_val[81:99]))*18
force_panel_3_C=np.max(np.abs(force_val[99:117]))*18
force_panel_2_CR=np.max(np.abs(force_val[117:135]))*18

shear_force_1_R=np.sum(shear_val[0:9])+np.sum(shear_val[135:144])
shear_force_2_TR=np.sum(shear_val[9:27])
shear_force_3_T=np.sum(shear_val[27:45])
shear_force_2_TL=np.sum(shear_val[45:63])
shear_force_1_L=np.sum(shear_val[63:81])
shear_force_2_CL=np.sum(shear_val[81:99])
shear_force_3_C=np.sum(shear_val[99:117])
shear_force_2_CR=np.sum(shear_val[117:135])

print(sig_boom)
print(sig_stiffner)

'''
shear_force_1_R=np.max(shear_val[0:9]) and np.max(shear_val[135:144])
shear_force_2_TR=np.max(shear_val[9:27])
shear_force_3_T=np.max(shear_val[27:45])
shear_force_2_TL=np.max(shear_val[45:63])
shear_force_1_L=np.max(shear_val[63:81])
shear_force_2_CL=np.max(shear_val[81:99])
shear_force_3_C=np.max(shear_val[99:117])
shear_force_2_CR=np.max(shear_val[117:135])
'''
print('shear force value on plate 1= ',shear_force_1_R)
print('Normal Force on plate 1=', force_panel_1_R)
print("shear force value on plate 3=", shear_force_3_C)
print("Normal Force", force_panel_3_C)


# total force on panel due to bending
Ftot= 551225.9 # N
shear_panel = 2907502.4922648575 # N/m
ds=0.13089*2 #m
a=0.5 #m
b=2.35619 #m
EA_stiffner= 0.006495 * 1e9
EI_stiffner= 301.4536406876552
T_skin=t_ply*len(Panel_3)

#
bucklepanel_3=stiffened_panel_buckling(Ftot,shear_panel,ds,a,b,EA_stiffner,EI_stiffner,A_mat,D_mat,T_skin)
    #inputs 
    #Ftot the total force on the panel 
    #ds: spacing between stiffeners
    #a panel length 
    #b panel width 
    #EA stiffness of stiffener E*A KEEP IN MIND
    #EI of stiffener E*I KEEP IN MINd
    #A_matrix  of the panel 
    #D_matrix of the panel
    #t skin thickness
    #only check that skin or stiffener buckle seperately
