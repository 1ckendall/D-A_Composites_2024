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
t_ply = 0.135 * 10 ** -3  # [m]

# Initializing 
#Laminate Design
Panel_3=[+45, -45, +45, -45, 0, 0, -45, +45, -45, +45,
                        90,
        +45, -45, +45, -45, 0, 0, -45, +45, -45, +45,
            0, 0, 0, +45, -45, 0, 0, 
            90, 90, 0, 90, 0, 90, 90,
            0, 0, +45, -45, 0, 0, 0,
        +45, -45, +45, -45, 0, 0, -45, +45, -45, +45,
                    90,
        +45, -45, +45, -45, 0, 0, -45, +45, -45, +45]
Panel_3_lamina = [Lamina(i, Ex, Ey, Gxy, vxy, Xt, Xc, Yc, Yt, S, t_ply) for i in Panel_3]
T_p3=t_ply*len(Panel_3)
# Creating Laminate
laminate_3=Laminate(Panel_3_lamina)

# calling A and D matrices
A_mat=laminate_3.A
D_mat=laminate_3.D
Ex = laminate_3.Ex

print('stiffnest value for panel in MPa=', Ex/1e6)

# testing
    # skin parameters
t_ply=0.135e-3
Eskin_val=[(45.23*10**9),(60.22*10**9),(71.21*10**9)]
T_skin=[T_p3-8,T_p3-4,T_p3]

#stiffner parameters
A_stiffner=[0,6.823e-5,0.000179,0.000179]

force_val,shear_val=compute_forces(T_skin,Eskin_val,A_stiffner)

force_panel_1_R=np.sum(force_val[0:9])+np.sum(force_val[135:144])
force_panel_2_TR=np.sum(force_val[9:27])
force_panel_3_T=np.sum(force_val[27:45])
force_panel_2_TL=np.sum(force_val[45:63])
force_panel_1_L=np.sum(force_val[63:81])
force_panel_2_CL=np.sum(force_val[81:99])
force_panel_3_C=np.sum(force_val[99:117])
force_panel_2_CR=np.sum(force_val[117:135])

shear_force_1_R=np.sum(shear_val[0:9])+np.sum(shear_val[135:144])
shear_force_2_TR=np.sum(shear_val[9:27])
shear_force_3_T=np.sum(shear_val[27:45])
shear_force_2_TL=np.sum(shear_val[45:63])
shear_force_1_L=np.sum(shear_val[63:81])
shear_force_2_CL=np.sum(shear_val[81:99])
shear_force_3_C=np.sum(shear_val[99:117])
shear_force_2_CR=np.sum(shear_val[117:135])

print('shear force value on plate 1= ',shear_force_1_R)
print('Normal Force on plate 1=', force_panel_1_R)
print("shear force value on plate 3=", shear_force_3_C)
print("Normal Force", force_panel_3_C)


# total force on panel due to bending
Ftot= 53287.9338935 # N
shear_panel = 292524.656885 # N/m
ds=0.13089 #m
a=1 #m
b=2.35619 #m
EA_stiffner= 0.008645395725365318e9
EI_stiffner= 326.4536406876552
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
