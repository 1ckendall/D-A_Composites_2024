from Class import Lamina, Laminate
from bucklingforpanels import stiffened_panel_buckling


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
Panel_3=[45, -45, 0, 0, 90, 90, 0, 0, -45, 45]
Panel_3_lamina = [Lamina(i, Ex, Ey, Gxy, vxy, Xt, Xc, Yc, Yt, S, t_ply) for i in Panel_3]
T_skin=t_ply*len(Panel_3)
# Creating Laminate
laminate_3=Laminate(Panel_3_lamina)

# calling A and D matrices
A_mat=laminate_3.A
D_mat=laminate_3.D

# total force on panel due to bending
Ftot= -1303685.336 # N/
shear_panel = 6186.8399 # N/m
ds=0.13089 #m
a=1 #m
b=2.35619 #m
EA_stiffner= 0
EI_stiffner=0




#
bucklepanel_3=stiffened_panel_buckling(Ftot,shear_panel,ds,a,b,EA,EI,A_mat,D_mat,T_skin)
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