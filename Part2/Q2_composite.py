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
Yc = 00*10**6 # [Pa]
S =	100*10**6 # [Pa]
rho  = 1610 # [kg/m^3]
t_ply = 0.135*10**-3 # [m]

# Coordinate system 
# x: pitch (+ve pointing outboard of left wing), 
# y: yaw (+ve pointing upwards), 
# z: roll (+ve pointing towards nose))

# Load case
M_x = 15000000 # [Nm]
V_y = 1500000 # [N]

# Geometry
D_outer = 6 # [m]
t = 0.002 # [m] TODO: placeholder
D_inner = D_outer - 2*t # [m]
R_outer = D_outer/2 # [m]
R_inner = D_inner/2 # [m]

# Knockdown/safety factors 
# TODO: find values from slides/literature
pass


# Black Aluminum design
# quasi isotropic laminate
quasi_isotropic_layup = [90, +45, -45, 0, 0, -45, +45, 90] # total thickness = 1.08 [mm], too thin
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




