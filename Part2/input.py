import numpy as np

# Aluminum material properties 
Ex = 69*10**3 # [MPa]
Ey = 69*10**3 # [MPa]
Gxy = 26*10**3 # [MPa]
νxy = 0.29
Xt = 410 # [MPa]
Xc = 430 # [MPa]
Yt = 400 # [MPa]
Yc = 430 # [MPa]
S =	230 # [MPa]
Ρ  = 2770 # [kg/m^3]

# Coordinate system 
# x: pitch (+ve pointing outboard of left wing), 
# y: yaw (+ve pointing upwards), 
# z: roll (+ve pointing towards nose))

# Load case
M_x = 15000000 # [Nm]
V_y = 1500000 # [N]


