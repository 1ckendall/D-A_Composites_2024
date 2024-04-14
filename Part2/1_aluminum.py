import numpy as np

# Aluminum material properties 
Ex = 69*10**9 # [Pa]
Ey = 69*10**9 # [Pa]
Gxy = 26*10**9 # [Pa]
νxy = 0.29
Xt = 410*10**6 # [Pa]
Xc = 430*10**6 # [Pa]
Yt = 400*10**6 # [Pa]
Yc = 430*10**6 # [Pa]
S =	230*10**6 # [Pa]
rho  = 2770 # [kg/m^3]

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

# area moment of inertia (circular; thin walled assumption NOT made)
I_xx = np.pi * (D_outer**4 - D_inner**4) / 64 # [m]
I_yy = I_xx # [m] # assume 

sf_bending = 1.5 # safety factor source: Airframe Stress Analysis and Sizing
sf_shear = 1.5 # safety factor source: Airframe Stress Analysis and Sizing


# Uniaxial Bending 
def calc_thickness_for_bending(M,R,Yt,sf,rho): # Momnent, Outer Radius, yeild strength, safety factor, rho; all values in SI units (Nm, m, Pa, Kg/m)
  # thin-walled assumption made, do note that tensile is limiting
  # justification for using Yt: sheet metal aluminum. Grain is aligned (Xt) along the circumferential/hoop direction of the fuselage, thus Yt (transverse grain) is along longitudinal direction 
  t =  M/(np.pi*(R**2)*(Yt/sf)) #M/(np.pi*(R*2)*(Yt/sf))
  weight_per_length = rho*2*np.pi*R*t
  return t,weight_per_length

thickness, wt = calc_thickness_for_bending(M = M_x, R = R_outer, Yt = Xt, sf = sf_bending, rho = rho) 

print(f"Bending thickness in m = {thickness}")
print(f"weight per length ={wt}")


# Shear 
def calc_thickness_for_shear(V, S, R, sf, rho, theta = np.pi/2):
    qs = V/(np.pi)*R*np.sin(theta)
    tau_allow = S/sf
    t = qs/tau_allow
    weight_per_length = rho*2*np.pi*R*t
    return t,weight_per_length

thickness, wt = calc_thickness_for_shear(V = V_y, S = S, R = R_outer, sf = sf_shear, rho = rho, theta = np.pi/2)
print(f"Shear thickness in m = {thickness}")
print(f"weight per length ={wt}")


# Bucklings









