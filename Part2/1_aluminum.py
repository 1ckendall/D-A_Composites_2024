import numpy as np
import matplotlib.pyplot as plt

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
# I_xx = np.pi * (D_outer**4 - D_inner**4) / 64 # [m]
# I_yy = I_xx # [m] # assume 
sf_global = 1.5
sf_bending = 1.5 # safety factor source: Airframe Stress Analysis and Sizing
sf_shear = 1.5 # safety factor source: Airframe Stress Analysis and Sizing


# Uniaxial Bending 
def calc_thickness_for_bending(M,R,Yt,sf,rho): # Momnent, Outer Radius, yeild strength, safety factor, rho; all values in SI units (Nm, m, Pa, Kg/m)
  # thin-walled assumption made, do note that tensile is limiting
  # justification for using Yt: sheet metal aluminum. Grain is aligned (Xt) along the circumferential/hoop direction of the fuselage, thus Yt (transverse grain) is along longitudinal direction 
  t =  M/(np.pi*(R**2)*(Yt/sf)) #M/(np.pi*(R*2)*(Yt/sf))
  weight_per_length = rho*2*np.pi*R*t
  return t,weight_per_length

thickness_bending, wt = calc_thickness_for_bending(M = M_x, R = R_outer, Yt = Xt, sf = sf_bending, rho = rho) 

print(f"Bending thickness in m = {thickness_bending}")
print(f"weight per length ={wt}")


# Shear 
def calc_thickness_for_shear(V, S, R, sf, rho, theta = np.pi/2):
    qs = V/(np.pi*R)*np.sin(theta)
    tau_allow = S/sf
    t = qs/tau_allow
    weight_per_length = rho*2*np.pi*R*t
    return t,weight_per_length

thickness_shear, wt = calc_thickness_for_shear(V = V_y, S = S, R = R_outer, sf = sf_shear, rho = rho, theta = np.pi/2)
print(f"Shear thickness in m = {thickness_shear}")
print(f"weight per length ={wt}")



def calc_vonMises(sigma_z, tau_xy):
    return np.sqrt(sigma_z**2 + 3*tau_xy**2)

def calc_variable_thickness(V, M, n_points=4):
    
    Ixx = np.pi*R_outer**3*t
    Iyy = Ixx
    
    y = np.linspace(0, R_outer, n_points)
    t_arr = np.ones(len(y))*thickness_bending
    sigma_z = np.zeros(len(y))
    tau_xy = np.zeros(len(y))
    sigma_a = Xt/sf_global # Yt taken as minimum
    
    for i in range(len(y)):
        qs = V/(np.pi*R_outer)*np.sqrt(1-(y[i]/R_outer)**2)
        t_arr[i] =   1/sigma_a * np.sqrt((M*y[i]/(np.pi*R_outer**3))**2 + 3*qs**2)  
        print(f'y = {y[i]} m, t = {t_arr[i]}m')
        
    return y, t_arr

y_bending, t_bending = calc_variable_thickness(V= 0, M= M_x)
y_shear, t_shear = calc_variable_thickness(V= V_y, M = 0)
y_combined, t_combined = calc_variable_thickness(V= V_y, M = M_x)

plt.step(y_bending, t_bending*10**3, where = 'pre', linestyle = '--', label = 'Bending Moment only')
plt.step(y_shear, t_shear*10**3,  where = 'pre', linestyle = '--',label = 'Shear force only')
plt.step(y_combined, t_combined*10**3,  where = 'pre', label = 'Combined Loading')
# plt.step(y_bending, t_bending*10**3, where = 'post',  label = 'Bending Moment only')
# plt.step(y_shear, t_shear*10**3,  where = 'post',label = 'Shear force only')
# plt.step(y_combined, t_combined*10**3,  where = 'post', label = 'Combined Loading')
# plt.plot(y_bending, t_bending*10**3, marker = 'x', label = 'Bending Moment only')
# plt.plot(y_shear, t_shear*10**3, marker = '*', label = 'Shear force only')
# plt.plot(y_combined, t_combined*10**3, marker = 'o', label = 'Combined Loading')
# plt.plot(y_bending, t_bending*10**3, linestyle = '--', label = 'Bending Moment only')
# plt.plot(y_shear, t_shear*10**3, linestyle = '--',  label = 'Shear force only')
# plt.plot(y_combined, t_combined*10**3, label = 'Combined Loading')
plt.legend()
plt.xlabel('Fuselage Height (y) [m]')
plt.ylabel('Skin thickness (t) [mm]')
plt.show()
    
    

# Bucklings









