import numpy as np
import matplotlib.pyplot as plt

# control runs
run0 = False

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

thickness_bending, weight_per_length_bending = calc_thickness_for_bending(M = M_x, R = R_outer, Yt = Xt, sf = sf_bending, rho = rho) 

print(f"Bending thickness in m = {thickness_bending}")
print(f"weight per length ={weight_per_length_bending}")


# Shear 
def calc_thickness_for_shear(V, S, R, sf, rho, theta = np.pi/2):
    qs = V/(np.pi*R)*np.sin(theta)
    tau_allow = S/sf
    t = qs/tau_allow
    weight_per_length = rho*2*np.pi*R*t
    return t,weight_per_length

thickness_shear, weight_per_length_shear = calc_thickness_for_shear(V = V_y, S = S, R = R_outer, sf = sf_shear, rho = rho, theta = np.pi/2)
print(f"Shear thickness in m = {thickness_shear}")
print(f"weight per length ={weight_per_length_shear}")



def calc_vonMises(sigma_z, tau_xy):
    return np.sqrt(sigma_z**2 + 3*tau_xy**2)

### this is incorrect, as it uses Ixx computed for a constant cross section
def calc_variable_thickness(V, M, n_points=3000):
    
    # Ixx = np.pi*R_outer**3*t
    # Iyy = Ixx
    
    y = np.linspace(0, R_outer, n_points)
    theta = np.zeros(len(y))
    t_arr = np.ones(len(y))*thickness_bending
    sigma_z = np.zeros(len(y))
    tau_xy = np.zeros(len(y))
    sigma_a = Xt/sf_global # Yt taken as minimum
    
    for i in range(len(y)):
        qs = V/(np.pi*R_outer)*np.sqrt(1-(y[i]/R_outer)**2)
        t_arr[i] =   1/sigma_a * np.sqrt((M*y[i]/(np.pi*R_outer**3))**2 + 3*qs**2)  
        if i > 0:
            theta[i] = np.arcsin(-(theta[0] - 1/R_outer*(y[i]-y[0])))
        # print(f'y = {y[i]} m, t = {t_arr[i]}m')
        #print(f'y = {y[i]} m, theta: {np.degrees(theta[i])}, t = {t_arr[i]}m')
        
    arc_length = theta*D_outer/2 # [m]
    mass_panel_segment = (arc_length[1:]-arc_length[:-1])*t_arr[1:]*rho # [kg/m]
    print(f'arc len: {arc_length}, mass_penl: {mass_panel_segment}')
    mass_unit_length_quarter_fuselage = np.sum(mass_panel_segment)
    mass_unit_length = 4*mass_unit_length_quarter_fuselage
    return y, t_arr, mass_unit_length

y_bending, t_bending, m_bending = calc_variable_thickness(V= 0, M= M_x)
y_shear, t_shear, m_shear = calc_variable_thickness(V= V_y, M = 0)
# y_combined, t_combined, m_combined = calc_variable_thickness(V= V_y, M = M_x)
# print(f'Mass of variable thickness panel: {m_combined}')

# # plt.step(y_bending, t_bending*10**3, where = 'pre', linestyle = '--', label = 'Bending Moment only')
# # plt.step(y_shear, t_shear*10**3,  where = 'pre', linestyle = '--',label = 'Shear force only')
# plt.step(y_combined, t_combined*10**3,  where = 'pre', label = 'Combined Loading')
# # plt.step(y_bending, t_bending*10**3, where = 'post',  label = 'Bending Moment only')
# # plt.step(y_shear, t_shear*10**3,  where = 'post',label = 'Shear force only')
# # plt.step(y_combined, t_combined*10**3,  where = 'post', label = 'Combined Loading')
# # plt.plot(y_bending, t_bending*10**3, marker = 'x', label = 'Bending Moment only')
# # plt.plot(y_shear, t_shear*10**3, marker = '*', label = 'Shear force only')
# # plt.plot(y_combined, t_combined*10**3, marker = 'o', label = 'Combined Loading')
# # plt.plot(y_bending, t_bending*10**3, linestyle = '--', label = 'Bending Moment only')
# # plt.plot(y_shear, t_shear*10**3, linestyle = '--',  label = 'Shear force only')
# # plt.plot(y_combined, t_combined*10**3, label = 'Combined Loading')
# plt.legend()
# plt.xlabel('Fuselage Height (y) [m]')
# plt.ylabel('Skin thickness (t) [mm]')
# plt.show()
    
    
# re-implementation in terms of theta (instead of y)
# correct for Ixx incorporating variable thickness
# solve for quarter fuselage (in first quadrant), assuming double symmetry

# thickness values are known a-prioro at midplane (t_min) and at top (t_max). Suppose that thickness has sinusuidal variation

# iteration 1 (initial)
def calc_linear_thicknes_var(t_min = thickness_shear, t_max = thickness_bending, isDiscrete = False, n_points = 4, interpolation_points = 4000):
    theta = np.linspace(0, np.pi/2, interpolation_points)
    y = R_outer*np.sin(theta)
    x = R_outer*np.cos(theta)
    
    print(f't_min: {t_min}; t_{max}: {t_max}')
    t = t_min + (t_max-t_min)*np.sin(theta)
    # for a sinosuidal thickness variation, the second moment of area is computed as a definite integral
    Ixx = 4*R_outer**3*(t_min * (np.pi/4) + (t_max-t_min)*(2/3))
    # print(f'Ixx: {Ixx}')
    
    sigma_z = M_x*y/Ixx
    qs = V_y/Ixx*t*(np.pi*R_outer)*np.cos(theta)
    
    sigma_vm = calc_vonMises(sigma_z= sigma_z, tau_xy= qs/t)
        
    arc_length = theta*D_outer/2 # [m]
    mass_panel_segment = (arc_length[1:]-arc_length[:-1])*t[1:]*rho # [kg/m]
    # print(f'theta: {np.degrees(theta)} arc len: {arc_length}, mass_penl: {mass_panel_segment}')
    mass_unit_length_quarter_fuselage = np.sum(mass_panel_segment)
    mass_unit_length = 4*mass_unit_length_quarter_fuselage
    SF = Xt/sigma_vm
    
    #### compute for a discrete case
    # theta_discrete = np.linspace(0, np.pi/2, n_points)
    theta_discrete = np.zeros(n_points-1)
    y_discrete = np.zeros(n_points-1)
    t_discrete = np.zeros(n_points-1)
    Ixx_discrete = np.zeros(n_points-1)
    theta_init = 0
    # interpolate (continuous to discrete)
    for i in range(n_points-1):
        theta_discrete[i] = theta[int(interpolation_points/n_points*(i+1))]
        y_discrete[i] = R_outer*np.sin(theta_discrete[i])
        # t_discrete[i] = t[int(interpolation_points/n_points*(i+1))]
        t_discrete[i] = t[np.argmin(np.abs(theta-theta_discrete[i]))]
        Ixx_discrete[i] = t_discrete[i]*R_outer**3*(0.5*(theta_discrete[i]-theta_init)-0.25*(np.sin(2*theta_discrete[i]) - np.sin(2*theta_init)))
        theta_init = theta_discrete[i]
    I_xx_discrete_tot = 4*np.sum(Ixx_discrete)
    
    
    # interpolate (discrete to continuous)
    sigma_z_discrete = M_x*y_discrete/I_xx_discrete_tot
    qs_discrete = V_y/I_xx_discrete_tot*t_discrete*(np.pi*R_outer)*np.cos(theta_discrete)
    # sigma_z_discrete = M_x*y/I_xx_discrete_tot
    # qs_discrete = V_y/I_xx_discrete_tot*t_discrete*(np.pi*R_outer)*np.cos(theta_discrete)
    
    sigma_vm_discrete = calc_vonMises(sigma_z= sigma_z_discrete, tau_xy= qs_discrete/t_discrete)
        
    arc_length_discrete = np.zeros(n_points)
    arc_length_discrete[1:] = theta_discrete*D_outer/2 # [m]
    mass_panel_segment_discrete = (arc_length_discrete[1:]-arc_length_discrete[:-1])*t_discrete*rho # [kg/m]
    # print(f'theta: {np.degrees(theta)} arc len: {arc_length}, mass_penl: {mass_panel_segment}')
    mass_unit_length_quarter_fuselage_discrete = np.sum(mass_panel_segment_discrete)
    mass_unit_length_discrete = 4*mass_unit_length_quarter_fuselage_discrete
    SF_discrete = Xt/sigma_vm_discrete
    
    if not isDiscrete:
        return y, t, SF, mass_unit_length
    elif isDiscrete:
        return y_discrete, t_discrete, SF_discrete, mass_unit_length_discrete
def calc_quadratic_thicknes_var(t_min = thickness_shear, t_max = thickness_bending, isDiscrete = False, n_points = 4, interpolation_points = 4000):
    print(f't_min: {t_min}; t_{max}: {t_max}')
    # alpha1 = 2
    # alpha2 = -2
    alpha1 = 1
    alpha2 = 0 #-0.5
    theta = np.linspace(0, np.pi/2, interpolation_points)
    y = R_outer*np.sin(theta)
    x = R_outer*np.cos(theta)
        
    t = t_min + alpha1*(t_max-t_min)*np.sin(theta) + alpha2*(t_max-t_min)*(np.sin(theta))**2
    #t = t_min + (t_max-t_min)*y/R_outer + alpha2*(t_max-t_min)*(y/R_outer)**2
    # plt.plot(y, t)
    # plt.show()
    # for a sinosuidal thickness variation, the second moment of area is computed as a definite integral
    Ixx = 4*R_outer**3*(t_min * (np.pi/4) + alpha1*(t_max-t_min)*(2/3) +  alpha2*(t_max-t_min)*(3/16))
    #print(f'Ixx: {Ixx}')
    
    sigma_z = M_x*y/Ixx
    qs = V_y/Ixx*t*(np.pi*R_outer)*np.cos(theta)
    
    sigma_vm = calc_vonMises(sigma_z= sigma_z, tau_xy= qs/t)
        
    arc_length = theta*D_outer/2 # [m]
    mass_panel_segment = (arc_length[1:]-arc_length[:-1])*t[1:]*rho # [kg/m]
    # print(f'theta: {np.degrees(theta)} arc len: {arc_length}, mass_penl: {mass_panel_segment}')
    mass_unit_length_quarter_fuselage = np.sum(mass_panel_segment)
    mass_unit_length = 4*mass_unit_length_quarter_fuselage
    SF = Xt/sigma_vm
    
    #### compute for a discrete case
    # theta_discrete = np.linspace(0, np.pi/2, n_points)
    theta_discrete = np.zeros(n_points-1)
    y_discrete = np.zeros(n_points-1)
    t_discrete = np.zeros(n_points-1)
    Ixx_discrete = np.zeros(n_points-1)
    theta_init = 0
    # interpolate (continuous to discrete)
    for i in range(n_points-1):
        theta_discrete[i] = theta[int(interpolation_points/n_points*(i+1))]
        y_discrete[i] = R_outer*np.sin(theta_discrete[i])
        # t_discrete[i] = t[int(interpolation_points/n_points*(i+1))]
        t_discrete[i] = t[np.argmin(np.abs(theta-theta_discrete[i]))]
        Ixx_discrete[i] = t_discrete[i]*R_outer**3*(0.5*(theta_discrete[i]-theta_init)-0.25*(np.sin(2*theta_discrete[i]) - np.sin(2*theta_init)))
        theta_init = theta_discrete[i]
    I_xx_discrete_tot = 4*np.sum(Ixx_discrete)
    
    
    # interpolate (discrete to continuous)
    sigma_z_discrete = M_x*y_discrete/I_xx_discrete_tot
    qs_discrete = V_y/I_xx_discrete_tot*t_discrete*(np.pi*R_outer)*np.cos(theta_discrete)
    # sigma_z_discrete = M_x*y/I_xx_discrete_tot
    # qs_discrete = V_y/I_xx_discrete_tot*t_discrete*(np.pi*R_outer)*np.cos(theta_discrete)
    
    sigma_vm_discrete = calc_vonMises(sigma_z= sigma_z_discrete, tau_xy= qs_discrete/t_discrete)
        
    arc_length_discrete = np.zeros(n_points)
    arc_length_discrete[1:] = theta_discrete*D_outer/2 # [m]
    mass_panel_segment_discrete = (arc_length_discrete[1:]-arc_length_discrete[:-1])*t_discrete*rho # [kg/m]
    # print(f'theta: {np.degrees(theta)} arc len: {arc_length}, mass_penl: {mass_panel_segment}')
    mass_unit_length_quarter_fuselage_discrete = np.sum(mass_panel_segment_discrete)
    mass_unit_length_discrete = 4*mass_unit_length_quarter_fuselage_discrete
    SF_discrete = Xt/sigma_vm_discrete
    
    if not isDiscrete:
        return y, t, SF, mass_unit_length
    elif isDiscrete:
        return y_discrete, t_discrete, SF_discrete, mass_unit_length_discrete


#### discrete thickness
def calc_discrete_thickness(t_arr, n_points = 4):
    pass
    
    
    
    


# print(f'Sigma_z: {sigma_z}, tau_xy = {qs/t}')
# print(f'SF: {Xt/sigma_vm}')
# print(f'thickness: {t}')
# print(f'weight per unit length: {mass_unit_length}')

# iteration 1
# conservative
y_var_1, t_var_1, SF_var_1, mass_unit_length_var_1 = calc_linear_thicknes_var(t_min = thickness_shear, t_max = thickness_bending)
#y_discrete_var_test, t_discrete_var_test, SF_discrete_var_test, mass_unit_length_discrete_var_test = calc_quadratic_thicknes_var(t_min = thickness_shear, t_max = thickness_bending, isDiscrete= True)

# iteration 2 
# reduce thickness at middle, and increase thickness at top 
delta_t_max_2 = 0.0002
delta_t_min_2 = -0.0002
#y_var_2, t_var_2, SF_var_2, mass_unit_length_var_2 = calc_linear_thicknes_var(t_min = thickness_shear + delta_t_min_2, t_max = thickness_bending + delta_t_max_2)
y_var_2, t_var_2, SF_var_2, mass_unit_length_var_2 = calc_quadratic_thicknes_var(t_min = thickness_shear + delta_t_min_2, t_max = thickness_bending + delta_t_max_2)
y_discrete_var_test, t_discrete_var_test, SF_discrete_var_test, mass_unit_length_discrete_var_test = calc_quadratic_thicknes_var(t_min = thickness_shear+ delta_t_min_2, t_max = thickness_bending+ delta_t_max_2, isDiscrete= True)
print(f'iteration 2: weight per unit length: {mass_unit_length_var_2}')

# iteration 3
# reduce thickness at middle, and increase thickness at top 
# quadratic 
delta_t_max_3 = 0.0005
delta_t_min_3 = -0.0004
y_var_3, t_var_3, SF_var_3, mass_unit_length_var_3 = calc_quadratic_thicknes_var(t_min = np.min(t_var_2) + delta_t_min_3, t_max = np.max(t_var_2) + delta_t_max_3)
print(f'iteration 3: weight per unit length: {mass_unit_length_var_3}')


# iteration 4
# discretised
y_discrete_var, t_discrete_var, SF_discrete_var, mass_unit_length_discrete_var = calc_quadratic_thicknes_var(t_min = np.min(t_var_2) + delta_t_min_3, t_max = np.max(t_var_2) + delta_t_max_3, isDiscrete= True)
print(f'Discrete case: thickness {t_discrete_var}')
print(f'iteration 4: weight per unit length: {mass_unit_length_discrete_var}')





# t_min = thickness_shear + delta_t_min
# t_max = thickness_bending + delta_t_max

# t = t_min + (t_max-t_min)*np.sin(theta)
# # for a sinosuidal thickness variation, the second moment of area is computed as a definite integral
# Ixx = 4*R_outer**3*(t_min * (np.pi/4) + (t_max-t_min)*(2/3))
# print(f'Ixx: {Ixx}')

# sigma_z = M_x*y/Ixx
# qs = V_y/Ixx*t*(np.pi*R_outer)*np.cos(theta)

# sigma_vm = calc_vonMises(sigma_z= sigma_z, tau_xy= qs/t)
# SF_arr = Xt/sigma_vm


# arc_length = theta*D_outer/2 # [m]
# mass_panel_segment = (arc_length[1:]-arc_length[:-1])*t[1:]*rho # [kg/m]
# #print(f'theta: {np.degrees(theta)} arc len: {arc_length}, mass_penl: {mass_panel_segment}')
# mass_unit_length_quarter_fuselage = np.sum(mass_panel_segment)
# mass_unit_length = 4*mass_unit_length_quarter_fuselage

# print(f'Sigma_z: {sigma_z}, tau_xy = {qs/t}')
# print(f'SF: {Xt/sigma_vm}')
# print(f'thickness: {t}')
# print(f'weight per unit length: {mass_unit_length}')


plt.figure('1')
plt.clf()
plt.plot(y_bending, t_bending*10**3, linestyle = '--', label = 'Bending Moment only')
plt.plot(y_shear, t_shear*10**3, linestyle = '--', label = 'Shear force only')
plt.plot(y_var_1, t_var_1*10**3, label = 'Combined Loading: Conservative')
plt.plot(y_var_2, t_var_2*10**3, label = 'Combined Loading: Linear')
#plt.plot(y_var_3, t_var_3*10**3,label = 'Combined Loading:  Quadratic')
plt.step(np.insert(y_discrete_var, 0, 0), np.insert(t_discrete_var, 0, t_discrete_var[0])*10**3,  color = 'black', marker = 'x',  label = 'Discrete')
plt.legend()
plt.xlabel('Fuselage Height (y) [m]')
plt.ylabel('Skin thickness (t) [mm]')


# plt.figure('1.1')
# plt.clf()
# plt.plot(y_var_3, t_var_3*10**3,label = 'Combined Loading:  Quadratic')
# plt.step(np.insert(y_discrete_var, 0, 0), np.insert(t_discrete_var, 0, t_discrete_var[0])*10**3,  color = 'black', marker = 'x',  label = 'Discrete')
# plt.legend()
# plt.xlabel('Fuselage Height (y) [m]')
# plt.ylabel('Skin thickness (t) [mm]')



plt.figure('2')
plt.clf()
plt.plot(y_var_1, SF_var_1, label = 'Actual: Conservative')
plt.plot(y_var_2, SF_var_2, label = 'Actual: Linear')
plt.plot(y_var_3, SF_var_3, label = 'Actual: Quadratic')
plt.plot(y_var_2, sf_global*np.ones(len(y_var_2)), linestyle = '--', dashes =(6,4), label = 'Constraint: Minimum')
plt.plot(y_discrete_var, sf_global*np.ones(len(y_discrete_var)), color = 'black', alpha = 0.5, marker = 'x', label = 'Discrete')
plt.legend()
plt.xlabel('Fuselage Height (y) [m]')
plt.ylabel('Skin thickness (t) [mm]')
plt.yticks(np.linspace(0, 3, 7))
plt.show()
plt.show()



