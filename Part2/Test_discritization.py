import numpy as np
import matplotlib.pyplot as plt
from Class import Lamina, Laminate
from scipy.interpolate import CubicSpline
'''
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

'''
#--------------------------------------------------------------Begin here-----------------------------------------------
# Geometry
D_outer = 6 # [m]
t = 0.002 # [m] TODO: placeholde
D_inner = D_outer - 2*t # [m]
R_outer = D_outer/2 # [m]
R_inner = D_inner/2 # [m]
Ixx_t=np.pi*(R_outer**3) # Ixx / thickness(unknown)



t_ply = 0.135*10**-3 # [m]

# Coordinate system 
# x: pitch (+ve pointing outboard of left wing), 
# y: yaw (+ve pointing upwards), 
# z: roll (+ve pointing towards nose)

#loading parameters
M_x = 15000000 # [Nm]
V_y = 1500000 # [N]

# ----------------------------------------------------------------------------------Functions --------------------------------------------------------
def calc_loadvector(V, M, n_points=4):  
    y = np.linspace(0, R_outer, n_points)    
    sigma_z_times_t_arr = np.zeros(len(y))
    qs_arr = np.zeros(len(y))
    for i in range(len(y)):
        sigma_z_times_t_arr[i] = M*y[i]/(np.pi*R_outer**3)
        qs_arr[i] = V/(np.pi*R_outer)*np.sqrt(1-(y[i]/R_outer)**2)
    return sigma_z_times_t_arr, qs_arr,y

def create_circle(num_points,R):
    # Angle range from 0 to 2*pi evenly spaced by the number of points
    Side_boom = np.linspace(np.pi/(num_points-1), 2*np.pi + np.pi/(num_points-1), num_points)
    Bending_boom = np.linspace(0, 2*np.pi, num_points)
    # Calculate x and y coordinates for each angle
    X = R*np.cos(Bending_boom)
    Y = R*np.sin(Bending_boom)

    x = R*np.cos(Side_boom)
    y = R*np.sin(Side_boom)

    arc_length = 2*np.pi*R/(num_points-1) #np.sqrt((x[i] - x[i-1])**2 + (y[i] - y[i-1])**2)

    interp_points = 1000  # Number of points after interpolation
    spline = CubicSpline(np.arange(num_points), np.column_stack([x, y]), axis=0)
    interp_indices = np.linspace(0, num_points - 1, interp_points)
    interp_circle = spline(interp_indices)

    # Extract interpolated x and y coordinates
    interp_x = interp_circle[:, 0]
    interp_y = interp_circle[:, 1]


    plt.figure(figsize=(8,8))
    plt.plot(x,y,'bo-')#,
    plt.axis('equal')
    plt.grid()
    plt.plot(interp_x, interp_y)  # Plot the circle with blue dots connected by lines
    
    return x, y, arc_length, X,Y


#------------------------------------------------------------ Initializing----------------------------------------------
def compute_forces(T_skin,Eskin_val,A_stiffner,Estiffner_val):

   # Geometry
    D_outer = 6 # [m]
    t = 0.002 # [m] TODO: placeholde
    D_inner = D_outer - 2*t # [m]
    R_outer = D_outer/2 # [m]
    R_inner = D_inner/2 # [m]
    Ixx_t=np.pi*(R_outer**3) # Ixx / thickness(unknown)
    
    num_points=144
    x,y,arc_length,X,Y=create_circle(num_points+1,3)
    x=x[:-1]
    y=y[:-1]
    X=X[:-1]
    Y=Y[:-1]
    #print(arc_length)
    #--------------------------------------------------------------Begin here-----------------------------------------------
 

    # Panel node breakdown
    Area_panel= [arc_length*t for t in T_skin]
    Panel_Areas=[]
    Panel_Areas.extend([Area_panel[0]]*9)
    Panel_Areas.extend([Area_panel[1]]*18)
    Panel_Areas.extend([Area_panel[2]]*18)
    Panel_Areas.extend([Area_panel[1]]*18)
    Panel_Areas.extend([Area_panel[0]]*18)
    Panel_Areas.extend([Area_panel[1]]*18)
    Panel_Areas.extend([Area_panel[2]]*18)
    Panel_Areas.extend([Area_panel[1]]*18)
    Panel_Areas.extend([Area_panel[0]]*9)
    # thickness array
    T_skinarr=[]
    T_skinarr.extend([T_skin[0]]*9)
    T_skinarr.extend([T_skin[1]]*18)
    T_skinarr.extend([T_skin[2]]*18)
    T_skinarr.extend([T_skin[1]]*18)
    T_skinarr.extend([T_skin[0]]*18)
    T_skinarr.extend([T_skin[1]]*18)
    T_skinarr.extend([T_skin[2]]*18)
    T_skinarr.extend([T_skin[1]]*18)
    T_skinarr.extend([T_skin[0]]*9)

    # Elasticity Array
    E_skin=[]
    E_skin.extend([Eskin_val[0]]*9)
    E_skin.extend([Eskin_val[1]]*18)
    E_skin.extend([Eskin_val[2]]*18)
    E_skin.extend([Eskin_val[1]]*18)
    E_skin.extend([Eskin_val[0]]*18)
    E_skin.extend([Eskin_val[1]]*18)
    E_skin.extend([Eskin_val[2]]*18)
    E_skin.extend([Eskin_val[1]]*18)
    E_skin.extend([Eskin_val[0]]*9)

    # Stiffner Area Addition
    stiff_loc=np.zeros(num_points)
    for i in range(num_points):
        if i % 2 == 0:  # Check if index is even
            stiff_loc[i]=A_stiffner[0]
        else:
            if i <= 9:
                stiff_loc[i]=A_stiffner[1]
            elif i <= 27:
                stiff_loc[i]=A_stiffner[2]
            elif i <= 45:
                stiff_loc[i]=A_stiffner[3]
            elif i <= 63:
                stiff_loc[i]=A_stiffner[2]
            elif i <= 81:
                stiff_loc[i]=A_stiffner[1]
            elif i <= 99:
                stiff_loc[i]=A_stiffner[2]
            elif i <= 117:
                stiff_loc[i]=A_stiffner[3]
            elif i <= 135:
                stiff_loc[i]=A_stiffner[2]
            elif i <= 144:
                stiff_loc[i]=A_stiffner[1]
    #print(stiff_loc)
    '''
    stiff_loc.extend([A_stiffner[1]]*8)
    stiff_loc.extend([A_stiffner[0]])
    stiff_loc.extend([A_stiffner[3]]*17)
    stiff_loc.extend([A_stiffner[0]])
    stiff_loc.extend([A_stiffner[2]]*17)
    stiff_loc.extend([A_stiffner[0]])
    stiff_loc.extend([A_stiffner[3]]*17)
    stiff_loc.extend([A_stiffner[0]])
    stiff_loc.extend([A_stiffner[2]]*17)
    stiff_loc.extend([A_stiffner[0]])
    stiff_loc.extend([A_stiffner[1]]*17)
    stiff_loc.extend([A_stiffner[0]])
    stiff_loc.extend([A_stiffner[3]]*17)
    stiff_loc.extend([A_stiffner[0]])
    stiff_loc.extend([A_stiffner[2]]*17)
    stiff_loc.extend([A_stiffner[0]])
    stiff_loc.extend([A_stiffner[1]]*9)
    '''
    # Elasticity matrices
    E_stiffner=np.zeros(num_points)
    for i in range(num_points):
        if stiff_loc[i]== A_stiffner[1]:
            E_stiffner[i]=Estiffner_val[1]
        elif stiff_loc[i]==A_stiffner[3]:
            E_stiffner[i]=Estiffner_val[3]
        elif stiff_loc[i]==A_stiffner[2]:
             E_stiffner[i]=Estiffner_val[2]
    '''
    E_stiffner.extend([Estiffner_val[1]]*8)
    E_stiffner.extend([Estiffner_val[0]])
    E_stiffner.extend([Estiffner_val[3]]*17)
    E_stiffner.extend([Estiffner_val[0]])
    E_stiffner.extend([Estiffner_val[2]]*17)
    E_stiffner.extend([Estiffner_val[0]])
    E_stiffner.extend([Estiffner_val[3]]*17)
    E_stiffner.extend([Estiffner_val[0]])
    E_stiffner.extend([Estiffner_val[1]]*17)
    E_stiffner.extend([Estiffner_val[0]])
    E_stiffner.extend([Estiffner_val[3]]*17)
    E_stiffner.extend([Estiffner_val[0]])
    E_stiffner.extend([Estiffner_val[2]]*17)
    E_stiffner.extend([Estiffner_val[0]])
    E_stiffner.extend([Estiffner_val[3]]*17)
    E_stiffner.extend([Estiffner_val[0]])
    E_stiffner.extend([Estiffner_val[1]]*9)
    '''
    #--------------------------------------------------------- setting discritization basis----------------------------------
    # Discretized boom parameters Stiffness and Area
    #creating boom at X,Y coordinate for discretized panel
    # calculating Boom Area
    B_bending=[]
    B_areas=[]
    E_booms=[]
    for i in range(num_points):
        sig_s1 = ((M_x*y[i-1])/Ixx_t)
        sig_s2 = ((M_x*y[i-num_points+1])/Ixx_t)
        sig_boom = ((M_x*y[i])/Ixx_t)

        # Boom Area = Stringer Area + skin 1 contribution + skin 2 contribution
        Sai = stiff_loc[i]
        #print(Sai)
        B1i = (Panel_Areas[i-1]/6)*(2+(sig_s1/sig_boom))
        B2i = (Panel_Areas[i]/6)*(2+(sig_s2/sig_boom))
        B_area = Sai+B1i+B2i
        Ixx_boom = (y[i]**2)*B_area
        B_bending.append(Ixx_boom)
        B_areas.append(B_area)
        E_boom=((E_skin[i]*Panel_Areas[i])+(E_skin[i-1]*Panel_Areas[i-1])+(E_stiffner[i]*stiff_loc[i]))/(stiff_loc[i]+Panel_Areas[i]+Panel_Areas[i-1])
        E_booms.append(E_boom)
        #print(E_boom)
        
    # calculate stress due to moment and shear on the discritised length
    sig_z=[]
    delQ=[]
    Ixx_full=np.sum(B_bending)

    for i in range(num_points):
        sig_boom=M_x*y[i]/Ixx_full
        delQi=-(V_y/Ixx_full)*B_areas[i]*y[i]
        sig_z.append(sig_boom)
        delQ.append(delQi)
    #print(sig_z)

    shear_val=np.zeros(num_points)
    sum=0
    for i in range(1,num_points):
        sum += delQ[i-int((3/4)*num_points)]
        #print(sum)
        shear_val[i-int((3/4)*num_points)+1]=sum
    #print(shear_val)
    
    sig_skin=[]
    skin1contri=[]
    skin2contri=[]
    sig_stiffner=[]
    for i in range(num_points):
        sig_boom=sig_z[i]
        sig_skin1 = (E_skin[i-1]/E_booms[i])*sig_boom
        sig_skin2 = (E_skin[i]/E_booms[i])*sig_boom
        sig_stiff = (E_stiffner[i]/E_booms[i])*sig_boom
        sig_stiffner.append(sig_stiff)
        skin1contri.append(sig_skin1)
        skin2contri.append(sig_skin2)
    
    for i in range(num_points):
        sigmaval=skin1contri[i]+((skin1contri[i]-skin2contri[i-1])/2)
        sig_skin.append(sigmaval)

    force_val=np.zeros(num_points)
    for i in range(num_points):
        force_val[i]= sig_z[i]*B_areas[i]

    
    #print(arc_length)
    node=-int(num_points/4)
    #print(sig_z[node])
    #print(shear_val[node])
    #print(sig_skin[node])
    #print(sig_stiffner[node])

    weight=(np.sum(Panel_Areas)*1610)+(np.sum(stiff_loc)*1610)

    plt.scatter(x[0:9],y[0:9],color='blue',marker='s',s=33)
    plt.scatter(x[9:27],y[9:27],color='red',marker='s',s=66)
    plt.scatter(x[27:45],y[27:45],color='black',marker='s',s=100)
    plt.scatter(x[45:63],y[45:63],color='red',marker='s',s=66)
    plt.scatter(x[63:81],y[63:81],color='blue',marker='s',s=33)
    plt.scatter(x[81:99],y[81:99],color='red',marker='s',s=66)
    plt.scatter(x[99:117],y[99:117],color='black',marker='s',s=100)
    plt.scatter(x[117:135],y[117:135],color='red',marker='s',s=66)
    plt.scatter(x[135:144],y[135:144],color='blue',marker='s',s=33)
    '''
    plt.scatter(x[node],y[node],color='black',marker='o',s=150)
    #plt.scatter(x[node-num_points+1],y[node-num_points+1],color='blue',marker='o',s=150)
    plt.scatter(x[node-1],y[node-1],color='red',marker='o',s=150)
    plt.scatter(X[node],Y[node],color='green',marker='x')
    '''
    plt.show()

    return force_val, shear_val, weight, sig_stiffner, sig_z

#-------------------------------------------------------------skin parameters-------------------------------------------

# skin parameters
t_ply=0.135e-3
Eskin_val=[(45.23*10**9),(60.22*10**9),(71.21*10**9)]
T_skin=[11*t_ply,13*t_ply,15*t_ply]

#stiffner parameters
A_stiffner=[0,6.823e-5,0.000230,0.000160]
Estiffner_val=[0, 48.216e9, 71.265e9, 74.07e9]

'''
force_val,shear_val, weight=compute_forces(T_skin,Eskin_val,A_stiffner)
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

print(shear_force_1_R)
print(force_panel_1_R)
print(shear_force_3_T)
print(force_panel_3_T)
'''