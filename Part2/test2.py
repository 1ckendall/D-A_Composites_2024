from Class import Lamina, Laminate
#from bucklingforpanels import stiffened_panel_buckling
from Test_discritization import compute_forces,create_circle
import numpy as np

# testing



y=create_circle(145,3)
# skin parameters
t_ply=0.135e-3
Panel1 = t_ply* 80
Panel2 = t_ply* 70
Panel3 = t_ply* 60
Eskin_val=[(45.23*10**9),(60.22*10**9),(71.21*10**9)]
T_skin=[Panel1,Panel2,Panel3]

#stiffner parameters
A_stiffner=[0,0.0001,0.00015,0.0002]
Estiffner_val=[0, 48.216e9, 60.265e9, 74.07e9]


force_val,shear_val,weight,sig_stiffner, sig_z=compute_forces(T_skin,Eskin_val,A_stiffner,Estiffner_val)
#print(shear_val)
print('weight=',weight)
'''
force_panel_1_R=(np.max(force_val[0:9])*9)+(np.sum(force_val[135:144])*9)
force_panel_2_TR=np.max(np.abs(force_val[9:27]))*18
force_panel_3_T=np.max(np.abs(force_val[27:45]))*18
force_panel_2_TL=np.max(np.abs(force_val[45:63]))*18
force_panel_1_L=np.max(np.abs(force_val[63:81]))*18
force_panel_2_CL=np.max(np.abs(force_val[81:99]))*18
force_panel_3_C=np.max(np.abs(force_val[99:117]))*18
force_panel_2_CR=np.max(np.abs(force_val[117:135]))*18
'''

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

#print(sig_boom)
stress_panel_3T_stiffners=sig_stiffner[27:45]
stress_panel_3C_stiffners=sig_stiffner[99:117]

stress_panel_2T_stiffners=sig_stiffner[9:27]
stress_panel_2C_stiffners=sig_stiffner[81:99]

stress_panel_1top_stiffners=sig_stiffner[0:9]
stress_panel_1bot_stiffners=sig_stiffner[135:144]
print(sig_z)
print(stress_panel_2T_stiffners)
print(stress_panel_3C_stiffners)