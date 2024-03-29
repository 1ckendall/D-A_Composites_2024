import matplotlib.pyplot as plt
import numpy as np
from pprint import pprint

class Lamina:
    def __init__(self, th, E1, E2, G12, v12, Xt, Xc, Yt, Yc, S, t,
                 Rfpt=None, Rfpc=None, Rtt=None, Rtc=None, Rls=None, Rts=None,
                 sigma_1 = None, sigma_2 = None, sigma_3 = None, tau_23 = None, tau_31 = None, tau_21 = None,
                 m_sigmaf = 1.3):
        
        # MATERIAL PROPERTIES
        self.th = np.radians(th)  # angle of fibres (degrees,
        # converted to radians for internal maths)
        self.m = np.cos(th)
        self.n = np.sin(th)
        self.E1 = E1  # E along fibres
        self.E2 = E2  # E across fibres
        self.G12 = G12  # Shear modulus
        self.v12 = v12  # Major Poisson
        self.v21 = self.v12 * self.E2 / self.E1  # Minor Poisson
        self.Xt = Xt  # Strength along fibres (tensile)
        self.Xc = Xc  # Strength along fibres (compression)
        self.Yt = Yt  # Strength across fibres (tensile)
        self.Yc = Yc  # Strength across fibres (compression)
        self.S = S  # Shear strength
        self.t = t  # Thickness
        self.Rfpt = Rfpt # Fibre parallel tensile strength
        self.Rfpc = Rfpc # Fibre parallel compressive strength
        self.Rtt = Rtt # Transverse tensile strength
        self.Rtc = Rtc # Transverse compressive strength
        self.Rls = Rls # Longitudinal shear strength
        self.Rts = Rts # Transverse shear strength
        
        # LOADING
        self.sigma_1 = sigma_1
        self.sigma_2 = sigma_2
        self.sigma_3 = sigma_3
        self.tau_23 = tau_23
        self.tau_31 = tau_31
        self.tau_21 = tau_21
        
        # Puck Parameters
        self.m_sigmaf = m_sigmaf 
        
        self.getQmat()
        self.getQbarmat()
        
        self.failed = False
        self.failuremode = None # 1 -> Tensile fibre failure, 2 -> Compressive fibre failure

    def getQmat(self):
        self.Q = 1-self.v12*self.v21
        self.Q11 = self.E1 / self.Q
        self.Q12 = self.v12 * self.E2 / self.Q 
        self.Q22 = self.E2 / self.Q
        self.Q66 = self.G12
        self.Qmat = np.array(
            [[self.Q11, self.Q12, 0], [self.Q12, self.Q22, 0], [0, 0, self.Q66]]
        )

    def getQbarmat(self):
        self.Qxx = self.Q11 * self.m**4 + 2 * (self.Q12 + self.Q66) * self.m**2 * self.n**2 + self.Q22 * self.n**4
        self.Qxy = (
            self.Q11 + self.Q22 - 4 * self.Q66
        ) * self.m**2 * self.n**2 + self.Q12 * (self.m**4 + self.n**4)
        self.Qyy = (
            self.Q11 * self.n**4
            + 2 * (self.Q12 + 2 * self.Q66) * self.m**2 * self.n**2
            + self.Q22 * self.m**4
        )
        self.Qxs = (self.Q11 - self.Q12 - 2 * self.Q66) * self.n * self.m**3 + (
            self.Q12 - self.Q22 + 2 * self.Q66
        ) * self.n**3 * self.m
        self.Qys = (self.Q11 - self.Q12 - 2 * self.Q66) * self.m * self.n**3 + (
            self.Q12 - self.Q22 + 2 * self.Q66
        ) * self.m**3 * self.n
        self.Qss = (
            self.Q11 + self.Q22 - 2 * self.Q12 - 2 * self.Q66
        ) * self.m**2 * self.n**2 + self.Q66 * (self.n**4 + self.m**4)
        self.Qbarmat = np.array(
            [
                [self.Qxx, self.Qxy, self.Qxs],
                [self.Qxy, self.Qyy, self.Qys],
                [self.Qxs, self.Qys, self.Qss],
            ]
        )
        
    def maxStressFibreFail(self):
        if self.sigma_1 > 0:
            if self.sigma_1 / self.Rfpt >= 1:
                self.failed = True
                self.failuremode = "Tensile Fibre Failure"
            if self.sigma_1 / -self.Rfpc >= 1:
                self.failed = True
                self.failuremode = "Compressive Fibre Failure"
        
    def PuckFibreFail(self):
        pass
        
    def PuckIFF(self):
        pass


class Laminate:
    def __init__(self, plys, Nx=0, Ny=0, Ns=0, Mx=0, My=0, Ms=0, midplane = True):
        self.plys = plys
        self.Nx = Nx
        self.Ny = Ny
        self.Ns = Ns
        self.Mx = Mx
        self.My = My
        self.Ms = Ms
        self.midplane = midplane
        self.getABD()
        self.abd = np.linalg.inv(self.ABD)

    def getABD(self):
        self.n_plys = len(self.plys)
        self.A = np.zeros([3, 3])
        self.B = np.zeros([3, 3])
        self.D = np.zeros([3, 3])
        self.z = np.zeros(self.n_plys+1)
        
        # define the z-position of laminae. Datum = bottom 
        for i in range(self.n_plys):
            self.z[i+1] = (i+1)*self.plys[i].t 
                    
        # change datum of z-position from bottom to midplane
        if self.midplane: 
            self.z -= np.max(self.z)/2
        
        # print(f'Z-arr: {self.z}')
    
        for i in range(self.n_plys):
            self.A += self.plys[i].Qbarmat * (self.z[i+1] - self.z[i])
            #print(f'Z-diff: {(self.z[i+1] - self.z[i])}, Qbarmat: {self.plys[i].Qbarmat}, A-contribution: { self.plys[i].Qbarmat * (self.z[i+1] - self.z[i])}')
            self.B += 1 / 2 * self.plys[i].Qbarmat * ((self.z[i+1]) ** 2 - self.z[i]**2)
            self.D += 1 / 3 * self.plys[i].Qbarmat * ((self.z[i+1]) ** 3 - self.z[i]**3)
        AB = np.concatenate((self.A, self.B), axis=0)
        BD = np.concatenate((self.B, self.D), axis=0)
        self.ABD = np.concatenate((AB, BD), axis=1)

if __name__ == "__main__":
    E1 = 140
    E2 = 10
    G12 = 5
    v12 = 0.3
    t = 0.125
    Xt = 1500
    Xc = 1200
    Yt = 50
    Yc = 250
    S = 70
    anglelist = [0, 45, -45, 90, 90, -45, 45, 0]
    plylist = []
    for angle in anglelist:
        plylist.append(Lamina(angle, E1, E2, G12, v12, Xt, Xc, Yt, Yc, S, t))
    L = Laminate(plylist)
    pprint(L.ABD)
    


#assignment part 
 #lamina properties for overall assignment DO NOT change this properties 
E1 = 145.3e3 #[MPa]
E2 = 8.5e3 #[MPa]
v12 = 0.31 #[--]
G12 = 4.58e3 #[MPa]
Xt = 1932 #[MPa]
Yt = 108 #[MPa]
Xc = 1480 #[MPa]
Yc = 220 #[MPa]
S = 132.8 #[MPa]
t= 0.125 #[mm] #this one is assumption can be changed 
#question 1 Daniel 
    #first defining angles for theta in question 

thetaquestion1  = list(range(0,1,1)) #input angle in the laminate definition reason for 0 to 90 range is that the properties after 
                                                                                               #90 would be the same as before 90 
nquestion1 = list(range(1,2,1)) # n can be only integers 
Exmatrix = np.zeros((len(nquestion1),len(thetaquestion1)))
Eymatrix = np.zeros((len(nquestion1),len(thetaquestion1)))
Vxymatrix= np.zeros((len(nquestion1),len(thetaquestion1)))
Vyxmatrix = np.zeros((len(nquestion1),len(thetaquestion1)))
Gxymatrix = np.zeros((len(nquestion1),len(thetaquestion1)))
#calculating the laminate properties for question 1 [15/+theta1/-theta1/75/75/75/75/-theta1/+theta/15] 
plylist= [] 
for j in range(len(nquestion1)):

    for i in range(len(thetaquestion1)) :

        layup  = [15,+thetaquestion1[i],-thetaquestion1[i],75,75,75,75,-thetaquestion1[i],+thetaquestion1[i],15]
        layup *= nquestion1[j]
        # print(layup)

        for angle in layup:
            plylist.append(Lamina(angle,E1,E2,G12,v12,Xt,Xc,Yt,Yc,S,t))
        Result = Laminate(plylist)
        h = t* len(layup)
        
        Amatrix = Result.A
        print(f'A_matrix: {Amatrix}')
        EX = (Amatrix[0][0] * Amatrix[1][1] - Amatrix[0][1]**2) / (h* Amatrix[1][1]) 
        EY = (Amatrix[0][0] * Amatrix[1][1] - Amatrix[0][1]**2) / (h* Amatrix[0][0]) 
        VXY  = Amatrix[0][1] / Amatrix[1][1]
        VYX  = Amatrix[0][1] / Amatrix[0][0]
        GXY = Amatrix[2][2] / h 
        
        print(f'Theta: {thetaquestion1[i]}, Ex: {EX}, Ey: {EY} \n ')
        
        Exmatrix[j][i] = EX
        Eymatrix [j][i] = EY
        Vxymatrix [j][i] = VXY
        Vyxmatrix[j][i] = VYX  
        Gxymatrix [j][i] = GXY  

print(Exmatrix)

# Create subplots
fig, ax = plt.subplots(figsize=(10, 8))

# Plot properties for each n value
for j in range(len(nquestion1)):
    # Plot EX for the current n value
    ax.plot(thetaquestion1, Exmatrix[j], label=f'n = {nquestion1[j]}')

# Set title and legend
ax.set_title('Laminate Properties for Different n Values')
ax.legend()

# Set labels
ax.set_xlabel('Angle (degrees)')
ax.set_ylabel('EX')

# Show plot
plt.show()
fig, ax = plt.subplots(figsize=(10, 8))

# Plot properties for each n value
for j in range(len(nquestion1)):
    # Plot EX for the current n value
    ax.plot(thetaquestion1, Eymatrix[j], label=f'n = {nquestion1[j]}')

# Set title and legend
ax.set_title('Laminate Properties for Different n Values')
ax.legend()

# # Set labels
# ax.set_xlabel('Angle (degrees)')
# ax.set_ylabel('Ey')

# # Show plot
# plt.show()
# fig, ax = plt.subplots(figsize=(10, 8))
# # Plot properties for each n value
# for j in range(len(nquestion1)):
#     # Plot EX for the current n value
#     ax.plot(thetaquestion1, Vxymatrix[j], label=f'n = {nquestion1[j]}')

# # Set title and legend
# ax.set_title('Laminate Properties for Different n Values')
# ax.legend()

# # Set labels
# ax.set_xlabel('Angle (degrees)')
# ax.set_ylabel('Vxy')

# # Show plot
# plt.show()
# fig, ax = plt.subplots(figsize=(10, 8))
# # Plot properties for each n value
# for j in range(len(nquestion1)):
#     # Plot EX for the current n value
#     ax.plot(thetaquestion1, Vyxmatrix[j], label=f'n = {nquestion1[j]}')

# # Set title and legend
# ax.set_title('Laminate Properties for Different n Values')
# ax.legend()

# # Set labels
# ax.set_xlabel('Angle (degrees)')
# ax.set_ylabel('vyx')

# # Show plot
# plt.show()
# fig, ax = plt.subplots(figsize=(10, 8))
# # Plot properties for each n value
# for j in range(len(nquestion1)):
#     # Plot EX for the current n value
#     ax.plot(thetaquestion1, Gxymatrix[j], label=f'n = {nquestion1[j]}')

# # Set title and legend
# ax.set_title('Laminate Properties for Different n Values')
# ax.legend()

# # Set labels
# ax.set_xlabel('Angle (degrees)')
# ax.set_ylabel('Gxy')

# # Show plot
# plt.show()







