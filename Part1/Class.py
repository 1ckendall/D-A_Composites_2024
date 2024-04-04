import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg as la
from pprint import pprint
from warnings import warn

class Lamina:
    def __init__(self, th, E1, E2, G12, v12, Xt, Xc, Yt, Yc, S, t,
                 Rfpt=None, Rfpc=None, Rtt=None, Rtc=None, Rls=None, Rts=None,
                 sigma_1 = None, sigma_2 = None, sigma_3 = None, tau_23 = None, tau_31 = None, tau_21 = None,
                 m_sigmaf = 1.3):
        
        # MATERIAL PROPERTIES
        self.th = np.radians(th)  # angle of fibres (degrees,
        # converted to radians for internal maths)
        self.m = np.cos(self.th)
        self.n = np.sin(self.th)
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
        
        self.failuremode = None

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
        ## version 1: Transformation Matrix (used as verification check)
        self.T = np.matrix([[self.m**2, self.n**2, 2*self.m*self.n],
                            [self.n**2, self.m**2, -2*self.m*self.n],
                            [-self.m*self.n, self.m*self.n, self.m**2-self.n**2]])
        self.Tinv = np.linalg.inv(self.T)
        
        assert np.allclose(self.Tinv, np.matrix([[self.m**2, self.n**2, -2*self.m*self.n],
                            [self.n**2, self.m**2, 2*self.m*self.n],
                            [self.m*self.n, -self.m*self.n, self.m**2-self.n**2]]))
        
        Qmat = self.Qmat
        Qmat[:,2] = 2*Qmat[:,2]
        
        Qbarmat = self.Tinv @ Qmat @ self.T
        
        Qbarmat = Qbarmat
        Qbarmat[:,2] -= 0.5*Qbarmat[:,2]
        
        ## version 2 (source: Engineering Mechanics of Composite Materials)
        self.Qxx = self.m**4 * self.Q11 + self.n**4 * self.Q22  + 2*self.m**2*self.n**2*self.Q12 + 4*self.m**2*self.n**2*self.Q66
        self.Qyy = self.n**4 * self.Q11 + self.m**4 * self.Q22 + 2*self.m**2*self.n**2*self.Q12 + 4*self.m**2*self.n**2*self.Q66
        self.Qxy = self.m**2 * self.n**2 * self.Q11 + self.m**2 * self.n**2 * self.Q22 + (self.m**4 + self.n**4)*self.Q12 - 4*self.m**2*self.n**2 *self.Q66
        self.Qxs = self.m **3 * self.n * self.Q11 - self.m * self.n **3 * self.Q22 - self.m *self.n *(self.m**2 - self.n**2) * self.Q12 - 2*self.m*self.n * (self.m**2 - self.n**2) * self.Q66
        self.Qys = self.n **3 * self.m * self.Q11 - self.n * self.m **3 * self.Q22 + self.m *self.n *(self.m**2 - self.n**2) * self.Q12 + 2*self.m*self.n * (self.m**2 - self.n**2) * self.Q66
        self.Qss = self.m**2 * self.n**2 * self.Q11 + self.m**2 * self.n**2 * self.Q22 - 2*self.m**2 *self.n**2 * self.Q12 + (self.m**2-self.n**2)**2 * self.Q66
        
        self.Qbarmat  = np.matrix([[self.Qxx, self.Qxy,  self.Qxs], 
                       [self.Qxy,  self.Qyy,  self.Qys],
                       [self.Qxs,  self.Qys,  self.Qss]])
        
        assert np.allclose(self.Qbarmat, Qbarmat)
        
    def maxStressFibreFail(self):
        if self.sigma_1 / self.Xt >= 1:
            self.failuremode = "Tensile Fibre"
        elif self.sigma_1 / -self.Xc >= 1:
            self.failuremode = "Compressive Fibre"
    
    def maxStressInterFibreFail(self):
        if self.sigma_2 / self.Yt >= 1: #or self.sigma_3 / self.Rtt >= 1:
            self.failuremode = "Tensile Inter-Fibre"
        elif self.sigma_2 / -self.Yc >= 1 :#or self.sigma_3 / self.Rtc >=1
            self.failuremode = "Compressive Inter-Fibre"
            
    def maxStressShearFail(self):
        if self.tau_21 / self.S >= 1 : #or self.tau_31 / self.Rls >= 1:
            self.failuremode = "Shear Parallel to Fibres"
        elif self.tau_21 / -self.S >=1:  
            self.failuremode = "Shear Parallel to Fibres" 
        #elif self.tau_23 / self.Rts >= 1: 
           # self.failuremode = "Shear Perpendicular to Fibres"
        
    def PuckFibreFail(self, sigma_2, gamma_21, epsilon1T, epsilon1C, epsilon1, m_sigmaf = 1.3):
        """Calculates the Puck failure of the fibres in the UD lamina

        Args:
            sigma_2 (_type_): Stress across fibres
            gamma_21 (_type_): Shear strain of the unidirectional layer
            epsilon1T (_type_): Tensile failure strain of a unidirectional layer in x1 direction
            epsilon1C (_type_): Compression failure strain of a unidirectional layer in x1 direction
            epsilon1 (_type_): Strain in x1 direction
            m_sigmaf (float, optional): _description_. Defaults to 1.3.
        """
        if sigma_2 < 0:
            failurecriterion = 1/epsilon1T*(epsilon1+self.v12/self.E1*m_sigmaf*sigma_2)
            if failurecriterion >= 1:
                self.failuremode = "FFT"
                #print(f"Ply failed at {sigma_2}: Tensile Fibre Failure")
        elif sigma_2 > 0:
            failurecriterion = 1/epsilon1C*(abs(epsilon1+self.v12/self.E1*m_sigmaf*sigma_2)) + 10*(gamma_21)**2
            if failurecriterion >= 1:
                self.failuremode = "FFC"
                #print(f"Ply failed at {sigma_2}: Compressive Fibre Failure (Kinking)")
    
    def get_RAperpperp(self, S21, p_perppara_minus, Yc):
        """RA⊥⊥: Fracture resistance of the action plane against its fracture due to transverse/transverse shear stressing"""
        return S21/2*p_perppara_minus*(np.sqrt(1+2*p_perppara_minus*Yc/S21)-1)
    
    def get_p_perpperp_minus(self, p_perppara_minus, RAperpperp, S21):
        return p_perppara_minus*RAperpperp/S21
    
    def get_tau_21c(self, S21, p_perpperp_minus):
        return S21*np.sqrt(1-2*p_perpperp_minus)
        
    def PuckIFF(self, sigma_2, sigma_1, sigma_1D, tau_21, tau_21c, S21, RAperpperp, p_perppara_plus, p_perppara_minus, p_perpperp_minus, Y_T, Y_C):
        """Calculates the Puck inter-fibre failure

        Args:
            sigma_2 (_type_): Stress across fibres
            sigma_1 (_type_): Stress along fibres
            sigma_1D (_type_): Stress value for linear degradation(sigma1D>0 for sigma1>0; sigma1D<0 for sigma1<0)
            tau_21 (_type_): Shear stresses of a unidirectional layer
            tau_21c (_type_): Shear stress at the `turning point' of the (sigma2, tau21) fracture curve
            S21 (_type_): Shear strength of a unidirectional layer transverse and parallel to the fibre direction
            RAperpperp (_type_): Fracture resistance of the action plane against its fracture due to transverse/transverse shear stressing
            p_perppara_plus (_type_): Slope of the sigma_n-tau_n1 fracture envelope for sigma_n geq 0 at sigma_n=0
            p_perppara_minus (_type_): Slope of the sigma_n-tau_n1 fracture envelope for sigma_n leq 0 at sigma_n=0
            p_perpperp_minus (_type_): Slope of the sigma_n-tau_nt fracture envelope for sigma_n leq 0 at sigma_n=0
            Y_T (_type_): Tensile strength of the unidirectional layer transverse to the fibre direction
            Y_C (_type_): Compressive strength of the unidirectional layer transverse to the fibre direction
        """
        if sigma_2 >= 0:
            # IFF A
            failurecriterion = np.sqrt((tau_21/S21)**2+(1+p_perppara_plus*Y_T/S21)**2*(sigma_2/Y_T)**2)+p_perppara_plus*sigma_2/S21+abs(sigma_1/sigma_1D)
            if failurecriterion >= 1:
                self.failuremode = "IFF A"
        elif sigma_2 < 0 and 0<=abs(sigma_2/tau_21)<=RAperpperp/abs(tau_21c):
            # IFF B
            failurecriterion = 1/S21*(np.sqrt(tau_21**2+(p_perppara_minus*sigma_2)**2)+p_perppara_minus*sigma_2)+abs(sigma_1/sigma_1D)
            if failurecriterion >= 1:
                self.failuremode = "IFF B"
        elif sigma_2 < 0 and 0<=abs(tau_21/sigma_2<=abs(tau_21c)/RAperpperp):
            # IFF C
            failurecriterion = ((tau_21/(2*(1+p_perpperp_minus)*S21))**2+(sigma_2/Y_C)**2)*Y_C/(-sigma_2)+abs(sigma_1/sigma_1D)
            if failurecriterion >=1:
                self.failuremode = "IFF C"
        else:
            # Something is funky
            warn("Something is wrong in IFF calculation")

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
        self.abd = la.inv(self.ABD)
        self.getEngineeringConst()

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
            
        for i in range(self.n_plys):
            self.A += self.plys[i].Qbarmat * (self.z[i+1] - self.z[i])
            self.B += 1 / 2 * self.plys[i].Qbarmat * ((self.z[i+1]) ** 2 - self.z[i]**2)
            self.D += 1 / 3 * self.plys[i].Qbarmat * ((self.z[i+1]) ** 3 - self.z[i]**3)
        AB = np.concatenate((self.A, self.B), axis=0)
        BD = np.concatenate((np.transpose(self.B), self.D), axis=0)
        self.ABD = np.concatenate((AB, BD), axis=1)
        
    def getEngineeringConst(self):
        self.Axx = self.A[0,0]
        self.Ayy = self.A[1,1]
        self.Ass = self.A[2,2]
        self.Axy = self.A[0,1]
        
        self.dxx = self.abd[3,3]
        self.dyy = self.abd[4,4]
        self.dss = self.abd[5,5]
        self.dxy = self.abd[3,4]
        
        self.A_const = self.Axx*self.Ayy - self.Axy*self.Axy
        self.h = np.max(self.z)-np.min(self.z)
        
        # in-plane engineering constants
        self.Ex = self.A_const/(self.h*self.Ayy)
        self.Ey = self.A_const/(self.h*self.Axx)
        self.Gxy = self.Ass/self.h
        self.vxy = self.Axy/self.Ayy
        self.vyx = self.Axy/self.Axx
        
        # flexural engineering constants
        self.Exb = 12/(self.h**3 * self.dxx)
        self.Eyb = 12/(self.h**3 * self.dyy)
        self.Gxyb = 12/(self.h**3 * self.dss)
        self.vxyb = -self.dxy/self.dxx
        self.vyxb = -self.dxy/self.dyy
        
    # def getStressStrain(self):
    #     # compute global strain (laminate level) on mid-plane
    #     self.load = np.array([self.Nx, self.Ny, self.Ns, self.Mx, self.My, self.Ms]).T
    #     self.strainMidplane = self.abd @ self.load # [e^0_x, e^0_y, e^0_s, kx, ky, ks]
    #     self.globalstrainVector = np.zeros((3*self.n_plys))
    #     self.globalstressVector = np.zeros((3*self.n_plys))
    #     self.localstrainVector = np.zeros((3*self.n_plys))
    #     self.localstressVector = np.zeros((3*self.n_plys))
        
    #     # compute loacal strain (lamina level)
    #     for i in range(self.n_plys):
    #         # global coordinate system (x,y)
    #         self.globalstrainVector[i:i+3] = self.strainMidplane[:3] + 0.5*(self.z[i+1]-self.z[i])*self.strainMidplane[3:] # [ex, ey, es], computed on the midplane of each lamina
    #         self.globalstressVector[i:i+3]  = self.plys[i].Qbarmat @ self.globalstrainVector[i:i+3]
            
    #         # local/principal coordinate system (1,2)
    #         self.localstrainVector[i:i+3]  = self.plys[i].T @ self.globalstrainVector[i:i+3]
    #         self.localstressVector[i:i+3]  = self.plys[i].T @ self.globalstressVector[i:i+3]
            
    #     return(self.localstrainVector, self.localstressVector)
    
    
    def getStressStrain(self):
        # compute global strain (laminate level) on mid-plane
        self.load = np.array([self.Nx, self.Ny, self.Ns, self.Mx, self.My, self.Ms]) #.reshape(-1,1)
        # self.strainMidplane = self.abd @ self.load # [e^0_x, e^0_y, e^0_s, kx, ky, ks]
        self.strainMidplane = np.linalg.solve(self.ABD.astype('float64'), self.load.astype('float64'))
        self.globalstrainVector = np.zeros((3*self.n_plys))
        self.globalstressVector = np.zeros((3*self.n_plys))
        self.localstrainVector = np.zeros((3*self.n_plys))
        self.localstressVector = np.zeros((3*self.n_plys))
        self.z_lamina_midplane = np.zeros(self.n_plys)
      
        # compute loacal strain (lamina level)
        for i in range(self.n_plys):
            self.z_lamina_midplane[i] = 0.5*(self.z[i+1]+self.z[i]) # z-values from laminate datum (ie: laminate midplane) to laminas' datums (ie: lamina midplane)
            
            # global coordinate system (x,y)
            self.globalstrainVector[3*i:3*(i+1)] = self.strainMidplane[:3] + self.z_lamina_midplane[i]*self.strainMidplane[3:] # [ex, ey, es], computed on the midplane of each lamina
            self.globalstressVector[3*i:3*(i+1)]  = self.plys[i].Qbarmat @ self.globalstrainVector[3*i:3*(i+1)]
            
            # local/principal coordinate system (1,2)
            self.localstrainVector[3*i:3*(i+1)]  = self.plys[i].T @ self.globalstrainVector[3*i:3*(i+1)]
            self.localstressVector[3*i:3*(i+1)]  = self.plys[i].T @ self.globalstressVector[3*i:3*(i+1)]
            
        self.e11 = self.localstrainVector[0::3]
        self.e22 = self.localstrainVector[1::3]
        self.e12 = self.localstrainVector[2::3]
        
        self.sigma11 = self.localstressVector[0::3]
        self.sigma22 = self.localstressVector[1::3]
        self.sigma12 = self.localstressVector[2::3]

        