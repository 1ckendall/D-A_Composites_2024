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
        self.v21 = self.v12 * self.E1 / self.E2  # Minor Poisson
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
        self.Q11 = self.E1**2 / (self.E1 - v12**2 * self.E2)
        self.Q12 = self.v12 * self.E1 * self.E2 / (self.E1 - v12**2 * self.E2)
        self.Q22 = self.E1 * self.E2 / (self.E1 - v12**2 * self.E2)
        self.Q66 = self.G12
        self.Qmat = np.array(
            [[self.Q11, self.Q12, 0], [self.Q12, self.Q22, 0], [0, 0, self.Q66]]
        )

    def getQbarmat(self):
        self.Qxx = self.Q11 * self.m**4 + 2 * (
            self.Q12 + 2 * self.Q66 * self.m**2 * self.n**2
        )
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
    def __init__(self, plys, Nx=0, Ny=0, Ns=0, Mx=0, My=0, Ms=0):
        self.plys = plys
        self.Nx = Nx
        self.Ny = Ny
        self.Ns = Ns
        self.Mx = Mx
        self.My = My
        self.Ms = Ms
        self.getABD()
        self.abd = np.linalg.inv(self.ABD)

    def getABD(self):
        self.A = np.zeros([3, 3])
        self.B = np.zeros([3, 3])
        self.D = np.zeros([3, 3])

        z = 0 # Define the bottom of the laminate as z = 0
        for ply in self.plys:
            self.A += ply.Qbarmat * ((z + ply.t) - z)
            self.B += 1 / 2 * ply.Qbarmat * ((z + ply.t) ** 2 - z**2)
            self.D += 1 / 3 * ply.Qbarmat * ((z + ply.t) ** 3 - z**3)
            z += ply.t
        AB = np.concatenate((self.A, self.B), axis=0)
        BD = np.concatenate((self.B, self.D), axis=0)
        self.ABD = np.concatenate((AB, BD), axis=1)

if __name__ == "__main__":
    E1 = 140e9
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
