import numpy as np
import math
import matplotlib.pyplot as plt

E1 = 140E9  #Pa
E2 = 10E9   #Pa
G12 = 5E9   #Pa
v12 = 0.3

def getm(th):
    th = math.radians(th)
    return math.cos(th)

def getn(th):
    th = math.radians(th)
    return math.sin(th)
    

def getEx(th):
    m = getm(th)
    n = getn(th)
    recipEx = m**4/E1+(1/G12-2*v12/E1)*m**2*n**2+n**4/E2
    return 1/recipEx

def getEy(th):
    m = getm(th)
    n = getn(th)
    recipEy = n**4/E1+(1/G12-2*v12/E1)*m**2*n**2+m**4/E2
    return 1/recipEy

def getvxy(th, Ex):
    m = getm(th)
    n = getn(th)
    return Ex*(v12/E1*(m**4+n**4)-(1/E1+1/E2-1/G12)*m**2*n**2)

def getGxy(th):
    m = getm(th)
    n = getn(th)
    recipGxy = 2*(2/E1+2/E2+4*v12/E1-1/G12)*m**2*n**2+1/G12*(m**4+n**4)
    return 1/recipGxy

def getnxs(th, Ex):
    m = getm(th)
    n = getn(th)
    return Ex*((2/E1+2*v12/E1-1/G12)*m**3*n-(2/E2+2*v12/E1-1/G12)*n**3*m)

def getnys(th, Ey):
    m = getm(th)
    n = getn(th)
    return Ey*((2/E1+2*v12/E1-1/G12)*n**3*m-(2/E2+2*v12/E1-1/G12)*m**3*n)

Exlist = []
Eylist = []
Gxylist = []
vxylist = []
nxslist = []
nyslist = []

thlist = np.arange(0, 90.1, 0.1)

for th in thlist:
    Ex = getEx(th)
    vxy = getvxy(th, Ex)
    Ey = getEy(th)
    Gxy = getGxy(th)
    nxs = getnxs(th, Ex)
    nys = getnys(th, Ey)
    Exlist.append(Ex)
    vxylist.append(vxy)
    Eylist.append(Ey)
    Gxylist.append(Gxy)
    nxslist.append(nxs)
    nyslist.append(nys)


fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10, 8))
   
ax[0,0].plot(thlist, Exlist, label = "Ex")
ax[0,0].plot(thlist, Eylist, label = "Ey")
ax[1,0].plot(thlist, vxylist, label = "vxy")
ax[0,1].plot(thlist, Gxylist, label = "Gxy")
ax[1,1].plot(thlist, nxslist, label = "nxs")
ax[1,1].plot(thlist, nyslist, label = "nys")
ax[0,0].legend()
ax[0,1].legend()
ax[1,0].legend()
ax[1,1].legend()
plt.show()
    

        