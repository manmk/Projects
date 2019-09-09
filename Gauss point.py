#Gauss point quadrature
#Mansour K

from __future__ import division
import numpy as np


def P1shapes(r,s):
    S = np.array([1-r-s,r,s])
    dSdr = np.array([-1,1,0])
    dSds = np.array([-1,0,1])
    return S,dSdr,dSds

def polyArea(x,y):
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

def a(x,y):
    return x + y

#Gauss point quadrature
qwgts=np.array([1/3,1/3,1/3])
rspts=np.array([[1/6,1/6],
[2/3,1/6],
[1/6,2/3]])

x = [-21.68467035,-21.17695462 ,-22.18700401]
y = [-9.94204652 ,-8.91258056 ,-9.12362242]

area,b,c = gradient(x,y)
xc = np.mean(x)
yc = np.mean(y)
a_centroid = a(xc,yc)
b = np.atleast_2d(b)
c = np.atleast_2d(c)
AK_simplified = (np.dot(b.T, b) + np.dot(c.T, c))*a_centroid*area

AK_quadrature = np.zeros((3,3))
for q in range(len(qwgts)):
    r = rspts[q,0]
    s = rspts[q,1]
    S,dSdx,dSdy,detJ = isopmap(x,y,r,s,P1shapes)
    x_physical = np.dot(x,S)
    y_physical = np.dot(y,S)
    a_gauss_point = a(x_physical,y_physical)

    dSdx =np.matrix(dSdx)
    dSdy =np.matrix(dSdy)
    wxarea = a_gauss_point*qwgts[q]*detJ/2
    AK_quadrature = AK_quadrature + (dSdx.T.dot(dSdx) + dSdy.T.dot(dSdy))* wxarea

# matrix
print AK_simplified
print " "
print AK_quadrature