# Mansour K
# A sample code in Python that solves the Poisson equation in 1D using the FEM
from numpy import *

# Problem definition
x0=0.0
xL=15.0
Nx=101
fi0=3.0   # Dirichlet condition
qL=13.0   # Neumann condition
Q0=5.0    # Heat load
km=1.0    # material

long=xL-x0
nnodes=5
nnod=2
nelems=nnodes-1
interv=long/(nnodes-1)
xi=zeros(nnodes,float)
xi[0]=x0
for i in range(1,nnodes):
    xi[i]=xi[i-1]+interv

f=zeros((nelems,nnod),float)
K=zeros((nelems,nnod,nnod),float)

for i in range(0,nelems):
    for j in range(0,nnod):
        f[i,j]=Q0*interv/2
        for k in range(0,nnod):
            K[i,j,k]=pow(-1,j+k)*(km/interv)
        

fg=zeros(nnodes,float)
Kg=zeros((nnodes,nnodes),float)
for i in range(0,nelems):
    for j in range(0,nnod):
        fg[i+j]=fg[i+j]+f[i,j]
        for k in range(0,nnod):
            Kg[i+j,i+k]=Kg[i+j,i+k]+K[i,j,k]
           
fg[nnodes-1]=fg[nnodes-1]-qL


Kf=Kg[1:,1:]
cfi=Kg[1:,0]*fi0
ff=fg[1:]-cfi

fi=inner(linalg.inv(Kf),ff)

fit=hstack((fi0,fi))
 realFi=-(Q0/(2*km))*pow(xi,2)+(-qL+Q0*(xL-x0))*xi/km+fi0;

print





