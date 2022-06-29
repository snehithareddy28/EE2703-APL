import sys,os
#import numpy as np
#import scipy
from pylab import *
import scipy.linalg as s
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3

if len((sys.argv)>1) :
    Nx = int(sys.argv[1])
    Ny = int(sys.argv[2])
    radius = int(sys.argv[3])
    Niter = int(sys.argv[4])

else :
    Nx = 25; # size along x
    Ny = 25; # size along y
    radius = 8; # radius of central lead
    Niter = 1500; # number of iterations to perform

phi = np.zeros((Ny,Nx))

x = np.linspace(-0.5,0.5,Nx)
y = np.linspace(-0.5,0.5,Ny)
Y,X = meshgrid(y,x)

    #print(X)
    #print(Y)

radius = radius * (1/2) * ((1/ Nx ) + (1/ Ny ) )
# phi = np.zeros (( Nx , Ny ) , dtype = float )
ii = np.where ( X **2 + Y **2 <= radius **2)
phi [ ii ] = 1.0

#plot potential
plt.xlabel("X")
plt.ylabel("Y")
plt.contourf(X,Y,phi)
plt.colorbar()
plt.show()

