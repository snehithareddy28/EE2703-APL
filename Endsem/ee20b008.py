"""
                                                END-SEMESTER EXAMINATION
                                                    ANDAPALLY SNEHITHA
                                                        EE20B008
                                            
"""
import cmath
import numpy as np
from pylab import *
import matplotlib.pyplot as plt

def getCurrents(l,N,Im,a):

    c = 2.9979e+8
    mu0 = np.pi*4e-7
    lamda = l*4
    f = c/lamda
    k = (2*np.pi)/lamda
    dz = l/N

    z = np.zeros(2*N+1)
    u = np.zeros(2*N-2)
    I = np.zeros(2*N+1)
    J = np.zeros(2*N-2)
    Iexp = np.zeros(2*N+1)

    i = np.arange(-N, N+1)
    z = i*dz
    i = np.arange(-N+1, N)
    u = i*dz
    u = np.delete(u, N-1)

    I[N] = Im
    for i in range(0,N):
        Iexp[i] = Im*np.sin(k*(l+z[i]))
    for i in range(N,2*N+1):
        Iexp[i] = Im*np.sin(k*(l-z[i]))

    #defining a Matrix which gives the Matrix M according to N
    M = np.identity(2*N-2, dtype = float)
    M = M/(2*np.pi*a)

    #computing vectors Rz and Ru
    Rz = np.sqrt(np.add(a*a, np.square(np.subtract.outer(z,z))))
    Ru = np.sqrt(np.add(a*a, np.square(np.subtract.outer(u,u))))
    print(Ru)

    list = np.arange(1,2*N)
    list = np.delete(list, N-1)
    
    #creating P and Pb matrix
    P = np.zeros((2*N-2, 2*N-2), dtype=np.complex128)
    for j in range(2*N-2):
        for i in range(2*N-2):
            P[i][j] = (mu0*cmath.exp(-1j*k*Ru[i][j])*dz)/(4*np.pi*Ru[i][j])

    Pb = np.zeros(2*N-2, dtype=np.complex128)

    c = 0
    for i in list:
        Pb[c] = (mu0*cmath.exp(-1j*k*Rz[i][N])*dz)/(4*np.pi*Rz[i][N])
        c+=1

    Q = np.zeros((2*N-2, 2*N-2), dtype=np.complex128)
    for j in range(2*N-2):
        for i in range(2*N-2):
            Q[i][j] = (P[i][j]*a*(1j*(k/Ru[i][j]) + 1/(Ru[i][j]*Ru[i][j])))/mu0

    Qb = np.zeros(2*N-2, dtype=np.complex128)
    c = 0
    for i in list:
        Qb[c] = (Pb[c]*a*(1j*(k/Rz[i][N]) + 1/(Rz[i][N]*Rz[i][N])))/mu0
        c+=1


    Y = Im*Qb
    X = M-Q
    X1= np.linalg.inv(X)
    J = np.matmul(X1,Y)
    #Inserting the boundary conditions
    I = np.insert(J,0,0)
    I = np.append(I,0)
    I = np.insert(I,N,complex(Im))


    #plotting both calculated and estimated graphs
    plt.plot(z,Iexp,'r', label = "Calculated current")
    plt.plot(z,np.abs(I),'g', label = "Theoretical current")
    plt.xlabel(r"z",size=14)
    plt.ylabel(r"Current",size=14)
    plt.title("Plots of calculated and theoretical current values")
    plt.legend(loc = "upper right")
    plt.grid(True)
    plt.show()

getCurrents(0.5,100,1,0.01)