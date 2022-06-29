import numpy as np
import math
import matplotlib.pyplot as py
import scipy.integrate as integrate

# Defining functions using numpy

def make_periodic(f_old):
    def f_new(x):
        return f_old(np.remainder(x,2*np.pi))
    return f_new

def expo(x):
    return np.exp(x)

expo_periodic = make_periodic(expo)

def coscos(x) :
    return np.cos(np.cos(x))

# Creating a range of x values from -2π to 4π.
x = np.linspace(-2*math.pi,4*math.pi,10000)

# Plotting functions
# $e^x$ in semi-log scale along with expected Fourier $e^x$}
py.semilogy(x, expo(x),color='r')
py.semilogy(x, expo_periodic(x),color='b')
py.xlabel('x',fontsize=15)
py.ylabel('$e^x$',fontsize=15)
py.grid()

l=[]
l.append('Actual value')
l.append('Periodic Extension')
py.legend(l)
py.show()

# Plot of cos(cos(x)) in linear scale
py.plot(x, coscos(x), color='b')
py.xlabel('x',fontsize=15)
py.ylabel('cos(cos(x))',fontsize=15)
py.grid()
py.show()

#this function returns the first n fourier coefficients of the function f using the integration method
def u(x,k):
    return expo(x)*np.cos(k*x)
def v(x,k):
    return expo(x)*np.sin(k*x)
def p(x,k):
    return coscos(x)*np.cos(k*x)
def q(x,k):
    return coscos(x)*np.sin(k*x)

# For expo(x)
expo_a0= (integrate.quad(u,0,2*np.pi,args=(0,))[0])/(2*np.pi)
expo_a = np.array([integrate.quad(u,0,2*np.pi,args=(i,))[0] for i in range(1,26)])/(np.pi)
expo_b = np.array([integrate.quad(v,0,2*np.pi,args=(i,))[0] for i in range(1,26)])/(np.pi)

# For cos(cos (x))
coscos_a0= (integrate.quad(p,0,2*np.pi,args=(0,))[0])/(2*np.pi)
coscos_a = np.array([integrate.quad(p,0,2*np.pi,args=(i,))[0] for i in range(1,26)])/(np.pi)
coscos_b = np.array([integrate.quad(q,0,2*np.pi,args=(i,))[0] for i in range(1,26)])/(np.pi)

    #print(expo_a0)
    #print(expo_a)
    #print(expo_b)

F=[None]*(1+len(expo_a)+len(expo_b))
F[0] = expo_a0
F[1::2] = expo_a
F[2::2] = expo_b
F = np.asarray(F)

Fmatrix = np.c_[F]
print(Fmatrix)

# Fourier Coefficients of $e^x$ by direct integration in semilog scale
py.semilogy(abs(F),'o',color='r',markersize=3)
py.xlabel('n',fontsize=15)
py.ylabel('Magnitude of Fourier co-efficients (in log scale)',fontsize=15)
py.grid()
py.show()

# Fourier Coefficients of $e^x$ by direct integration in loglog scale
py.loglog(abs(F),'o',color='r',markersize=3)
py.xlabel('n',fontsize=15)
py.ylabel('Magnitude of Fourier co-efficients (in log scale)',fontsize=15)
py.grid()
py.show()


G=[None]*(1+len(coscos_a)+len(coscos_b))
G[0] = coscos_a0
G[1::2] = coscos_a
G[2::2] = coscos_b
G = np.asarray(G)

Gmatrix = np.c_[G]
print(Gmatrix)

#Fourier Coefficients of cos(cos(x)) by direct integration in semilog scale
py.semilogy(abs(G),'o',color='r',markersize=3)
py.xlabel('n',fontsize=15)
py.ylabel('Magnitude of Fourier co-efficients',fontsize=15)
py.grid()
py.show()

#Fourier Coefficients of cos(cos(x)) by direct integration in loglog scale
py.loglog(abs(G),'o',color='r',markersize=3)
py.xlabel('n',fontsize=15)
py.ylabel('Magnitude of Fourier co-efficients',fontsize=15)
py.grid()
py.show()

X = np.linspace(0,2*np.pi,401)
X=X[:-1]
b1 = expo(X)
b2 = coscos(X)
A = np.zeros((400,51))
A[:,0] = 1

for k in range(1,26):
    A[:,2*k-1] = np.cos(k*X)
    A[:,2*k] = np.sin(k*X)

c1 = np.linalg.lstsq(A,b1,rcond=-1)[0]
c2 = np.linalg.lstsq(A,b2,rcond=-1)[0]

#Fourier co-efficients of $e^x$ in semilog scale
py.semilogy(abs(c1),'o', color='g', markersize = 3, label = 'using lstsq method')
py.semilogy(abs(F),'o', color='r', markersize = 3, label = 'using direct intergration method')
py.xlabel('n',fontsize=15)
py.ylabel('Fourier co-efficients',fontsize=15)
py.legend([" Using lstsq method ", " Using direct intergration method "])
py.grid()
py.show()

#Fourier co-efficients of $e^x$ in loglog scale
py.loglog(abs(c1),'o', color='g', markersize = 3, label = 'using lstsq method')
py.loglog(abs(F),'o', color='r', markersize = 3, label = 'using direct intergration method')
py.xlabel('n',fontsize=15)
py.ylabel('Fourier co-efficients',fontsize=15)
py.legend([" Using lstsq method ", " Using direct intergration method "])
py.grid()
py.show()

#Fourier co-efficients of cos(cos(x)) in semilog scale
py.semilogy(abs(c2),'o', color='g', markersize = 3, label = 'using lstsq method')
py.semilogy(abs(G),'o', color='r', markersize = 3, label = 'using direct intergration method')
py.xlabel('n',fontsize=15)
py.ylabel('Fourier co-efficients',fontsize=15)
py.legend([" Using lstsq method ", " Using direct intergration method "])
py.grid()
py.show()

#Fourier co-efficients of cos(cos(x)) in loglog scale
py.loglog(abs(c2),'o', color='g', markersize = 3, label = 'using lstsq method')
py.loglog(abs(G),'o', color='r', markersize = 3, label = 'using direct intergration method')
py.xlabel('n',fontsize=15)
py.ylabel('Fourier co-efficients',fontsize=15)
py.legend([" Using lstsq method ", " Using direct intergration method "])
py.grid()
py.show()

f_final = A.dot(c1)
#Original function $e^x$ and the Reconstructed function from Fourier Coefficients
py.semilogy(X,f_final,'o',color = 'r',markersize=0.5)
py.semilogy(X,expo(X),'o',color ='g',markersize=0.5)
py.xlabel('x',fontsize=15)
py.ylabel('$e^x$',fontsize=15)
py.legend(['Reconstructed from Fourier Coefficients', 'Original Function'])
py.grid()
py.show()

g_final = A.dot(c2)
#Original function cos(cos(x)) and the Reconstructed function from Fourier Coefficients
py.plot(X,g_final,'-',color = 'r',markersize=1)
py.plot(X,coscos(X),'o',color ='g',markersize=1)
py.xlabel('x',fontsize=15)
py.ylabel('cos(cos(x))',fontsize=15)
py.legend(['Reconstructed from Fourier Coefficients', 'Original Function'])
py.grid()
py.show()


print (" Absolute difference between Fourier Coefficients in expo (x) = exp (x) is {}". format ( np .sum (abs ( c1 - F ) ) ) )
print (" Largest Deviation between Fourier Coefficients in expo(x) = exp(x) is {}". format ( np . amax (abs( c1 - F ) ) ) )
print (" Mean Error in expo (x) = exp(x) is {}". format ( np . mean (abs( c1 - F ) ) ) )

print (" Absolute difference between Fourier Coefficients in coscos (x) = cos (cos (x)) is {}". format ( np .sum (abs ( c2 - G ) ) ) )
print (" Largest Deviation between Fourier Coefficients in coscos (x) = cos (cos (x)) is {}". format ( np . amax (abs ( c2 - G ) ) ) )
print (" Mean Error in coscos (x) = cos(cos(x)) is {}". format ( np . mean (abs( c2 - G ) ) ) )