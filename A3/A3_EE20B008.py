from pickle import TRUE
import numpy as np
import scipy.special as sc
import scipy.linalg as sl
import matplotlib.pyplot as py

######### q2 #########
try :
    f = np.loadtxt("fitting.dat",dtype="float")
    # Loading the data to a numpy array
except OSError:
    print("Error: Couldn't open file")
    exit()

sigma = np.logspace(-1,-3,9)
# Sigma array
#print(sigma)


######### q4 #########
def g(t,A,B):
    return A*sc.jn(2,t)+B*t
# Function to generate f(t)    


######### q3 #########
py.plot(f[:,0],f[:,1:])
py.plot(f[:,0],g(f[:,0],1.05,-0.105),color='k')

# Plotting the data from fitting.dat

r=[]
'''for x in range(len(sigma)):
    r.append(sigma[x])'''
for x in range(len(sigma)):
    r.append(f"$\sigma${x+1}={sigma[x]: .3f}")
r.append("True Value")

py.legend(r)
py.xlabel("Time")
py.ylabel("f(t)+Noise")

py.show()


######### q5 #########
py.errorbar(f[::5,0],f[::5,1],sigma[0],fmt='ro')
py.plot(f[:,0],g(f[:,0],1.05,-0.105))

# Plotting with error bars
p=[]
p.append("Error bar")
p.append("True value")

py.legend(p)
py.show()


######### q6 #########
# Forming the matrix M
j = sc.jn(2,f[:,0])
'''
ls=[]
for t in range(len(f[:,0])):
    ls.append(sc.jn(2,f[t][0]))
j = np.array(ls)
'''
M = np.c_[j,f[:,0]]

a = np.array([1.05,-0.105])
N = np.c_[a]
s = np.matmul(M,N)
o = np.array(g(f[:,0],1.05,-0.105))
np.array_equal(o,s)
if TRUE:
    print("Yes! Arrays are equal.")
else:
    print("Arrays are not equal :(")


######### q7 #########
A1 = np.linspace(0,2,21)
B1 = np.linspace(-0.2,0,21)
# Creating the error matrix
E = np.zeros([21,21])
for j in range(len(A1)):
    for k in range(len(B1)):
        for l in range(len(f[:,0])):
            E[j][k] += (f[l][1]- g(f[l][0], A1[j], B1[k]))**2/101
#print(E)


######### q8 #########
# Plotting the contour
cnt=py.contour(A1,B1,E, levels = 20)
py.clabel(cnt)
py.xlabel('A')
py.ylabel('B')
py.show()


######### q9 #########
# Finding the best fit for column 1
H = sl.lstsq(M,f[:,1])
print(H)
tp=H[0]
print(tp[0])
# Repeating the same method for all columns and finding MS error

# Using scipy.linalg.lstsq to find least squares solution
def lstsq(input):
    return([sl.lstsq(M,input)[0][0], sl.lstsq(M,input)[0][1]])

# Estimate A and B for all sigma values
def manydata():
    error_a = []
    error_b =[]
    for s in sigma:  # Iterating over all sigma values
        est_a = []
        est_b = []
        for i in range(1000): # 1000 data files
            y = np.random.normal(scale=s, size = (101))
            est = lstsq(y+g(f[:,0],1.05,-0.105))
            est_a.append(est[0])
            est_b.append(est[1])
        error_a.append(np.square(np.subtract(1.05,est_a)/1.05).mean())   
        error_b.append(np.square(np.subtract(-0.105,est_b)/-0.105).mean())   # Calculating normalised MSE
 
    return [error_a, error_b]

mse = manydata()
######### q10 #########
# Plot in linear
py.plot(sigma, mse[0], label = 'A')
py.plot(sigma, mse[1], label = 'B')
py.xlabel('sigma')
py.ylabel('MSE Error')
py.legend()
py.title('Error in estimating A, B over 1000 data files (linear)')
py.show()


######### q11 #########
# Plot in loglog
py.loglog(sigma, mse[0], label = 'A')
py.loglog(sigma, mse[1], label = 'B')
py.xlabel('sigma')
py.ylabel('MSE Error')
py.legend()
py.title('Error in estimating A, B over 1000 data files (loglog)')
py.show()