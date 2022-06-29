import numpy as np
from pylab import *
import scipy.signal as sp

def ForcedOscillation(a,n):
    t = np.linspace(0,50,1000)
    X = sp.lti([1,a],np.polymul([1,0,2.25], np.polyadd(np.polymul([1,a],[1,a]),[2.25])))
    t, x = sp.impulse(X, None, t)

    figure(n)
    title("x(t) time domain")
    ylabel("$x(t)\\rightarrow$")
    xlabel("$t$ (in sec)$\\rightarrow$")
    plot(t, x, label = "Decay = " + str(a))
    grid(True)
    legend()
    show()

ForcedOscillation(0.5,1)
ForcedOscillation(0.05,2)

t = np.linspace(0,50,1000)
H = sp.lti([1],[1,0,2.25])
t, h = sp.impulse(H, None, t)

figure(3)
title("Output with varying frequency. decay coefficient = 0.05 ")
ylabel("$x(t)\\rightarrow$")
xlabel("$t$ (in sec)$\\rightarrow$")
grid(True)

for f in np.linspace(1.4,1.6,5):
    u = np.cos(f*t)*np.exp(-0.05*t)
    t, y, svec = sp.lsim(H, u, t)
    plot(t, y, label = "Frequency = "+ str(f))
legend()
show()


t = np.linspace(0, 20, 1000)
X = sp.lti(np.polymul([1, 0], [0.5, 0, 1]), np.polyadd(np.polymul([1, 0, 1], [0.5, 0, 1]), [-1]))
Y = sp.lti([1, 0], np.polyadd(np.polymul([1, 0, 1], [0.5, 0, 1]), [-1]))

#Time domain function values of X(s) and Y(s)
t, x = sp.impulse(X, None, t)
t, y = sp.impulse(Y, None, t)

#Plotting functions
figure(4)
plot(t, x, label="x(t)")
plot(t, y, label="y(t)")
title("x(t) & y(t) in time domain")
ylabel("$Signal\\rightarrow$")
xlabel("$t$ (in sec)$\\rightarrow$")
legend()
grid(True)
show() 

# Transfer functions for H(s)
H = sp.lti([1], [1e-12, 1e-4, 1])

#Obtaining bode plot
w, S, phi = H.bode()

#Plotting functions
figure(num=5,figsize=(7,10))
plot1 = subplot(2, 1, 1)
plot1.semilogx(w, S)
plot1.set_title("Magnitude Plot")
plot1.set_ylabel("Magnitude(in dB)$\\rightarrow$")
plot1.set_xlabel("$\\omega\\rightarrow$")
grid(True)

plot2 = subplot(2, 1, 2)
plot2.semilogx(w, phi)
plot2.set_title("Phase Plot")
plot2.set_ylabel("Phase$\\rightarrow$")
plot2.set_xlabel("$\\omega\\rightarrow$")
grid(True)
show()

def LCRresponse(t,n):
    #Intitial short term response 
    Vin = np.cos(1e3*t)-np.cos(1e6*t)
    t, Vout, svec = sp.lsim(H, Vin, t)  #Convolution for initial response

    #Plotting functions
    figure(n)
    plot(t, Vout)
    title("$V_{out}$ vs $t$ (Initial response)")
    xlabel("$t$ (in sec)$\\rightarrow$")
    ylabel("$V_{out}\\rightarrow$")
    grid(True)
    show()

#Intitial short term response
t = np.linspace(0, 30e-6, 10000)
LCRresponse(t,6)

#Long term response
t = np.linspace(0, 1e-2, 10000)
LCRresponse(t,7)