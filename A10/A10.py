import numpy as np
import math

p = 1.03
q = 3.09

Ts = math.sqrt(1/(2*np.pi*p*p))
print(np.exp(-(q*q)/(4*np.pi*Ts*Ts))/Ts)



Ts = math.sqrt(1/(2*np.pi*p*p))
w = q - 2*np.pi*1
print(np.exp(-(w*w)/(4*np.pi*Ts*Ts))/Ts)