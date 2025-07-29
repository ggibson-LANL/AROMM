import numpy as np
import math

c = np.arange(0.01,1.01,0.01)
print(c)

d = np.log10(c)-1.0
print(d)

#Efficiency
E = 10**(d)
print(E)

#Life_time=1/(Ek(agg)) k(agg)=1/100 s
Lt=100/(E**2)
print(Lt)

#-dt/Lt
dLt=(-10000)/Lt
print(dLt)

#-dt/Lt, Lt=10**(10)
dLt1=(-10000)/(10**10)
print(dLt1)

#exponential value
print(np.exp(dLt))
print(np.exp(dLt1))
print(10*np.exp(dLt))
print(10*np.exp(dLt1))
print(np.log10(10*np.exp(dLt)))
print(np.log10(10*np.exp(dLt1)))

