import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import triang
from pyDOE import lhs
from scipy.optimize import minimize


#                                      ***    Fixed values    ***
#                                 -------------------------------------
#Times step length (in seconds)
dt = 10000

#creating numpy array for the chemicals
cmax =30

#                                    ***    No of sample runs    ***
#                                 -------------------------------------
# Define the number of Monte Carlo samples
num_samples = 3000
#Creating 2D arrys of num of samples
raws = 100000


#                        ***    Initial Paramters and standard deviations    ***
#                                 -------------------------------------

#      Rivers
#River name acronyms:            Nodes:Dist       Velocities        Dil.frac       Ini.DOC 
#    Teslin - Tes                    x                x                x              x
#    Pelly - Pel                     x                x                x              x
#    White+Donjec - W_D              x                x                x              x
#    Stewart - Ste                   x                x                x              x 
#    Porcupine - Por                 x                x                x              x
#    Tanana - Tan                    x                x                x              x
#    Koyukuk - Koy                   x                x                x              x
#    Yukon - Yuk                     x                x                x           (Teslin)   
#Rivers variables:
#    nodes - Distances 
#    velocities
#    dilution fraction
#    Initial DOC concentrations

def triangular_lhs(num_samples, bounds):
    """
    Generate random samples from a triangular distribution using Latin Hypercube Sampling (LHS).

    Parameters:
    - num_samples: Number of samples to generate.
    - bounds: A list of tuples, where each tuple contains the lower and upper bounds for each parameter.

    Returns:
    - An array representing the Latin Hypercube Samples from a triangular distribution.
    """
    num_parameters = len(bounds)
    samples = lhs(num_parameters, samples=num_samples)
    #print(samples)
    #print(samples.shape)
    #print(samples[0],samples[0].shape)
    
  
    
    # Transform LHS samples to triangular distribution
    sum = []
    for i in range(num_parameters):
        samples[:, i] = triang.ppf(samples[:, i], c=(bounds[i][2] - bounds[i][0]) / (bounds[i][1] - bounds[i][0]),
                                   loc=bounds[i][0], scale=bounds[i][1] - bounds[i][0])
        #print(samples[:, i],i)
        
        #print(samples[:, i].shape)

        #print(samples[0, i],i)
        
    return samples




#Chemical fraction:fraction of DOC 
#Protein-0
prot_f_mean = 4.92/100
prot_f_std = 4.92/100*(1/10)
#Heteroploycondensate(HPC)-3
hpc_f_mean = 83.66/100
hpc_f_std = 83.66/100*(1/10)
#Polysaccharides-4
poly_f_mean = 4.92/100
poly_f_std = 4.92/100*(1/10)
#Lipids-9
lip_f_mean = 4.92/100
lip_f_std = 4.92/100*(1/10)
#Pigments-10
pig_f_mean = 0.98/100
pig_f_std = 0.98/100*(1/10)
#VPhenol-13
V_Phe_f_mean = 0.098/100
V_Phe_f_std = 0.098/100*(1/10)
#SPhenol-14
S_Phe_f_mean = 0.196/100
S_Phe_f_std = 0.196/100*(1/10)
#CPhenol-15
C_Phe_f_mean = 0.294/100
C_Phe_f_std = 0.294/100*(1/10)



#                        ***   Chemical fraction  ***
num_parameters = 8  # Number of parameters or dimensions
 # Bounds for each parameter: (min, max, mode)
bounds = [(prot_f_mean-prot_f_std,prot_f_mean+prot_f_std,prot_f_mean),  
          (hpc_f_mean-hpc_f_std,hpc_f_mean+hpc_f_std,hpc_f_mean),
          (poly_f_mean-poly_f_std,poly_f_mean+poly_f_std,poly_f_mean),
          (lip_f_mean-lip_f_std,lip_f_mean+lip_f_std,lip_f_mean),
          (pig_f_mean-pig_f_std,pig_f_mean+pig_f_std,pig_f_mean),
          (V_Phe_f_mean-V_Phe_f_std,V_Phe_f_mean+V_Phe_f_std,V_Phe_f_mean),
          (S_Phe_f_mean-S_Phe_f_std,S_Phe_f_mean+S_Phe_f_std,S_Phe_f_mean),
          (C_Phe_f_mean-C_Phe_f_std,C_Phe_f_mean+C_Phe_f_std,C_Phe_f_mean)]  # Bounds for each parameter

#print(bounds)
# Generate Latin Hypercube Samples from triangular distribution
samples = triangular_lhs(num_samples, bounds)
#print(samples.shape,"****")
#print(samples)



for i in range(num_samples):
    sum_1 =np.sum(samples[i,:])
    print(sum_1)  
    if np.isclose(sum_1,1):
        print('yes')
    else:
        print('no')
    

#  if np.isclose(np.sum(result.x), target_sum):
#                success[i] = True


for i in range(num_parameters):
    plt.hist(samples[:,i], density=True, alpha=0.6, color='black')
#    plt.hist(adjusted_samples[:,i], density=True, alpha=0.6, color='blue')
    plt.xlabel('Value')
    plt.ylabel('Probability')
    plt.title('Histogram of Triangular Distribution')
    plt.show()










a = np.ones((2,8))
print(a)

x,y = a.shape

success = np.zeros(x, dtype=bool) 
print(success.shape)
print(success) 

for i in range(x):
    print(success[i])


np.all(success)
print(np.all(success))
print(not np.all(success))
#for i in range(8):
#    a[:,i]=i

   
#while not np.all(success):
#   for i in range(x):
#       if success[i]:
#           continue
        
  


#for i in range(x):
#    np.sum(a[i,:])
#    print(np.sum(a[i,:]))    
    
#cons = {'type': 'eq', 'fun': constraint}
bounds = [(0, 1) for _ in range(8)]
print(bounds)
