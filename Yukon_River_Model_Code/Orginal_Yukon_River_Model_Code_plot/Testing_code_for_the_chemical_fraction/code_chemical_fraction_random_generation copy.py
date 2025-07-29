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


#      Chemical properties
#Chemical acronyms                     Chem.fraction  CDOM.fraction  Chem.Tau     Chem.Production(initial)
#    Protein(0) - prot                      X               X            X              X
#    Polypeptide - polypep                                               X              
#    Amino acid - amino_acid                                             X              
#    Heteropolycondensate(HPC)-3 - hpc      X               X            X              -
#    Polysaccharides-4 - poly               X                            X              X
#    Oligosaccharides-5 - olig                                           X              
#    Humic acid-6 - humic_acid                                           X              
#    Monosaccharides-7 - mono                                            X              
#    Inorganic carbon-8 - ino_C                                          X              
#    Lipids-9 - lip                         X                            X              X
#    Pigments-10 - pig                      X               X            X              X
#    Porphyrin-11 - por                                     X            X              
#    CDOM-12 - CDOM                   
#    V-Phenol-13 - V_Phe                    X               X            X              X
#    S-Phenol-14 - S_Phe                    X               X            X              X
#    C-Phenol-15 - C_Phe                    X               X            X              X
#chemical variables:
#    chemical fraction:fraction of DOC
#    CDOM fraction: frations of CDOM components 
#    Tau values - reaction time in seconds 
#    Prod values - initial production values


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
  
    
    # Transform LHS samples to triangular distribution
    sum = []
    for i in range(num_parameters):
        samples[:, i] = triang.ppf(samples[:, i], c=(bounds[i][2] - bounds[i][0]) / (bounds[i][1] - bounds[i][0]),
                                   loc=bounds[i][0], scale=bounds[i][1] - bounds[i][0])
        #print(samples[:, i],i)
        
        #print(samples[:, i].shape)

        #print(samples[0, i],i)
        
    return samples


def adjust_samples_to_sum(samples, target_sum=1):
    num_samples, num_parameters = samples.shape
    adjusted_samples = np.zeros_like(samples)
    success = np.zeros(num_samples, dtype=bool)

    while not np.all(success):
        for i in range(num_samples):
            if success[i]:
                continue

            def constraint(x):
                return np.sum(x) - target_sum

            cons = {'type': 'eq', 'fun': constraint}
            bounds = [(0, 1) for _ in range(num_parameters)]

            result = minimize(
                lambda x: np.sum((samples[i] - x) ** 2),
                samples[i],
                method='SLSQP',
                bounds=bounds,
                constraints=cons
            )
            adjusted_samples[i] = result.x
            if np.isclose(np.sum(result.x), target_sum):
                success[i] = True

        if not np.all(success):
            samples = triangular_lhs(num_samples, bounds)

    return adjusted_samples


#Chemical acronyms                     
#    Protein(0) - prot                                                
#    Polypeptide - polypep                                              
#    Amino acid - amino_acid                                             
#    Heteroploycondensate(HPC)-3 - hpc                                   
#    Polysaccharides-4 - poly              
#    Oligosaccharides-5 - olig                                           
#    Humic acid-6 - humic_acid                                           
#    Monosaccharides-7 - mono                                             
#    Inorganic carbon-8 - ino_C                                          
#    Lipids-9 - lip                        
#    Pigments-10 - pig                     
#    Porphyrin-11 - por                                               
#    CDOM-12 - CDOM                   
#    V-Phenol-13 - V_Phe                    
#    S-Phenol-14 - S_Phe                   
#    C-Phenol-15 - C_Phe 

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

# Generate Latin Hypercube Samples from triangular distribution
samples = triangular_lhs(num_samples, bounds)
#print(samples.shape)

# Adjust samples to ensure their sum equals 1
adjusted_samples = adjust_samples_to_sum(samples, target_sum=1)
#print(adjusted_samples.shape)
print(np.allclose(np.sum(adjusted_samples, axis=1), 1),"Sum is adding to the 1")

for i in range(num_parameters):
    plt.hist(samples[:,i], density=True, alpha=0.6, color='black')
    plt.xlabel('Value')
    plt.ylabel('Probability')
    plt.title('Histogram of Triangular Distribution')
    plt.show()
    
#Chemical fraction
#Protein
prot_f_samples = adjusted_samples[:,0]
#Heteroplycondenstate
hpc_f_samples = adjusted_samples[:,1]
#Polysaccharide
poly_f_samples = adjusted_samples[:,2]
#Lipids
lip_f_samples = adjusted_samples[:,3]
#Pigments
pig_f_samples = adjusted_samples[:,4]
#V-Phenol
V_Phe_f_samples = adjusted_samples[:,5]
#S-Phenol
S_Phe_f_samples = adjusted_samples[:,6]
#C-Phenol
C_Phe_f_samples = adjusted_samples[:,7]


#Creating input CSV file 
InputData_samples = {
    'Protein_fraction':samples[:,0],
    'Heteroploycondensate(HPC)_fraction':samples[:,1],
    'Polysaccharides_fraction':samples[:,2],
    'Lipids_fraction':samples[:,3],
    'Pigments_fraction':samples[:,4],
    'V-Phenol_fraction':samples[:,5],    
    'S-Phenol_fraction':samples[:,6],   
    'C-Phenol_fraction':samples[:,7]     
          }

df = pd.DataFrame(InputData_samples)
df.to_csv('InputData_samples.csv',index=False)

#Creating input CSV file 
InputData_adjusted = {
    'Protein_fraction':prot_f_samples,
    'Heteroploycondensate(HPC)_fraction':hpc_f_samples,
    'Polysaccharides_fraction':poly_f_samples,
    'Lipids_fraction':lip_f_samples,
    'Pigments_fraction':pig_f_samples,
    'V-Phenol_fraction':V_Phe_f_samples,    
    'S-Phenol_fraction':S_Phe_f_samples,   
    'C-Phenol_fraction':C_Phe_f_samples,     
          }

df = pd.DataFrame(InputData_adjusted)
df.to_csv('InputData_adjusted.csv',index=False)