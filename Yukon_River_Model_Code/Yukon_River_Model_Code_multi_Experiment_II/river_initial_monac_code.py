import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import triang
from scipy.stats import norm
from scipy.stats import uniform
from pyDOE import lhs
from scipy.optimize import minimize
import copy


Trib_Name = ['Teslin','Pelly','White+Donjec','Stewart','Porcupine','Tanana','Koyukuk','Yukon']

#                                      ***    Fixed values    ***
#                                 -------------------------------------
#Times step length (in seconds)
dt = 10000

#creating numpy array for the chemicals
cmax =30

#                                    ***    No of sample runs    ***
#                                 -------------------------------------
# Define the number of Monte Carlo samples
num_samples = 10
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

#Nodes: distance in km
#for tribuatries:1st 100 km, total distance
#for the main river: 1st 100 km, distance where tributaries meet the main stem, and total distance
#Teslin
Tes_nodes = [0, 100.0, 374.0]
#Pelly   
Pel_nodes = [0, 100.0, 475.0] 
#White+Donjec  
W_D_nodes = [0, 100.0, 316.7] 
#Stewart  
Ste_nodes = [0, 100.0, 583.5] 
#Porcupine  
Por_nodes = [0, 100.0, 783.6] 
#Tanana  
Tan_nodes = [0, 100.0, 884.5]
#Koyukuk   
Koy_nodes = [0, 100.0, 859.9] 
#Yukon  
#Yuk_nodes = [Tes(MP):0, 100.0, Pel(CP1):234.0, W_D(CP2):383.8, Ste(CP3):399.5, 
#             Por(CP4):1009.6, Tan(CP5):1473.4, Koy(CP6)1799.4, (DR)2457.8, 2613.3] 
Yuk_nodes = [0, 100.0, 234.0, 383.8, 399.5, 1009.6, 1473.4, 1799.4, 2457.8, 2613.3]   

#Traingular distribution
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


#Normal distribution #Scaling the lhs to avoid negative values
def normal_lhs(num_samples, bounds,CDF_min,CDF_max):
    """
    Generate random samples from a normal distribution using Latin Hypercube Sampling (LHS).

    Parameters:
    - num_samples: Number of samples to generate.
    - bounds: A list of tuples, where each tuple contains the mean and standard deviation for each parameter.

    Returns:
    - An array representing the Latin Hypercube Samples from a normal distribution.
    """
    num_parameters = len(bounds)
    #Scaling the lhs . a = minimum value that trauncate the normal distribution from the 0. max = 0.99
    samples = CDF_min + lhs(num_parameters, samples=num_samples)*(CDF_max - CDF_min)

    # Transform LHS samples to normal distribution
    for i in range(num_parameters):
        mean, std = bounds[i]
        samples[:, i] = norm.ppf(samples[:, i], loc=mean, scale=std)
        # Truncate negative values to zero
        #samples[:, i] = np.maximum(samples[:, i], 10**(-2))

    return samples

#Unifrom distribution
def uniform_lhs(num_samples, bounds):
    """
    Generate random samples from a uniform distribution using Latin Hypercube Sampling (LHS).

    Parameters:
    - num_samples: Number of samples to generate.
    - bounds: A list of tuples, where each tuple contains the lower and upper bounds for each parameter.

    Returns:
    - An array representing the Latin Hypercube Samples from a uniform distribution.
    """
    num_parameters = len(bounds)
    samples = lhs(num_parameters, samples=num_samples)

    # Transform LHS samples to uniform distribution
    for i in range(num_parameters):
        lower, upper = bounds[i]
        samples[:, i] = uniform.ppf(samples[:, i], loc=lower, scale=upper - lower)
        # Truncate negative values to zero
        samples[:, i] = np.maximum(samples[:, i], 0)

    return samples



#Changing the standard deviation values to the actual deviation values from the calculation
#Normal distribution considered
#Velocities: in ms-1
#V100 for the 1st 100 km 
#Teslin
Tes_V100_mean = 0.1
Tes_V100_std = 0.00
Tes_V1_mean = 0.56                   
Tes_V1_std = 0.50
#Pelly
Pel_V100_mean = 0.1
Pel_V100_std = 0.000
Pel_V1_mean = 0.73
Pel_V1_std = 0.57
#White+Donjec
W_D_V100_mean = 0.1
W_D_V100_std = 0.000
W_D_V1_mean = 0.52
W_D_V1_std = 0.63
#Stewart
Ste_V100_mean = 0.1
Ste_V100_std = 0.000
Ste_V1_mean = 0.55
Ste_V1_std = 0.55
#Porcupine
Por_V100_mean = 0.1
Por_V100_std = 0.000
Por_V1_mean = 0.72
Por_V1_std = 0.76
#Tanana
Tan_V100_mean = 0.1
Tan_V100_std = 0.000
Tan_V1_mean = 0.91
Tan_V1_std = 0.55
#Koyukuk
Koy_V100_mean = 0.1
Koy_V100_std = 0.000
Koy_V1_mean = 0.72
Koy_V1_std = 0.53
#Yukon # For all the sections in the Yukon river we assumes they have a same velocity value
Yuk_V100_mean = 0.1   #first 100 km
Yuk_V100_std = 0.000
Yuk_V1_mean = 0.87     #from 100 km to Pelly R meeting point(CP1) 
Yuk_V1_std = 0.56
Yuk_V2_mean = 0.87     #from CP1 to White+D R meeting point(CP2) 
Yuk_V2_std = 0.56
Yuk_V3_mean = 0.87     #from CP2 to Stewart R meeting point(CP3) 
Yuk_V3_std = 0.56
Yuk_V4_mean = 0.87     #from CP3 to Porcupine R meeting point(CP4) 
Yuk_V4_std = 0.56
Yuk_V5_mean = 0.87     #from CP4 to Tanana R meeting point(CP5) 
Yuk_V5_std = 0.56
Yuk_V6_mean = 0.87     #from CP5 to Koyukuk R meeting point(CP6) 
Yuk_V6_std = 0.56
Yuk_V7_mean = 0.87     #from CP6 to Delta region starting point(DR) 
Yuk_V7_std = 0.56
Yuk_V8_mean = 0.1     #from DR to river end point (Delta region)
Yuk_V8_std = 0.00

#Dilution fraction(DF) for the tributaries
#normal distribution considered
#Teslin
Tes_Df_mean = 0.0     # no dilution happen at the merging point 
Tes_Df_std = 0.0
#Pelly   
Pel_Df_mean = 0.29
Pel_Df_std = 0.12
#White+Donjec  
W_D_Df_mean = 0.064
W_D_Df_std = 0.043
#Stewart  
Ste_Df_mean = 0.17
Ste_Df_std = 0.032 
#Porcupine  
Por_Df_mean = 0.090
Por_Df_std = 0.063
#Tanana  
Tan_Df_mean = 0.18
Tan_Df_std = 0.047
#Koyukuk   
Koy_Df_mean = 0.060
Koy_Df_std = 0.041

#Initial DOC concentration: in micromolar C (micromol C / L)
#Max of 5% of the baseline was considered as standard deviation
#Normal distribution considered
#Teslin
Tes_DOC_mean = 387.21  #169.4 
Tes_DOC_std = 387.21*(25/100)
#Pelly   
Pel_DOC_mean = 680.86
Pel_DOC_std = 680.86*(25/100)
#White+Donjec  
W_D_DOC_mean = 662.41
W_D_DOC_std = 662.41*(25/100)
#Stewart  
Ste_DOC_mean = 678.31
Ste_DOC_std = 678.31*(25/100)
#Porcupine  
Por_DOC_mean = 502.63
Por_DOC_std = 502.63*(25/100)
#Tanana  
Tan_DOC_mean = 315.94
Tan_DOC_std = 315.94*(25/100)
#Koyukuk   
Koy_DOC_mean = 198.54
Koy_DOC_std = 198.54*(25/100)

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
#Remove the fraction bound rest of the chemicals except hpc
#Traingular distribution considered 
#10% of the baseline
#Protein-0
prot_f_mean = 4.92/100
prot_f_std = 4.92/100*(1/10)
#Heteroploycondensate(HPC)-3
hpc_f_mean = 83.66/100               # this is the one that actually changing by 10% of the mean
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

#CDOM component fractions:fraction of CDOM components to create CDOM
#CDOM = prot_f_CDOM *Protein 
#       + hpc_f_CDOM * Heteroploycondensate(HPC) 
#       + pig_f_CDOM * Pigments 
#       + por_f_CDOM * Porphyrin
#       + V_Phe_f_CDOM * V-Phenol 
#       + S_Phe_f_CDOM * S-Phenol 
#       + C_Phe_f_CDOM * C-Phenol 
#Protein-0
prot_f_CDOM_mean = 0.1  
prot_f_CDOM_std = 0.00
#Heteroploycondensate(HPC)-3
hpc_f_CDOM_mean = 0.1
hpc_f_CDOM_std = 0.00
#Pigments-10
pig_f_CDOM_mean = 1.0
pig_f_CDOM_std = 0.00
#Porphyrin-11
por_f_CDOM_mean = 1.0     # new addition
por_f_CDOM_std = 0.00
#VPhenol-13
V_Phe_f_CDOM_mean = 0.33
V_Phe_f_CDOM_std = 0.00
#SPhenol-14
S_Phe_f_CDOM_mean = 0.33
S_Phe_f_CDOM_std = 0.00
#CPhenol-15
C_Phe_f_CDOM_mean = 0.33
C_Phe_f_CDOM_std = 0.00


#Tau values - reaction time in seconds
#10% of the baseline was considered
#Uniform distribution considered.
#protein-0
prot_Tau_mean =8640000              
prot_Tau_std = 8640000*(1/10)          
#Polypeptide-1
polypep_Tau_mean = 864000*(3/13)*(10)     
polypep_Tau_std = 864000*(3/13)*(10)*(1/10)               
#Amino acid-2
amino_acid_Tau_mean = 864000*(1/10)*(10)  
amino_acid_Tau_std = 864000*(1/10)*(10)*(1/10)                 
#Heteroploycondensate(HPC)-3
hpc_Tau_mean = 864000*(1000000)
hpc_Tau_std = 864000*(1000000)*(1/10) 
#Polysaccharides-4
poly_Tau_mean = 864000*3*(10)
poly_Tau_std = 864000*3*(10)*(1/10)
#Oligosaccharides-5
olig_Tau_mean = (10/11)*864000*(10)
olig_Tau_std = (10/11)*864000*(10)*(1/10)
#Humic acid-6
humic_acid_Tau_mean = 864000*3*(10)
humic_acid_Tau_std = 864000*3*(10)*(1/10)
#Monosaccharides-7
mono_Tau_mean = 864000*3*(10)
mono_Tau_std = 864000*3*(10)*(1/10)
#Inorganic carbon-8
ino_C_Tau_mean = 86400*3*(10)
ino_C_Tau_std =86400*3*(10)*(1/10)
#Lipids-9
lip_Tau_mean = 864000*(3/13)*(10)
lip_Tau_std = 864000*(3/13)*(10)*(1/10)
#Pigment-10
pig_Tau_mean = 864000*(10)
pig_Tau_std = 864000*(10)*(1/10)
#Porphyrin-11
por_Tau_mean = 864000*(10)
por_Tau_std = 864000*(10)*(1/10)
#VPhenol-13
V_Phe_Tau_mean = 86400*3*(10)
V_Phe_Tau_std = 86400*3*(10)*(1/10)
#SPhenol-14
S_Phe_Tau_mean = 86400*3*(10)
S_Phe_Tau_std = 86400*3*(10)*(1/10)
#CPhenol-15
C_Phe_Tau_mean = 86400*3*(10)
C_Phe_Tau_std = 86400*3*(10)*(1/10)

#Prod values - Initial production values
#Protein-0
prot_Prod_mean = 0.00
prot_Prod_std = 0.00
#Polysaccharides-4
poly_Prod_mean = 0.00
poly_Prod_std = 0.00
#Lipids-9
lip_Prod_mean = 0.00
lip_Prod_std = 0.00
#Pigments-10
pig_Prod_mean = 0.00
pig_Prod_std = 0.00
#VPhenol-13
V_Phe_Prod_mean = 0.00
V_Phe_Prod_std = 0.00
#SPhenol-14
S_Phe_Prod_mean = 0.000
S_Phe_Prod_std = 0.00
#CPhenol-15
C_Phe_Prod_mean = 0.00
C_Phe_Prod_std = 0.00






#                        ***   Random Sample genration   ***
#                          -------------------------------

#Variables                     Random Sample generation
#   Nodes:Dist                             
#   Velocities                          x
#   Dil.frac                            x
#   Ini.DOC                             x
#   Chem.fraction                       x
#   Chem.Tau                            x
#   Chem.Production(initial)            x (this is mostly zero)
#River acronyms 
#    Teslin - Tes                
#    Pelly - Pel                   
#    White+Donjec - W_D             
#    Stewart - Ste                  
#    Porcupine - Por                
#    Tanana - Tan                  
#    Koyukuk - Koy                  
#    Yukon - Yuk


#                        ***   Velocity  ***
num_parameters = 8  # Number of parameters or dimensions
 # Bounds for each parameter: (min, max, mode)
bounds = [(Tes_V1_mean, Tes_V1_std),
          (Pel_V1_mean, Pel_V1_std),
          (W_D_V1_mean, W_D_V1_std),
          (Ste_V1_mean, Ste_V1_std),
          (Por_V1_mean, Por_V1_std),
          (Tan_V1_mean, Tan_V1_std),
          (Koy_V1_mean, Koy_V1_std),
          (Yuk_V1_mean, Yuk_V1_std)]  # Bounds for each parameter


CDF_min = np.zeros(len(bounds))
CDF_max = np.zeros(len(bounds))
for i in range(num_parameters):
    x = np.linspace(-5,5, 1000)
# Compute PDF values
    pdf_values = norm.pdf(x, loc=bounds[i][0], scale=bounds[i][1])
    # Plot PDF
    plt.plot(x, pdf_values, 'k-', label='PDF')
    plt.title(Trib_Name[i]+'- Normal Distribution PDF - Velocity')
    plt.axvline(x=0, color='green', linestyle='--', label='x =0')
    plt.axvline(x=bounds[i][0]-bounds[i][1], color='k', linestyle='--', label=f'x = {bounds[i][0]-bounds[i][1]}')          #
    plt.axvline(x=bounds[i][0]+bounds[i][1], color='r', linestyle='--', label=f'x = {bounds[i][0]+bounds[i][1]}')          #
    plt.xlabel('x')
    plt.ylabel('Probability Density')
    plt.legend()
    plt.show()    

    #Computing CDF value at 0
    cdf_at_zero = norm.cdf(0, loc=bounds[i][0], scale=bounds[i][1])        #
    plt.plot(0, cdf_at_zero, 'ro', label=f'cdf_val = {cdf_at_zero}')
    plt.axvline(x=0, color='green', linestyle='--', label='x = 0') 


    #Compution CDF value at the lower limit
    if bounds[i][0]-bounds[i][1] < 0:                                              #
        # Compute CDF values
        cdf_low_limit = norm.cdf(0, loc=bounds[i][0], scale=bounds[i][1])          #
        #print(i,cdf_low_limit,'cdf_low_limit')
        CDF_min[i] = cdf_low_limit  # Correct assignment                           #
        #CDF_min[i] = cdf_low_limit
        plt.plot(0, cdf_low_limit, 'ro', label=f'cdf_val = {cdf_low_limit}')
        plt.axvline(x=0, color='k', linestyle='--', label='x=0')
    else:
        # Compute CDF values
        cdf_low_limit = norm.cdf(bounds[i][0]-bounds[i][1], loc=bounds[i][0], scale=bounds[i][1])     #
        #print(i,cdf_low_limit,'cdf_low_limit')
        CDF_min[i] = cdf_low_limit  # Correct assignment                                              #
        #CDF_min[i] = cdf_low_limit 
        plt.plot(bounds[i][0]-bounds[i][1], cdf_low_limit, 'ro', label=f'cdf_val = {cdf_low_limit}')  #
        plt.axvline(x=bounds[i][0]-bounds[i][1], color='k', linestyle='--', label=f'x= {bounds[i][0]-bounds[i][1]}') #      

    #Compute the CDF value to the lower limit
    # Compute CDF values
    cdf_upper_limit = norm.cdf(bounds[i][0]+bounds[i][1], loc=bounds[i][0], scale=bounds[i][1])     #
    #print(i,cdf_upper_limit,'cdf_upper_limit')
    #CDF_max[i] = cdf_upper_limit  # Correct assignment 
    # changing the CDF upper limit
    #changing the Cummulative functions upper limit
    CDF_max[i] = 0.9999  # Correct assignment                                              #
    #CDF_min[i] = cdf_low_limit 
    plt.plot(bounds[i][0]+bounds[i][1], cdf_upper_limit, 'ro', label=f'cdf_val = {cdf_upper_limit}')  #
    plt.axvline(x=bounds[i][0]+bounds[i][1], color='r', linestyle='--', label=f'x= {bounds[i][0]+bounds[i][1]}') #     

    #print("Difference",CDF_max[i]-CDF_min[i])

    # Plot CDF
    cdf_values = norm.cdf(x, loc=bounds[i][0], scale=bounds[i][1]) 
    plt.plot(x, cdf_values, 'k-', label='CDF')
    plt.title(Trib_Name[i]+'-Normal Distribution CDF - Velocity')         #
    plt.xlabel('x')
    plt.ylabel('Cumulative Probability')
    plt.legend()
    plt.show()

# Generate Latin Hypercube Samples from normal distribution. Lower limit CDF_min to CDF_max,trunc the negative values by setting lower value to zero
samples = normal_lhs(num_samples, bounds,CDF_min,CDF_max)
print(samples.shape)


for i in range(num_parameters):
    # Generate x values
    x = np.linspace(-5,5, 1000)
    figure = plt.figure()
    # Compute PDF values & plot PDF    
    pdf_values = norm.pdf(x, loc=bounds[i][0], scale=bounds[i][1])
    plt.plot(x, pdf_values, 'k-', label='PDF')
    #
    plt.hist(samples[:,i], density=True, alpha=0.6, color='red')
    plt.xlabel('Value')
    plt.ylabel('Probability')    
    plt.title(Trib_Name[i]+' - Histogram of Normal Distribution - Velocity')
    plt.axvline(bounds[i][0]-bounds[i][1], color='k', linestyle='-.', label=f'$\mu$ - $\sigma$ = {bounds[i][0]-bounds[i][1]:.3f}')
    plt.axvline(bounds[i][0], color='green', linestyle='--', label=f'$\mu$ = {bounds[i][0]}')
    plt.axvline(bounds[i][0]+bounds[i][1], color='k', linestyle='--', label=f'$\mu$ + $\sigma$ = {bounds[i][0]+bounds[i][1]:.3f}') 
    plt.legend(loc='lower center', bbox_to_anchor = (0.5, -0.3),ncol =4)
    plt.tight_layout()
    plt.show()
    figure.savefig(f"output/analysis/Figures/1_Velocity_Histogram of Normal Distribution_{Trib_Name[i]}.png")
    #figure.savefig(f"/Users/amadini/Documents/Current Working folder/MONACO/01_Jan_Monaco/MONACO_Transfer_Docs_Codes/03_24_2024/TESTS_06_10_2024/Annual/0_30_cm_chemical_fraction_code_changed/Annual_YR_Sites_Avg_3(Prod_Coe)_chemical_fraction_code_changed/Final_Monte_Carlo_Analysis_Code/Exp_02_multi%_derivations_norm_uniform/TEST_01_New_Code_fractions_multi%_annual_chemical_fraction_code_changed_GG_suggestion_for_fraction_added_norm_code_changed_AGAIN_DOC_5%/output/analysis{Trib_Name[i]}+ - Histogram of Normal Distribution - Velocity.png")

#Velocities
#Teslin
Tes_V100_samples = np.random.normal(Tes_V100_mean, Tes_V100_std, num_samples)
Tes_V1_samples = samples[:,0]
#print(Tes_V1_samples)
#Pelly
Pel_V100_samples = np.random.normal(Pel_V100_mean, Pel_V100_std, num_samples)
Pel_V1_samples = samples[:,1]
#White+Donjec
W_D_V100_samples = np.random.normal(W_D_V100_mean, W_D_V100_std, num_samples)
W_D_V1_samples = samples[:,2]
#Stewart  
Ste_V100_samples = np.random.normal(Ste_V100_mean, Ste_V100_std, num_samples)
Ste_V1_samples = samples[:,3]
#Porcupine 
Por_V100_samples = np.random.normal(Por_V100_mean, Por_V100_std, num_samples)
Por_V1_samples = samples[:,4]
#Tanana 
Tan_V100_samples = np.random.normal(Tan_V100_mean, Tan_V100_std, num_samples)
Tan_V1_samples = samples[:,5]
#Koyukuk
Koy_V100_samples = np.random.normal(Koy_V100_mean, Koy_V100_std, num_samples)
Koy_V1_samples = samples[:,6]
#Yukon - Yuk
Yuk_V100_samples = np.random.normal(Yuk_V100_mean, Yuk_V100_std, num_samples)
Yuk_V1_samples = samples[:,7]
Yuk_V2_samples = samples[:,7]
Yuk_V3_samples = samples[:,7]
Yuk_V4_samples = samples[:,7]
Yuk_V5_samples = samples[:,7]
Yuk_V6_samples = samples[:,7]
Yuk_V7_samples = samples[:,7]
Yuk_V8_samples = samples[:,7]

#                        ***   Dilution_factor  ***
num_parameters = 6  # Number of parameters or dimensions
 # Bounds for each parameter: (min, max, mode)
bounds = [(Pel_Df_mean,Pel_Df_std),
          (W_D_Df_mean,W_D_Df_std),
          (Ste_Df_mean,Ste_Df_std),
          (Por_Df_mean,Por_Df_std),
          (Tan_Df_mean,Tan_Df_std),
          (Koy_Df_mean,Koy_Df_std)]  # Bounds for each parameter


CDF_min = np.zeros(len(bounds))
CDF_max = np.zeros(len(bounds))
for i in range(num_parameters):
    x = np.linspace(-0.2,0.8, 1000)
# Compute PDF values
    pdf_values = norm.pdf(x, loc=bounds[i][0], scale=bounds[i][1])
    # Plot PDF
    plt.plot(x, pdf_values, 'k-', label='PDF')
    plt.title(Trib_Name[i+1]+'- Normal Distribution PDF - Dilution factor')
    plt.axvline(x=0, color='green', linestyle='--', label='x =0')
    plt.axvline(x=bounds[i][0]-bounds[i][1], color='k', linestyle='--', label=f'x = {bounds[i][0]-bounds[i][1]}')          #
    plt.axvline(x=bounds[i][0]+bounds[i][1], color='r', linestyle='--', label=f'x = {bounds[i][0]+bounds[i][1]}')          #
    plt.xlabel('x')
    plt.ylabel('Probability Density')
    plt.legend()
    plt.show()    

    #Computing CDF value at 0
    cdf_at_zero = norm.cdf(0, loc=bounds[i][0], scale=bounds[i][1])        #
    plt.plot(0, cdf_at_zero, 'ro', label=f'cdf_val = {cdf_at_zero}')
    plt.axvline(x=0, color='green', linestyle='--', label='x = 0') 


    #Compution CDF value at the lower limit
    if bounds[i][0]-bounds[i][1] < 0:                                              #
        # Compute CDF values
        cdf_low_limit = norm.cdf(0, loc=bounds[i][0], scale=bounds[i][1])          #
        #print(i,cdf_low_limit,'cdf_low_limit')
        CDF_min[i] = cdf_low_limit  # Correct assignment                           #
        #CDF_min[i] = cdf_low_limit
        plt.plot(0, cdf_low_limit, 'ro', label=f'cdf_val = {cdf_low_limit}')
        plt.axvline(x=0, color='k', linestyle='--', label='x=0')
    else:
        # Compute CDF values
        cdf_low_limit = norm.cdf(bounds[i][0]-bounds[i][1], loc=bounds[i][0], scale=bounds[i][1])     #
        #print(i,cdf_low_limit,'cdf_low_limit')
        CDF_min[i] = cdf_low_limit  # Correct assignment                                              #
        #CDF_min[i] = cdf_low_limit 
        plt.plot(bounds[i][0]-bounds[i][1], cdf_low_limit, 'ro', label=f'cdf_val = {cdf_low_limit}')  #
        plt.axvline(x=bounds[i][0]-bounds[i][1], color='k', linestyle='--', label=f'x= {bounds[i][0]-bounds[i][1]}') #      

    #Compute the CDF value to the lower limit
    # Compute CDF values
    cdf_upper_limit = norm.cdf(bounds[i][0]+bounds[i][1], loc=bounds[i][0], scale=bounds[i][1])     #
    #print(i,cdf_upper_limit,'cdf_upper_limit')
    #CDF_max[i] = cdf_upper_limit  # Correct assignment    
    CDF_max[i] = 0.9999                                           #
    #CDF_min[i] = cdf_low_limit 
    plt.plot(bounds[i][0]+bounds[i][1], cdf_upper_limit, 'ro', label=f'cdf_val = {cdf_upper_limit}')  #
    plt.axvline(x=bounds[i][0]+bounds[i][1], color='r', linestyle='--', label=f'x= {bounds[i][0]+bounds[i][1]}') #     

    #print("Difference",CDF_max[i]-CDF_min[i])

    # Plot CDF
    cdf_values = norm.cdf(x, loc=bounds[i][0], scale=bounds[i][1]) 
    plt.plot(x, cdf_values, 'k-', label='CDF')
    plt.title(Trib_Name[i]+'-Normal Distribution CDF - Velocity')         #
    plt.xlabel('x')
    plt.ylabel('Cumulative Probability')
    plt.legend()
    plt.show()

# Generate Latin Hypercube Samples from normal distribution. Lower limit CDF_min to CDF_max,trunc the negative values by setting lower value to zero
samples = normal_lhs(num_samples, bounds,CDF_min,CDF_max)
print(samples.shape)


for i in range(num_parameters):
    # Generate x values
    x = np.linspace(-0.2,0.8, 1000)
    figure = plt.figure()
    # Compute PDF values & plot PDF
    pdf_values = norm.pdf(x, loc=bounds[i][0], scale=bounds[i][1])
    plt.plot(x, pdf_values, 'k-', label='PDF')
    #
    plt.hist(samples[:,i], density=True, alpha=0.6, color='red')
    plt.xlabel('Value')
    plt.ylabel('Probability')    
    plt.title(Trib_Name[i+1]+' - Histogram of Normal Distribution - Dilution factor')
    plt.axvline(bounds[i][0]-bounds[i][1], color='k', linestyle='-.', label=f'$\mu$ - $\sigma$ = {bounds[i][0]-bounds[i][1]:.3f}')
    plt.axvline(bounds[i][0], color='green', linestyle='--', label=f'$\mu$ = {bounds[i][0]}')
    plt.axvline(bounds[i][0]+bounds[i][1], color='k', linestyle='--', label=f'$\mu$ + $\sigma$ = {bounds[i][0]+bounds[i][1]:.3f}') 
    plt.legend(loc='lower center', bbox_to_anchor = (0.5, -0.3),ncol =4)
    plt.tight_layout()
    plt.show()
    figure.savefig(f"output/analysis/Figures/2_Dilution_factor_Histogram of Normal Distribution_{Trib_Name[i+1]}.png")

#Dilution fractions
#Teslin
Tes_Df_samples = np.random.normal(Tes_Df_mean, Tes_Df_std, num_samples)
#Pelly
Pel_Df_samples = samples[:,0]
#White+Donjec
W_D_Df_samples = samples[:,1]
#Stewart  
Ste_Df_samples = samples[:,2]
#Porcupine - Por
Por_Df_samples = samples[:,3]
#Tanana - Tan 
Tan_Df_samples = samples[:,4]
#Koyukuk - Koy
Koy_Df_samples = samples[:,5]


#                        ***   DOC  ***
num_parameters = 7  # Number of parameters or dimensions
 # Bounds for each parameter: (min, max, mode)
bounds = [(Tes_DOC_mean,Tes_DOC_std),
          (Pel_DOC_mean,Pel_DOC_std),
          (W_D_DOC_mean,W_D_DOC_std),
          (Ste_DOC_mean,Ste_DOC_std),
          (Por_DOC_mean,Por_DOC_std),
          (Tan_DOC_mean,Tan_DOC_std),
          (Koy_DOC_mean,Koy_DOC_std)]  # Bounds for each parameter


CDF_min = np.zeros(len(bounds))
CDF_max = np.zeros(len(bounds))
for i in range(num_parameters):
    x = np.linspace(-0.5,1000, 1000)
# Compute PDF values
    pdf_values = norm.pdf(x, loc=bounds[i][0], scale=bounds[i][1])
    # Plot PDF
    plt.plot(x, pdf_values, 'k-', label='PDF')
    plt.title(Trib_Name[i]+'- Normal Distribution PDF - DOC')
    plt.axvline(x=0, color='green', linestyle='--', label='x =0')
    plt.axvline(x=bounds[i][0]-bounds[i][1], color='k', linestyle='--', label=f'x = {bounds[i][0]-bounds[i][1]}')          #
    plt.axvline(x=bounds[i][0]+bounds[i][1], color='r', linestyle='--', label=f'x = {bounds[i][0]+bounds[i][1]}')          #
    plt.xlabel('x')
    plt.ylabel('Probability Density')
    plt.legend()
    plt.show()    

    #Computing CDF value at 0
    cdf_at_zero = norm.cdf(0, loc=bounds[i][0], scale=bounds[i][1])        #
    plt.plot(0, cdf_at_zero, 'ro', label=f'cdf_val = {cdf_at_zero}')
    plt.axvline(x=0, color='green', linestyle='--', label='x = 0') 


    #Compution CDF value at the lower limit
    if bounds[i][0]-bounds[i][1] < 0:                                              #
        # Compute CDF values
        cdf_low_limit = norm.cdf(0, loc=bounds[i][0], scale=bounds[i][1])          #
        #print(i,cdf_low_limit,'cdf_low_limit')
        CDF_min[i] = cdf_low_limit  # Correct assignment                           #
        #CDF_min[i] = cdf_low_limit
        plt.plot(0, cdf_low_limit, 'ro', label=f'cdf_val = {cdf_low_limit}')
        plt.axvline(x=0, color='k', linestyle='--', label='x=0')
    else:
        # Compute CDF values
        cdf_low_limit = norm.cdf(bounds[i][0]-bounds[i][1], loc=bounds[i][0], scale=bounds[i][1])     #
        #print(i,cdf_low_limit,'cdf_low_limit')
        CDF_min[i] = cdf_low_limit  # Correct assignment                                              #
        #CDF_min[i] = cdf_low_limit 
        plt.plot(bounds[i][0]-bounds[i][1], cdf_low_limit, 'ro', label=f'cdf_val = {cdf_low_limit}')  #
        plt.axvline(x=bounds[i][0]-bounds[i][1], color='k', linestyle='--', label=f'x= {bounds[i][0]-bounds[i][1]}') #      

    #Compute the CDF value to the lower limit
    # Compute CDF values
    cdf_upper_limit = norm.cdf(bounds[i][0]+bounds[i][1], loc=bounds[i][0], scale=bounds[i][1])     #
    #print(i,cdf_upper_limit,'cdf_upper_limit')
    CDF_max[i] = cdf_upper_limit  # Correct assignment  
    CDF_max[i] = 0.9999                                            #
    #CDF_min[i] = cdf_low_limit 
    plt.plot(bounds[i][0]+bounds[i][1], cdf_upper_limit, 'ro', label=f'cdf_val = {cdf_upper_limit}')  #
    plt.axvline(x=bounds[i][0]+bounds[i][1], color='r', linestyle='--', label=f'x= {bounds[i][0]+bounds[i][1]}') #     

    #print("Difference",CDF_max[i]-CDF_min[i])

    # Plot CDF
    cdf_values = norm.cdf(x, loc=bounds[i][0], scale=bounds[i][1]) 
    plt.plot(x, cdf_values, 'k-', label='CDF')
    plt.title(Trib_Name[i]+'-Normal Distribution CDF - Velocity')         #
    plt.xlabel('x')
    plt.ylabel('Cumulative Probability')
    plt.legend()
    plt.show()

# Generate Latin Hypercube Samples from normal distribution. Lower limit CDF_min to CDF_max,trunc the negative values by setting lower value to zero
samples = normal_lhs(num_samples, bounds,CDF_min,CDF_max)
print(samples.shape)


for i in range(num_parameters):
    # Generate x values
    x = np.linspace(-0.5,1000, 1000)
    figure = plt.figure()
    # Compute PDF values & plot PDF
    pdf_values = norm.pdf(x, loc=bounds[i][0], scale=bounds[i][1])
    plt.plot(x, pdf_values, 'k-', label='PDF')
    #
    plt.hist(samples[:,i], density=True, alpha=0.6, color='red')
    plt.xlabel('Value')
    plt.ylabel('Probability')    
    plt.title(Trib_Name[i]+' - Histogram of Normal Distribution - DOC')
    plt.axvline(bounds[i][0]-bounds[i][1], color='k', linestyle='-.', label=f'$\mu$ - $\sigma$ = {bounds[i][0]-bounds[i][1]:.3f}')
    plt.axvline(bounds[i][0], color='green', linestyle='--', label=f'$\mu$ = {bounds[i][0]}')
    plt.axvline(bounds[i][0]+bounds[i][1], color='k', linestyle='--', label=f'$\mu$ + $\sigma$ = {bounds[i][0]+bounds[i][1]:.3f}') 
    plt.legend(loc='lower center', bbox_to_anchor = (0.5, -0.3),ncol =4)
    plt.tight_layout()
    plt.show()
    figure.savefig(f"output/analysis/Figures/3_DOC_Histogram of Normal Distribution_{Trib_Name[i]}.png")
    
#Initial DOC concentrations
#Teslin
Tes_DOC_samples = samples[:,0]
#Pelly
Pel_DOC_samples = samples[:,1]
#White+Donjec
W_D_DOC_samples = samples[:,2]
#Stewart  
Ste_DOC_samples = samples[:,3]
#Porcupine - Por
Por_DOC_samples = samples[:,4]
#Tanana - Tan 
Tan_DOC_samples = samples[:,5]
#Koyukuk - Koy
Koy_DOC_samples = samples[:,6]

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

#Remove the fraction bound rest of the chemicals except hpc. Only the HPC will be change in the bounds
num_parameters = 1  # Number of parameters or dimensions
 # Bounds for each parameter: (min, max, mode)
bounds = [#(prot_f_mean-prot_f_std,prot_f_mean+prot_f_std,prot_f_mean),  
          (hpc_f_mean-hpc_f_std, hpc_f_mean+hpc_f_std)
          #(poly_f_mean-poly_f_std,poly_f_mean+poly_f_std,poly_f_mean),
          #(lip_f_mean-lip_f_std,lip_f_mean+lip_f_std,lip_f_mean),
          #(pig_f_mean-pig_f_std,pig_f_mean+pig_f_std,pig_f_mean),
          #(V_Phe_f_mean-V_Phe_f_std,V_Phe_f_mean+V_Phe_f_std,V_Phe_f_mean),
          #(S_Phe_f_mean-S_Phe_f_std,S_Phe_f_mean+S_Phe_f_std,S_Phe_f_mean),
          #(C_Phe_f_mean-C_Phe_f_std,C_Phe_f_mean+C_Phe_f_std,C_Phe_f_mean)
          ]  # Bounds for each parameter


# Generate Latin Hypercube Samples from uniform distribution as we don't know its exact distribution
samples = uniform_lhs(num_samples, bounds)
print(samples.shape,"****")
print(samples)


# Adjust samples to ensure their sum equals 1
#adjusted_samples = sum_test(samples, bounds)
print(samples.shape)


for i in range(num_parameters):
    figure = plt.figure()
    plt.hist(samples[:,i], density=True, alpha=0.6, color='grey')
    plt.xlabel('Value')
    plt.ylabel('Probability')
    plt.title('f_hpc - Histogram of Unifrom Distribution - Chemical fraction') #
    #mean - SD
    plt.axvline(bounds[i][0], color='k', linestyle='-.', label=f'lower = {bounds[i][0]:.3f}')
    #mean values
    plt.axvline((bounds[i][0]+bounds[i][1])/2, color='green', linestyle='--', label=f'$\mu$ = {(bounds[i][0]+bounds[i][1])/2:.3f}')
    #mean +SD
    plt.axvline(bounds[i][1], color='k', linestyle='--', label=f'upper = {bounds[i][1]:.3f}')   
    plt.legend(loc='lower center', bbox_to_anchor = (0.5, -0.3),ncol =3)
    plt.tight_layout()
    plt.show()
    #figure.savefig(f"output/analysis/Figures/4_Chemical_fraction_Histogram of Uniform Distribution_f_hpc.png")
    figure.savefig(f"output/analysis/Figures/4_Chemical_fraction_Histogram of Uniform Distribution_f_hpc.png")
    
#Chemical fraction
#Chemical fraction sum 
Chem_fraction_sum = prot_f_mean + hpc_f_mean + poly_f_mean + lip_f_mean + pig_f_mean + V_Phe_f_mean + S_Phe_f_mean + C_Phe_f_mean
Chem_fraction_sum_no_HPC = prot_f_mean + poly_f_mean + lip_f_mean + pig_f_mean + V_Phe_f_mean + S_Phe_f_mean + C_Phe_f_mean
print(Chem_fraction_sum)
print("HPC fraction",samples[:,0])
print("Rest of the chemical fraction",(Chem_fraction_sum - samples[:,0]))

#Protein
prot_f_samples = (Chem_fraction_sum - samples[:,0])*(prot_f_mean/Chem_fraction_sum_no_HPC)
print(prot_f_samples[0])
#Heteroplycondenstate
hpc_f_samples = samples[:,0]
#Polysaccharide
poly_f_samples = (Chem_fraction_sum - samples[:,0])*(poly_f_mean/Chem_fraction_sum_no_HPC)
#Lipids
lip_f_samples = (Chem_fraction_sum - samples[:,0])*(lip_f_mean/Chem_fraction_sum_no_HPC)
#Pigments
pig_f_samples = (Chem_fraction_sum - samples[:,0])*(pig_f_mean/Chem_fraction_sum_no_HPC)
#V-Phenol
V_Phe_f_samples = (Chem_fraction_sum - samples[:,0])*(V_Phe_f_mean/Chem_fraction_sum_no_HPC)
#S-Phenol
S_Phe_f_samples = (Chem_fraction_sum - samples[:,0])*(S_Phe_f_mean/Chem_fraction_sum_no_HPC)
#C-Phenol
C_Phe_f_samples = (Chem_fraction_sum - samples[:,0])*(C_Phe_f_mean/Chem_fraction_sum_no_HPC)

sum_check = prot_f_samples + hpc_f_samples + poly_f_samples + lip_f_samples + pig_f_samples + V_Phe_f_samples + S_Phe_f_samples + C_Phe_f_samples
print(sum_check)

Chem_sam = [prot_f_samples,poly_f_samples,lip_f_samples,pig_f_samples,V_Phe_f_samples,S_Phe_f_samples,C_Phe_f_samples]
Chem_frac = ['f_prot','f_poly','f_lip','f_pig','f_V_Phe','f_S_Phe','f_C_Phe']

bounds_1 = [(prot_f_mean-prot_f_std,prot_f_mean+prot_f_std,prot_f_mean),  
          #(hpc_f_mean-hpc_f_std, hpc_f_mean+hpc_f_std)
          (poly_f_mean-poly_f_std,poly_f_mean+poly_f_std,poly_f_mean),
          (lip_f_mean-lip_f_std,lip_f_mean+lip_f_std,lip_f_mean),
          (pig_f_mean-pig_f_std,pig_f_mean+pig_f_std,pig_f_mean),
          (V_Phe_f_mean-V_Phe_f_std,V_Phe_f_mean+V_Phe_f_std,V_Phe_f_mean),
          (S_Phe_f_mean-S_Phe_f_std,S_Phe_f_mean+S_Phe_f_std,S_Phe_f_mean),
          (C_Phe_f_mean-C_Phe_f_std,C_Phe_f_mean+C_Phe_f_std,C_Phe_f_mean)
          ]  # Bounds for each parameter

for i in range(len(Chem_sam)):
    figure = plt.figure()
    plt.hist(Chem_sam[i], density=True, alpha=0.6, color='grey')
    plt.xlabel('Value')
    plt.ylabel('Probability')
    plt.title(Chem_frac[i]+' - Histogram of Uniform Distribution - Chemical fraction')  
    #mean - SD
    plt.axvline(bounds_1[i][0], color='k', linestyle='-.', label=f'lower = {bounds_1[i][0]:.3f}')
    #mean values
    plt.axvline(bounds_1[i][2], color='green', linestyle='--', label=f'$\mu$ = {bounds_1[i][2]:.3f}')
    #mean +SD
    plt.axvline(bounds_1[i][1], color='k', linestyle='--', label=f'upper = {bounds_1[i][1]:.3f}')
    plt.legend(loc='lower center', bbox_to_anchor = (0.5, -0.3),ncol =3)
    plt.tight_layout()
    plt.show()
    figure.savefig(f"output/analysis/Figures/4_Chemical_fraction_Histogram of Uniform Distribution_{Chem_frac[i]}.png")

             
#CDOM component fractions
#Protein
prot_f_CDOM_samples = np.random.normal(prot_f_CDOM_mean, prot_f_CDOM_std, num_samples)
#Heteroplycondenstate
hpc_f_CDOM_samples = np.random.normal(hpc_f_CDOM_mean,hpc_f_CDOM_std, num_samples)
#Pigments
pig_f_CDOM_samples = np.random.normal(pig_f_CDOM_mean,pig_f_CDOM_std, num_samples)
#Porpherin
por_f_CDOM_samples = np.random.normal(pig_f_CDOM_mean,pig_f_CDOM_std, num_samples)
#V-Phenol
V_Phe_f_CDOM_samples = np.random.normal(V_Phe_f_CDOM_mean,V_Phe_f_CDOM_std, num_samples)
#S-Phenol
S_Phe_f_CDOM_samples = np.random.normal(S_Phe_f_CDOM_mean,S_Phe_f_CDOM_std, num_samples)
#C-Phenol
C_Phe_f_CDOM_samples = np.random.normal(C_Phe_f_CDOM_mean,C_Phe_f_CDOM_std, num_samples)


#                        ***   Chemical turnover time  ***
num_parameters = 14  # Number of parameters or dimensions
 # Bounds for each parameter: (min, max, mode)
bounds = [(prot_Tau_mean-prot_Tau_std,prot_Tau_mean+prot_Tau_std),
          (polypep_Tau_mean-polypep_Tau_std,polypep_Tau_mean+polypep_Tau_std),
          (amino_acid_Tau_mean-amino_acid_Tau_std,amino_acid_Tau_mean+amino_acid_Tau_std),
          (hpc_Tau_mean-hpc_Tau_std,hpc_Tau_mean+hpc_Tau_std),
          (poly_Tau_mean-poly_Tau_std,poly_Tau_mean+poly_Tau_std),
          (olig_Tau_mean-olig_Tau_std,olig_Tau_mean+olig_Tau_std),
          (humic_acid_Tau_mean-humic_acid_Tau_std,humic_acid_Tau_mean+humic_acid_Tau_std),
          (mono_Tau_mean-mono_Tau_std,mono_Tau_mean+mono_Tau_std),
          (ino_C_Tau_mean-ino_C_Tau_std,ino_C_Tau_mean+ino_C_Tau_std),
          (lip_Tau_mean-lip_Tau_std,lip_Tau_mean+lip_Tau_std),
          (pig_Tau_mean-pig_Tau_std,pig_Tau_mean+pig_Tau_std),
          (por_Tau_mean-por_Tau_std,por_Tau_mean+por_Tau_std),
          (V_Phe_Tau_mean-V_Phe_Tau_std,V_Phe_Tau_mean+V_Phe_Tau_std),
          (S_Phe_Tau_mean-S_Phe_Tau_std,S_Phe_Tau_mean+S_Phe_Tau_std),
          (C_Phe_Tau_mean-C_Phe_Tau_std,C_Phe_Tau_mean+C_Phe_Tau_std)]  # Bounds for each parameter


#print(bounds)
# Generate Latin Hypercube Samples from uniform distribution as we don't know its exact distribution
samples = uniform_lhs(num_samples, bounds)
print(samples.shape,"****")
#print(samples)

chem_tau = ['prot_Tau','polypep_Tau','amino_acid_Tau','hpc_Tau','poly_Tau','olig_Tau',
            'humic_acid_Tau','mono_Tau','ino_C_Tau','lip_Tau','pig_Tau','por_Tau',
            'V_Phe_Tau','S_Phe_Tau','C_Phe_Tau']

   
for i in range(len(chem_tau)):
    figure = plt.figure()
    plt.hist(samples[:,i], density=True, alpha=0.6, color='pink')
    plt.xlabel('Value')
    plt.ylabel('Probability')
    plt.title(chem_tau[i]+'Histogram of Uniform Distribution - Tau')
    #mean - SD
    plt.axvline(bounds[i][0], color='k', linestyle='-.', label=f'lower = {bounds[i][0]:.3f}')
    #mean values
    plt.axvline((bounds[i][0]+bounds[i][1])/2, color='green', linestyle='--', label=f'$\mu$ = {(bounds[i][0]+bounds[i][1])/2:.3f}')
    #mean +SD
    plt.axvline(bounds[i][1], color='k', linestyle='--', label=f'upper = {bounds[i][1]:.3f}')   
    plt.legend(loc='lower center', bbox_to_anchor = (0.5, -0.3),ncol =3)
    plt.tight_layout()
    plt.show()
    figure.savefig(f"output/analysis/Figures/5_Tau_Histogram of Unifrom Distribution_{chem_tau[i]}.png")
    
    
    
    
#Tau values - in seconds
#Protein
prot_Tau_samples = samples[:,0]
#Polypeptide
polypep_Tau_samples = samples[:,1]
#Amino acid
amino_acid_Tau_samples = samples[:,2]
#Heteroplycondenstate
hpc_Tau_samples = samples[:,3]
#Polysaccharide
poly_Tau_samples = samples[:,4]
#Oligosaccharide
olig_Tau_samples = samples[:,5]
#Humic acid
humic_acid_Tau_samples = samples[:,6]
#Monosaccharide
mono_Tau_samples = samples[:,7]
#Inorganic C
ino_C_Tau_samples = samples[:,8]
#Lipids
lip_Tau_samples = samples[:,9]
#Pigments
pig_Tau_samples = samples[:,10]
#Porphyrin
por_Tau_samples = samples[:,11]
#V-Phenol
V_Phe_Tau_samples = samples[:,12]
#S-Phenol
S_Phe_Tau_samples = samples[:,13]
#C-Phenol
C_Phe_Tau_samples = samples[:,14]
    
#Prod values - Initial Production Values
#Proteins
prot_Prod_samples = np.random.normal(prot_Prod_mean,prot_Prod_std, num_samples)
#Polysaccharide
poly_Prod_samples = np.random.normal(poly_Prod_mean,poly_Prod_std, num_samples)
#Lipids
lip_Prod_samples = np.random.normal(lip_Prod_mean,lip_Prod_std, num_samples)
#Pigments
pig_Prod_samples = np.random.normal(pig_Prod_mean,pig_Prod_std, num_samples)
#V-Phenol
V_Phe_Prod_samples = np.random.normal(V_Phe_Prod_mean,V_Phe_Prod_std, num_samples)
#S-Phenol
S_Phe_Prod_samples = np.random.normal(S_Phe_Prod_mean,S_Phe_Prod_std, num_samples)
#C-Phenol
C_Phe_Prod_samples = np.random.normal(C_Phe_Prod_mean,C_Phe_Prod_std, num_samples)

#                                      *** Creating CSV files     ***
#                                 -------------------------------------
#Creating CSV file/es for all initial input values for "n" no of sample runs

Velocity_Values = {
    'Teslin_V100 / ms-1':Tes_V100_samples,
    'Teslin_V1 / ms-1':Tes_V1_samples,

    'Pelly_V100 / ms-1':Pel_V100_samples,
    'Pelly_V1 / ms-1':Pel_V1_samples,

    'White+Donjec_V100 / ms-1':W_D_V100_samples,
    'White+Donjec_V1 / ms-1':W_D_V1_samples,

    'Stewart_V100 / ms-1':Ste_V100_samples,
    'Stewart_V1 / ms-1':Ste_V1_samples,

    'Porcupine_V100 / ms-1':Por_V100_samples,
    'Porcupine_V1 / ms-1':Por_V1_samples,

    'Tanana_V100 / ms-1':Tan_V100_samples,
    'Tanana_V1 / ms-1':Tan_V1_samples,

    'Koyukuk_V100 / ms-1':Koy_V100_samples,
    'Koyukuk_V1 / ms-1':Koy_V1_samples,
 
    'Yukon_V100 / ms-1':Yuk_V100_samples,
    'Yukon_V1 / ms-1':Yuk_V1_samples,
    'Yukon_V2 / ms-1':Yuk_V2_samples,
    'Yukon_V3 / ms-1':Yuk_V3_samples,
    'Yukon_V4 / ms-1':Yuk_V4_samples,
    'Yukon_V5 / ms-1':Yuk_V5_samples,
    'Yukon_V6 / ms-1':Yuk_V6_samples,
    'Yukon_V7 / ms-1':Yuk_V7_samples,
    'Yukon_V8 / ms-1':Yuk_V8_samples,  
    
          }
df = pd.DataFrame(Velocity_Values)
#print(df)
df.to_csv('output/random_initial_val/Velocity_Values.csv')
df.to_csv('output/analysis/Velocity_Values.csv')

Dilution_fractions = {
    'Teslin_Df':Tes_Df_samples,

    'Pelly_Df':Pel_Df_samples,

    'White+Donjec_Df':W_D_Df_samples,

    'Stewart_Df':Ste_Df_samples,

    'Porcupine_Df':Por_Df_samples,

    'Tanana_Df':Tan_Df_samples,

    'Koyukuk_Df':Koy_Df_samples,
          }
df = pd.DataFrame(Dilution_fractions)
#print(df)
df.to_csv('output/random_initial_val/Dilution_fractions.csv')
df.to_csv('output/analysis/Dilution_fractions.csv')

Initial_DOC_Values = {
    'Teslin_DOC / micromolar':Tes_DOC_samples,

    'Pelly_DOC / micromolar':Pel_DOC_samples,

    'White+Donjec_DOC/ micromolar':W_D_DOC_samples,

    'Stewart_DOC / micromolar':Ste_DOC_samples,

    'Porcupine_DOC / micromolar':Por_DOC_samples,

    'Tanana_DOC / micromolar':Tan_DOC_samples,

    'Koyukuk_DOC / micromolar':Koy_DOC_samples,
          }
df = pd.DataFrame(Initial_DOC_Values)
#print(df)
df.to_csv('output/random_initial_val/Initial_DOC_Values.csv')
df.to_csv('output/analysis/Initial_DOC_Values.csv')

Chemical_fraction = {
    'Protein':prot_f_samples,

    'Polypeptide':"",

    'Amino acid':"",

    'Heteroploycondensate(HPC)':hpc_f_samples,

    'Polysaccharides':poly_f_samples,

    'Oligosaccharides':"",

    'Humic acid':"",

    'Monosaccharides':"",

    'Inorganic carbon':"",

    'Lipids':lip_f_samples,

    'Pigments':pig_f_samples,

    'Porphyrin':"",

    'V-Phenol':V_Phe_f_samples,    

    'S-Phenol':S_Phe_f_samples,   

    'C-Phenol':C_Phe_f_samples,   
          }
df = pd.DataFrame(Chemical_fraction)
#print(df)
df.to_csv('output/random_initial_val/Chemical_fraction.csv')
df.to_csv('output/analysis/Chemical_fraction.csv')

CDOM_component_fractions = {
    'Protein':prot_f_CDOM_samples,

    'Polypeptide':"",

    'Amino acid':"",

    'Heteroploycondensate(HPC)':hpc_f_CDOM_samples,

    'Polysaccharides':"",

    'Oligosaccharides':"",

    'Humic acid':"",

    'Monosaccharides':"",

    'Inorganic carbon':"",

    'Lipids':"",

    'Pigments':pig_f_CDOM_samples,

    'Porphyrin':por_f_CDOM_samples,

    'V-Phenol':V_Phe_f_CDOM_samples,    

    'S-Phenol':S_Phe_f_CDOM_samples,   

    'C-Phenol':C_Phe_f_CDOM_samples,   
          }
df = pd.DataFrame(CDOM_component_fractions)
#print(df)
df.to_csv('output/random_initial_val/CDOM_component_fractions.csv')
df.to_csv('output/analysis/CDOM_component_fractions.csv')

Tau_Values = {
    'Protein / seconds':prot_Tau_samples,

    'Polypeptide / seconds':polypep_Tau_samples,

    'Amino acid / seconds':amino_acid_Tau_samples,

    'Heteroploycondensate(HPC) / seconds':hpc_Tau_samples,

    'Polysaccharides / seconds':poly_Tau_samples,

    'Oligosaccharides / seconds':olig_Tau_samples,

    'Humic acid / seconds':humic_acid_Tau_samples,

    'Monosaccharides / seconds':mono_Tau_samples,

    'Inorganic carbon/ seconds':ino_C_Tau_samples,

    'Lipids / seconds':lip_Tau_samples,

    'Pigments / seconds':pig_Tau_samples,

    'Porphyrin / seconds':por_Tau_samples,

    'V-Phenol / seconds':V_Phe_Tau_samples,    

    'S-Phenol / seconds':S_Phe_Tau_samples,   

    'C-Phenol / seconds':C_Phe_Tau_samples,   
          }
df = pd.DataFrame(Tau_Values)
#print(df)
df.to_csv('output/random_initial_val/Tau_Values.csv')
df.to_csv('output/analysis/Tau_Values.csv')

Production_Values = {
    'Protein':prot_Prod_samples,

    'Polypeptide':"",

    'Amino acid':"",

    'Heteroploycondensate(HPC)':"",

    'Polysaccharides':poly_Prod_samples,

    'Oligosaccharides':"",

    'Humic acid':"",

    'Monosaccharides':"",

    'Inorganic carbon':"",

    'Lipids':lip_Prod_samples,

    'Pigments':pig_Prod_samples,

    'Porphyrin':"",

    'V-Phenol':V_Phe_Prod_samples,    

    'S-Phenol':S_Phe_Prod_samples,   

    'C-Phenol':C_Phe_Prod_samples,   
          }
df = pd.DataFrame(Production_Values)
#print(df)
df.to_csv('output/random_initial_val/Production_Values.csv')
df.to_csv('output/analysis/Production_Values.csv')

