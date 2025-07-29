import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

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
raws = 5000


#                        ***    Initial Paramters and standard deviations    ***
#                                 -------------------------------------

#      Rivers
#River name acronyms:            Nodes:Dist       Velocities         Ini.DOC
#    Teslin - Tes                    x                x                x
#    Pelly - Pel                     x                x                x
#    White+Donjec - W_D              x                x                x
#    Stewart - Ste                   x                x                x
#    Porcupine - Por                 x                x                x
#    Tanana - Tan                    x                x                x
#    Koyukuk - Koy                   x                x                x
#    Yukon - Yuk                     x                x             (Teslin)   
#Rivers variables:
#    nodes - Distances 
#    velocities
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

#Velocities: in ms-1
#V100 for the 1st 100 km 
#Teslin
Tes_V100_mean = 0.1
Tes_V100_std = 0.000
Tes_V1_mean = 0.4
Tes_V1_std = 0.04
#Pelly
Pel_V100_mean = 0.1
Pel_V100_std = 0.000
Pel_V1_mean = 0.4
Pel_V1_std = 0.17
#White+Donjec
W_D_V100_mean = 0.1
W_D_V100_std = 0.000
W_D_V1_mean = 0.4
W_D_V1_std = 0.17
#Stewart
Ste_V100_mean = 0.1
Ste_V100_std = 0.000
Ste_V1_mean = 0.4
Ste_V1_std = 0.17
#Porcupine
Por_V100_mean = 0.1
Por_V100_std = 0.000
Por_V1_mean = 0.4
Por_V1_std = 0.17
#Tanana
Tan_V100_mean = 0.1
Tan_V100_std = 0.000
Tan_V1_mean = 0.4
Tan_V1_std = 0.17
#Koyukuk
Koy_V100_mean = 0.1
Koy_V100_std = 0.000
Koy_V1_mean = 0.4
Koy_V1_std = 0.17
#Yukon
Yuk_V100_mean = 0.1   #first 100 km
Yuk_V100_std = 0.000
Yuk_V1_mean = 0.4     #from 100 km to Pelly R meeting point(CP1) 
Yuk_V1_std = 0.17
Yuk_V2_mean = 0.4     #from CP1 to White+D R meeting point(CP2) 
Yuk_V2_std = 0.17
Yuk_V3_mean = 0.4     #from CP2 to Stewart R meeting point(CP3) 
Yuk_V3_std = 0.17
Yuk_V4_mean = 0.4     #from CP3 to Porcupine R meeting point(CP4) 
Yuk_V4_std = 0.17
Yuk_V5_mean = 0.4     #from CP4 to Tanana R meeting point(CP5) 
Yuk_V5_std = 0.17
Yuk_V6_mean = 0.4     #from CP5 to Koyukuk R meeting point(CP6) 
Yuk_V6_std = 0.17
Yuk_V7_mean = 0.4     #from CP6 to Delta region starting point(DR) 
Yuk_V7_std = 0.17
Yuk_V8_mean = 0.4     #from DR to river end point (Delta region)
Yuk_V8_std = 0.17

#Initial DOC concentration: in micromolar C (micromol C / L)
#Teslin
Tes_DOC_mean = 600
Tes_DOC_std = 30
#Pelly   
Pel_DOC_mean = 600
Pel_DOC_std = 30
#White+Donjec  
W_D_DOC_mean = 600
W_D_DOC_std = 30
#Stewart  
Ste_DOC_mean = 600
Ste_DOC_std = 30 
#Porcupine  
Por_DOC_mean = 600
Por_DOC_std = 30
#Tanana  
Tan_DOC_mean = 600
Tan_DOC_std = 30
#Koyukuk   
Koy_DOC_mean = 600
Koy_DOC_std = 30

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
prot_f_mean = 5/100
prot_f_std = 0.00
#Heteroploycondensate(HPC)-3
hpc_f_mean = 85/100
hpc_f_std = 0.00
#Polysaccharides-4
poly_f_mean = 5/100
poly_f_std = 0.00
#Lipids-9
lip_f_mean = 5/100
lip_f_std = 0.00
#Pigments-10
pig_f_mean = 1/100
pig_f_std = 0.00
#VPhenol-13
V_Phe_f_mean = 0.1/100
V_Phe_f_std = 0.00
#SPhenol-14
S_Phe_f_mean = 0.2/100
S_Phe_f_std = 0.00
#CPhenol-15
C_Phe_f_mean = 0.3/100
C_Phe_f_std = 0.00

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
#protein-0
prot_Tau_mean =8640000
prot_Tau_std = 300
#Polypeptide-1
polypep_Tau_mean = 864000*(3/13)*(10)
polypep_Tau_std = 0
#Amino acid-2
amino_acid_Tau_mean = 864000*(3/13)*(10)
amino_acid_Tau_std = 0
#Heteroploycondensate(HPC)-3
hpc_Tau_mean = 864000*(3/13)*(10)
hpc_Tau_std = 0
#Polysaccharides-4
poly_Tau_mean = 864000*(3/13)*(10)
poly_Tau_std = 0
#Oligosaccharides-5
olig_Tau_mean = 864000*(3/13)*(10)
olig_Tau_std = 0
#Humic acid-6
humic_acid_Tau_mean = 864000*(3/13)*(10)
humic_acid_Tau_std = 0
#Monosaccharides-7
mono_Tau_mean = 864000*(3/13)*(10)
mono_Tau_std = 0
#Inorganic carbon-8
ino_C_Tau_mean = 864000*(3/13)*(10)
ino_C_Tau_std = 0
#Lipids-9
lip_Tau_mean = 864000*(3/13)*(10)
lip_Tau_std = 0
#Pigment-10
pig_Tau_mean = 864000*(3/13)*(10)
pig_Tau_std = 0
#Porphyrin-11
por_Tau_mean = 864000*(3/13)*(10)
por_Tau_std = 0
#VPhenol-13
V_Phe_Tau_mean = 864000*(3/13)*(10)
V_Phe_Tau_std = 0
#SPhenol-14
S_Phe_Tau_mean = 864000*(3/13)*(10)
S_Phe_Tau_std = 0
#CPhenol-15
C_Phe_Tau_mean = 864000*(3/13)*(10)
C_Phe_Tau_std = 0

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

#Velocities
#Teslin
Tes_V100_samples = np.random.normal(Tes_V100_mean, Tes_V100_std, num_samples)
Tes_V1_samples = np.random.normal(Tes_V1_mean, Tes_V1_std, num_samples)
#Pelly
Pel_V100_samples = np.random.normal(Pel_V100_mean, Pel_V100_std, num_samples)
Pel_V1_samples = np.random.normal(Pel_V1_mean, Pel_V1_std, num_samples)
#White+Donjec
W_D_V100_samples = np.random.normal(W_D_V100_mean, W_D_V100_std, num_samples)
W_D_V1_samples = np.random.normal(W_D_V1_mean, W_D_V1_std, num_samples)
#Stewart  
Ste_V100_samples = np.random.normal(Ste_V100_mean, Ste_V100_std, num_samples)
Ste_V1_samples = np.random.normal(Ste_V1_mean, Ste_V1_std, num_samples)
#Porcupine 
Por_V100_samples = np.random.normal(Por_V100_mean, Por_V100_std, num_samples)
Por_V1_samples = np.random.normal(Por_V1_mean, Por_V1_std, num_samples)
#Tanana 
Tan_V100_samples = np.random.normal(Tan_V100_mean, Tan_V100_std, num_samples)
Tan_V1_samples = np.random.normal(Tan_V1_mean, Tan_V1_std, num_samples)
#Koyukuk
Koy_V100_samples = np.random.normal(Koy_V100_mean, Koy_V100_std, num_samples)
Koy_V1_samples = np.random.normal(Koy_V1_mean, Koy_V1_std, num_samples)
#Yukon - Yuk
Yuk_V100_samples = np.random.normal(Yuk_V100_mean, Yuk_V100_std, num_samples)
Yuk_V1_samples = np.random.normal(Yuk_V1_mean, Yuk_V1_std, num_samples)
Yuk_V2_samples = np.random.normal(Yuk_V2_mean, Yuk_V2_std, num_samples)
Yuk_V3_samples = np.random.normal(Yuk_V3_mean, Yuk_V3_std, num_samples)
Yuk_V4_samples = np.random.normal(Yuk_V4_mean, Yuk_V4_std, num_samples)
Yuk_V5_samples = np.random.normal(Yuk_V5_mean, Yuk_V5_std, num_samples)
Yuk_V6_samples = np.random.normal(Yuk_V6_mean, Yuk_V6_std, num_samples)
Yuk_V7_samples = np.random.normal(Yuk_V7_mean, Yuk_V7_std, num_samples)
Yuk_V8_samples = np.random.normal(Yuk_V8_mean, Yuk_V8_std, num_samples)

#Initial DOC concentrations
#Teslin
Tes_DOC_samples = np.random.normal(Tes_DOC_mean, Tes_DOC_std, num_samples)
#Pelly
Pel_DOC_samples = np.random.normal(Pel_DOC_mean, Pel_DOC_std, num_samples)
#White+Donjec
W_D_DOC_samples = np.random.normal(W_D_DOC_mean, W_D_DOC_std, num_samples)
#Stewart  
Ste_DOC_samples = np.random.normal(Ste_DOC_mean, Ste_DOC_std, num_samples)
#Porcupine - Por
Por_DOC_samples = np.random.normal(Por_DOC_mean, Por_DOC_std, num_samples)
#Tanana - Tan 
Tan_DOC_samples = np.random.normal(Tan_DOC_mean, Tan_DOC_std, num_samples)
#Koyukuk - Koy
Koy_DOC_samples = np.random.normal(Koy_DOC_mean, Koy_DOC_std, num_samples)


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

#Chemical fraction
#Protein
prot_f_samples = np.random.normal(prot_f_mean, prot_f_std, num_samples)
#Heteroplycondenstate
hpc_f_samples = np.random.normal(hpc_f_mean,hpc_f_std, num_samples)
#Polysaccharide
poly_f_samples = np.random.normal(poly_f_mean,poly_f_std, num_samples)
#Lipids
lip_f_samples = np.random.normal(lip_f_mean,lip_f_std, num_samples)
#Pigments
pig_f_samples = np.random.normal(pig_f_mean,pig_f_std, num_samples)
#V-Phenol
V_Phe_f_samples = np.random.normal(V_Phe_f_mean,V_Phe_f_std, num_samples)
#S-Phenol
S_Phe_f_samples = np.random.normal(S_Phe_f_mean,S_Phe_f_std, num_samples)
#C-Phenol
C_Phe_f_samples = np.random.normal(C_Phe_f_mean,C_Phe_f_std, num_samples)

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

#Tau values - in seconds
#Protein
prot_Tau_samples = np.random.normal(prot_Tau_mean,prot_Tau_std, num_samples)
#Polypeptide
polypep_Tau_samples = np.random.normal(polypep_Tau_mean,polypep_Tau_std, num_samples)
#Amino acid
amino_acid_Tau_samples = np.random.normal(amino_acid_Tau_mean,amino_acid_Tau_std, num_samples)
#Heteroplycondenstate
hpc_Tau_samples = np.random.normal(hpc_Tau_mean,hpc_Tau_std, num_samples)
#Polysaccharide
poly_Tau_samples = np.random.normal(poly_Tau_mean,poly_Tau_std, num_samples)
#Oligosaccharide
olig_Tau_samples = np.random.normal(olig_Tau_mean,olig_Tau_std, num_samples)
#Humic acid
humic_acid_Tau_samples = np.random.normal(humic_acid_Tau_mean,humic_acid_Tau_std, num_samples)
#Monosaccharide
mono_Tau_samples = np.random.normal(mono_Tau_mean,mono_Tau_std, num_samples)
#Inorganic C
ino_C_Tau_samples = np.random.normal(ino_C_Tau_mean,ino_C_Tau_std, num_samples)
#Lipids
lip_Tau_samples = np.random.normal(lip_Tau_mean,lip_Tau_std, num_samples)
#Pigments
pig_Tau_samples = np.random.normal(pig_Tau_mean,pig_Tau_std, num_samples)
#Porphyrin
por_Tau_samples = np.random.normal(por_Tau_mean,por_Tau_std, num_samples)
#V-Phenol
V_Phe_Tau_samples = np.random.normal(V_Phe_Tau_mean,V_Phe_Tau_std, num_samples)
#S-Phenol
S_Phe_Tau_samples = np.random.normal(S_Phe_Tau_mean,S_Phe_Tau_std, num_samples)
#C-Phenol
C_Phe_Tau_samples = np.random.normal(C_Phe_Tau_mean,C_Phe_Tau_std, num_samples)

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







