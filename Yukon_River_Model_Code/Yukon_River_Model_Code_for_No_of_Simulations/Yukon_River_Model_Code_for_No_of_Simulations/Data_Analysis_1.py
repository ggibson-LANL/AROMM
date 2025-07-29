import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import river_initial_monac_code as rv
import Yukon_river as yr

#Chemicals
yr.Conc_Che0       #Protein 
yr.Conc_Che1       #Polypeptide
yr.Conc_Che2         #Amino acid
yr.Conc_Che3         #Heteropolycondensate
yr.Conc_Che4         #Polysaccharide
yr.Conc_Che5         #Oligosaccharide
yr.Conc_Che6         #Humic acid
yr.Conc_Che7         #Monosaccharides
yr.Conc_Che8         #Inorganic carbon
yr.Conc_Che9         #Lipids
yr.Conc_Che10        #Pigments
yr.Conc_Che11        #Porphyrin
yr.Conc_Che12        #CDOM
yr.Conc_Che13        #V-Phenol
yr.Conc_Che14        #S-Phenol
yr.Conc_Che15        #C-Phenol
yr.Conc_Che16        #Phenol
yr.Conc_Che17       #DOC

#Velocities
#Teslin
rv.Tes_V100_samples 
rv.Tes_V1_samples 
#Pelly
rv.Pel_V100_samples 
rv.Pel_V1_samples 
#White+Donjec
rv.W_D_V100_samples 
rv.W_D_V1_samples 
#Stewart  
rv.Ste_V100_samples 
rv.Ste_V1_samples 
#Porcupine 
rv.Por_V100_samples 
rv.Por_V1_samples 
#Tanana 
rv.Tan_V100_samples 
rv.Tan_V1_samples 
#Koyukuk
rv.Koy_V100_samples 
rv.Koy_V1_samples 
#Yukon - Yuk
rv.Yuk_V100_samples 
rv.Yuk_V1_samples 
rv.Yuk_V2_samples 
rv.Yuk_V3_samples 
rv.Yuk_V4_samples 
rv.Yuk_V5_samples 
rv.Yuk_V6_samples 
rv.Yuk_V7_samples 
rv.Yuk_V8_samples 

#Dilution fractions
#Teslin
rv.Tes_Df_samples 
#Pelly
rv.Pel_Df_samples 
#White+Donjec
rv.W_D_Df_samples 
#Stewart  
rv.Ste_Df_samples 
#Porcupine - Por
rv.Por_Df_samples 
#Tanana - Tan 
rv.Tan_Df_samples 
#Koyukuk - Koy
rv.Koy_Df_samples

#Initial DOC concentrations
#Teslin
rv.Tes_DOC_samples 
#Pelly
rv.Pel_DOC_samples 
#White+Donjec
rv.W_D_DOC_samples
#Stewart  
rv.Ste_DOC_samples 
#Porcupine - Por
rv.Por_DOC_samples 
#Tanana - Tan 
rv.Tan_DOC_samples
#Koyukuk - Koy
rv.Koy_DOC_samples 


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


#Tau values - in seconds
#Protein
rv.prot_Tau_samples
#Polypeptide
rv.polypep_Tau_samples 
#Amino acid
rv.amino_acid_Tau_samples 
#Heteroplycondenstate
rv.hpc_Tau_samples
#Polysaccharide
rv.poly_Tau_samples
#Oligosaccharide
rv.olig_Tau_samples
#Humic acid
rv.humic_acid_Tau_samples 
#Monosaccharide
rv.mono_Tau_samples
#Inorganic C
rv.ino_C_Tau_samples
#Lipids
rv.lip_Tau_samples
#Pigments
rv.pig_Tau_samples
#Porphyrin
rv.por_Tau_samples
#V-Phenol
rv.V_Phe_Tau_samples
#S-Phenol
rv.S_Phe_Tau_samples
#C-Phenol
rv.C_Phe_Tau_samples

#Prod values - Initial Production Values
#Proteins
rv.prot_Prod_samples 
#Polysaccharide
rv.poly_Prod_samples 
#Lipids
rv.lip_Prod_samples 
#Pigments
rv.pig_Prod_samples 
#V-Phenol
rv.V_Phe_Prod_samples 
#S-Phenol
rv.S_Phe_Prod_samples 
#C-Phenol
rv.C_Phe_Prod_samples 





# ---------------------------------    Creating plots ----------------------------------------------

print(Conc_Che17) 
            
#Chemical conc vs distance
n  = str(i)
plt.plot(dist_,C_Conc[0:j,3] ,'k--',label="DOC - "+n+ " th sample")
plt.scatter([400,800],[200,700], color='k', marker='o', label="DOC - " + str(j) + "th sample")   
plt.legend(loc='upper right', bbox_to_anchor=(0.9, 0.95), fontsize='medium')
plt.grid(True)
plt.ylabel("[DOC] / micromolar")
plt.xlabel("Distance / km")
plt.title("Teslin River")
plt.show()
#Chemical conc vs time
plt.plot(time_,C_Conc[0:j,17] ,'r--',label="DOC - "+n+ " th sample")
plt.legend(loc='upper right', bbox_to_anchor=(0.9, 0.95), fontsize='medium')
plt.grid(True)
plt.ylabel("[DOC] / micromolar")
plt.xlabel("Time / days")
plt.title("Teslin River")
plt.show()


