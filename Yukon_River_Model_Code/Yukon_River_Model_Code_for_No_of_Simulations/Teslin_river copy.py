import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import river_initial_monac_code as rv

#                               ***    function of the chemical decay   ***
#                                 -------------------------------------
def f1(dt,C0,Tau,P):
    term1 = C0 * math.exp((-1*dt)/Tau)
    term2 = P * Tau * (1-math.exp((-1*dt)/Tau))
    #print("C0",C0)
    #print(Tau)
    #print(dt/Tau)
    #print(math.exp((-1*dt)/Tau))
    #print("Decay value", term1)
    #print("Production term", term2)
    Conc = term1 + term2
    # print("The concentration is " , Conc)
    return Conc

#                               ***  loop for the monte-carlo analysis   ***
#                                   -------------------------------------

Tes_imax =[]
Tes_tmax = []
results =[]
for i in range(rv.num_samples):
    
                    # Code to generate the imax and tmax values #
                    
     Tes_velo = [rv.Tes_V100_samples[i],rv.Tes_V1_samples[i]]
     #print(Tes_velo)
     
     Tes_bounds =[]
     for j in range(len(rv.Tes_nodes)-1):
         #print(j)
         a = rv.Tes_nodes[j]/(Tes_velo[j]*rv.dt*(1/1000))
         b= rv.Tes_nodes[j+1]/(Tes_velo[j]*rv.dt*(1/1000))
         Tes_bounds.append(a)
         Tes_bounds.append(b)
     #print(Tes_bounds) 
     #print(len(Tes_bounds))
     
     Tes_ts_val = [] # np.arange() values - all the time steps
     Tes_i_val =[] # no of timesetps in each segments
     x=0
     while x < len(Tes_bounds):
         a = np.arange(Tes_bounds[x],Tes_bounds[x+1])
         Tes_ts_val.append(a)
         b =len(a) 
         Tes_i_val.append(b)
         x = x+2
    #print(len(ts_val))
    #print(len(i_val))
     print("##Tes_i_values##",Tes_i_val)
     
    #The most important part
     Tes_i_sum=[]
    #values for the while loop
     for y in range(len(Tes_i_val)):
         #print(i)
         sum_ = np.sum(Tes_i_val[0:y+1])
         Tes_i_sum.append(sum_)
     print("Values need to run the river code",Tes_i_sum) #"for the last i value,additional +1 has added")
      #print("i_values",len(i_values))
    
     #Teslin_imax values = the last item of the Tes_i_sum
     Tes_imax.append(Tes_i_sum[len(Tes_i_sum)-1])
     #print(Tes_imax)
     
     #time_calculation
     #print("Tes_imax_values:",Tes_imax[i])
     t_ts = np.arange(0,Tes_imax[i])*rv.dt*(1/86400) # time is in days
     time_ = t_ts.tolist()
     #print("time_length",len(time_))
     Tes_tmax.append(time_[len(time_)-1])
     #print(Tes_tmax)

     imax = Tes_imax[i]
     #print(imax)
                         # Creating arrays for C_Conc,Tau and Prod  #
     
     C_Conc = np.ones((imax,rv.cmax)) 
     Tau = np.ones((imax,rv.cmax))
     Prod = np.ones((imax,rv.cmax))
     
            #  Getting initial values for C_Conc, Tau and Prod for the chemical  #
                      
     #Initial concentration of the species ia an array - units micrmolar C (micromol C / L)
     DOC=rv.Tes_DOC_samples[i]
     #Protein
     C_Conc[:,0]=DOC*rv.prot_f_samples[i]
     #Poltpeptide
     C_Conc[:,1]=0
     #Amino acid
     C_Conc[:,2]=0
     #Heteroploycondensate
     C_Conc[:,3]=DOC*rv.hpc_f_samples[i]
     #Polysaccharides
     C_Conc[:,4]=DOC*rv.poly_f_samples[i]
     #Oligosaccharides
     C_Conc[:,5]=0
     #Humic acid
     C_Conc[:,6]=0
     #Monosaccharides
     C_Conc[:,7]=0
     #Inorganic carbon
     C_Conc[:,8]=0
     #Lipids
     C_Conc[:,9]=DOC*rv.lip_f_samples[i]
     #Pigments
     C_Conc[:,10]=DOC*rv.pig_f_samples[i]
     #Porpherin
     C_Conc[:,11]=0
     #VPhenol
     C_Conc[:,13]=DOC*rv.V_Phe_f_samples[i]
     #SPhenol
     C_Conc[:,14]=DOC*rv.S_Phe_f_samples[i]
     #CPhenol
     C_Conc[:,15]=DOC*rv.C_Phe_f_samples[i]
     #CDOM
     C_Conc[:,12] = (rv.prot_f_CDOM_samples[i] * C_Conc[:,0]) + \
               (rv.hpc_f_CDOM_samples[i] * C_Conc[:,3]) + \
               (rv.pig_f_CDOM_samples[i] * C_Conc[:,10]) + \
               (rv.por_f_CDOM_samples[i] * C_Conc[:,11]) + \
               (rv.V_Phe_f_CDOM_samples[i] * C_Conc[:,13]) + \
               (rv.S_Phe_f_CDOM_samples * C_Conc[:,14]) + \
               (rv.C_Phe_f_CDOM_samples * C_Conc[:,15])
     #Phenol
     C_Conc[:,16] = C_Conc[:,13]+C_Conc[:,14]+C_Conc[:,15]      
     #TDOC # Does not have Inorganic Carbon-8
     TDOC = (C_Conc[:,0] + C_Conc[:,1] + C_Conc[:,2] + C_Conc[:,3] +
             C_Conc[:,4] + C_Conc[:,5] + C_Conc[:,6] + C_Conc[:,7] +
             C_Conc[:,9] + C_Conc[:,10] + C_Conc[:,11] + C_Conc[:,13] +
             C_Conc[:,14] + C_Conc[:,15])
       
     #Reaction time (Tau) values of chemical species in an array  
     #Protein
     Tau[:,0] =rv.prot_Tau_samples[i]
     #Polypeptide
     Tau[:,1] =rv.polypep_Tau_samples[i]
     #Amino acid
     Tau[:,2] =rv.amino_acid_Tau_samples[i]
     #Heteroploycondensate(HPC)
     Tau[:,3] =rv. hpc_Tau_samples[i]
     #Polysaccharides
     Tau[:,4] =rv.poly_Tau_samples[i]
     #Oligosaccharides
     Tau[:,5] =rv.olig_Tau_samples[i]
     #Humic acid
     Tau[:,6] =rv.humic_acid_Tau_samples[i]
     #Monosaccharides
     Tau[:,7] =rv.mono_Tau_samples[i]
     #Inorganic carbon
     Tau[:,8] =rv.ino_C_Tau_samples[i]
     #Lipids
     Tau[:,9] =rv.lip_Tau_samples[i]
     #Pigments
     Tau[:,10] =rv.pig_Tau_samples[i]
     #Porphyrin
     Tau[:,11] =rv.por_Tau_samples[i]
     #CDOM
     #Tau[:,12] = -
     #V-Phenol
     Tau[:,13] =rv.V_Phe_Tau_samples[i]
     #S-Phenol
     Tau[:,14] =rv.S_Phe_Tau_samples[i]
     #C-Phenol
     Tau[:,15] =rv.C_Phe_Tau_samples[i]
     
     #Initial Production (Prod) values of chemical species in an array 
     #Protein
     Prod[:,0] = rv.prot_Prod_samples
     #Polypeptide
     #Prod[0,1] = -
     #Amino acid
     #Prod[:,2] = -
     #Heteroploycondensate(HPC)
     #Prod[:,3] = -
     #Polysaccharides
     Prod[:,4] = rv.poly_Prod_samples
     #Oligosaccharides
     #Prod[:,5] = -
     #Humic acid
     #Prod[:,6] = -
     #Monosaccharides
     #Prod[:,7] = -
     #Inorganic carbon
     #Prod[:,8] = -
     #Lipids
     Prod[:,9] = rv.lip_Prod_samples
     #Pigments
     Prod[:,10] = rv.pig_Prod_samples
     #Porphyrin
     #Prod[:,11] = -
     #CDOM
     #Prod[:,12] = -
     #V-Phenol
     Prod[:,13] = rv.V_Phe_Prod_samples
     #S-Phenol
     Prod[:,14] = rv.S_Phe_Prod_samples
     #C-Phenol
     Prod[:,15] = rv.C_Phe_Prod_samples
    
                                             #  Entering to the river code  #
     j=0
     while j < Tes_imax[i]:
     #while j < imax:
        #Protein 
        Conc0 = f1(dt, C_Conc[j,0], Tau[j,0], Prod[j,0])
        #Chem1:Peptides, two chanells
        Prod[j,1] = C_Conc[j,0]/Tau[j,0]
        Conc1 =  f1(dt,C_Conc[j,1],Tau[j,1],Prod[j,1])
        # Run the rest of your simulation code for this set of parameters
        if j < imax - 1:
            C_Conc[j+1,0]=Conc0
            C_Conc[j+1,1]=Conc1
        j=j+1
     #print("j",j)
     #C_Conc[j,0] =((D1_samples[i])*C_Conc[j,0])+((D2_samples[i])*tr_C_Conc_samples[i])
     #while j < imax:
     #   Conc0 = f1(dt, C_Conc[j,0], Tau[j,0], Prod[j,0])
        # Run the rest of your simulation code for this set of parameters
     #   if j < imax - 1:
      #      C_Conc[j+1,0]=Conc0
      #  j=j+1        
     #print(j)
     #print(Conc0)
     #print(C_Conc[j-1,0])
     # Store the result 
     results.append(Conc0)

#print(time_) 
#print(len(time_))
#print(len(time_))

#print("Teslin_imax",Tes_imax)
#print("Teslin_tmax",Tes_tmax)
#print(results)
#print(len(results))

#Creating a CSV file for values of the chemcial concentration
Sample_Values = {
    'Teslin_V100 / ms-1':Tes_V100_samples,
    'Teslin_V1 / ms-1':Tes_V1_samples,
    'Teslin_imax':Tes_imax,
    'Teslin_tmax / days':Tes_tmax,
    'Teslin_DOC':Tes_DOC_samples,
    'Protein_fraction':prot_f_samples,
    'Protein_Tau_samples':prot_Tau_samples,
    'Protein_Prod_samples / days':prot_Prod_samples,
    'Polypeptide_Tau_samples':polypep_Tau_samples,
    'Output_Conc':results,
    
          }
df = pd.DataFrame(Sample_Values)
#print(df)
df.to_csv('output_m/Sample values.csv')
    
# Example data
#x =DOC_samples 
#print(len(x))
#x1=Tau_samples
#print(len(x1))
#x2=Prod_samples
#print(len(x2))
#x3 = Tes_V100_samples 
#x4 = Tes_V1_samples
#x3=imax_int
#x4=n1_samples
#x5=D1_samples
#x6=D2_samples
#x7=tr_C_Conc_samples
y = results

# Create scatter plot
#plt.scatter(x, y)
plt.title('Scatter Plot Example')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.show

#plt.scatter(x1, y)
plt.title('Scatter Plot Example')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.show

#plt.scatter(x2, y)
plt.title('Scatter Plot Example')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.show

#plt.scatter(x3, y)
plt.title('Scatter Plot Example')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.show

#plt.scatter(x4, y)
plt.title('Scatter Plot Example')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.show
# Create scatter plot
#plt.scatter(x1, y)
# Create scatter plot
#plt.scatter(x2, y)
# Create scatter plot
#plt.scatter(x3, y)
# Create scatter plot
#plt.scatter(x4, y)
# Create scatter plot
#plt.scatter(x5, y)
# Create scatter plot
#plt.scatter(x6, y)
# Create scatter plot
#plt.scatter(x7, y)
