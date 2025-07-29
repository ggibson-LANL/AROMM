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
#Creating lists for the end values
imax =[]  #Total number of time steps at the end of each iteration cycle
tmax = [] #Total time to reach the coast at the end of each iteration cycle - in days
dmax = [] #Total river distance at the end of each iteration cycle.
#Chemicals
Conc_Che0 = []        #Protein 
Conc_Che1 = []        #Polypeptide
Conc_Che2 = []        #Amino acid
Conc_Che3 = []        #Heteropolycondensate
Conc_Che4 = []        #Polysaccharide
Conc_Che5 = []        #Oligosaccharide
Conc_Che6 = []        #Humic acid
Conc_Che7 = []        #Monosaccharides
Conc_Che8 = []        #Inorganic carbon
Conc_Che9 = []        #Lipids
Conc_Che10 = []       #Pigments
Conc_Che11 = []       #Porphyrin
Conc_Che12 = []       #CDOM
Conc_Che13 = []       #V-Phenol
Conc_Che14 = []       #S-Phenol
Conc_Che15 = []       #C-Phenol
Conc_Che16 = []       #Phenol
Conc_Che17 = []       #DOC

#Creating 2D arrays for store all iteration values for n no of samples
i_iter = np.full((rv.raws, rv.num_samples), np.nan)  # i iteration timesteps per cycles for no of samples 
t_iter = np.full((rv.raws, rv.num_samples), np.nan)  # time for iteration steps per cycles for no of samples
d_iter = np.full((rv.raws, rv.num_samples), np.nan)  # distance for iteration steps per cycles for no of samples
#Chemicals # chemical concentration values per iteration cycle for no of samples
Che0_iter = np.full((rv.raws, rv.num_samples), np.nan)
Che1_iter = np.full((rv.raws, rv.num_samples), np.nan)
Che2_iter = np.full((rv.raws, rv.num_samples), np.nan)
Che3_iter = np.full((rv.raws, rv.num_samples), np.nan)
Che4_iter = np.full((rv.raws, rv.num_samples), np.nan)
Che5_iter = np.full((rv.raws, rv.num_samples), np.nan)
Che6_iter = np.full((rv.raws, rv.num_samples), np.nan)
Che7_iter = np.full((rv.raws, rv.num_samples), np.nan)
Che8_iter = np.full((rv.raws, rv.num_samples), np.nan)
Che9_iter = np.full((rv.raws, rv.num_samples), np.nan)
Che10_iter = np.full((rv.raws, rv.num_samples), np.nan)
Che11_iter = np.full((rv.raws, rv.num_samples), np.nan)
Che12_iter = np.full((rv.raws, rv.num_samples), np.nan)
Che13_iter = np.full((rv.raws, rv.num_samples), np.nan)
Che14_iter = np.full((rv.raws, rv.num_samples), np.nan)
Che15_iter = np.full((rv.raws, rv.num_samples), np.nan)
Che16_iter = np.full((rv.raws, rv.num_samples), np.nan)
Che17_iter = np.full((rv.raws, rv.num_samples), np.nan)


for i in range(rv.num_samples):
    
                               # Code to generate the imax,tmax and dmax values #
                    
     velo = [rv.Tes_V100_samples[i],rv.Tes_V1_samples[i]]
     #print(Tes_velo)
     
     bounds =[]
     for j in range(len(rv.Tes_nodes)-1):
         #print(j)
         a = rv.Tes_nodes[j]/(velo[j]*rv.dt*(1/1000))
         b= rv.Tes_nodes[j+1]/(velo[j]*rv.dt*(1/1000))
         bounds.append(a)
         bounds.append(b)
     #print(Tes_bounds) 
     #print(len(Tes_bounds))
     
                                        # time step calculations #
     ts_val = [] # np.arange() values - all the time steps
     i_val =[] # no of timesetps in each segments
     x=0
     while x < len(bounds):
         a = np.arange(bounds[x],bounds[x+1])
         ts_val.append(a)
         b =len(a) 
         i_val.append(b)
         x = x+2
     #print((ts_val))
     #print(len(i_val))
     #print(i_val)
     #print("##i_values##",i_val)
     
     ts_no=[] # all timestep values in one list
     for xx in range(len(ts_val)):
         ts_no.extend(ts_val[xx]) 
     #print(ts_no) 
     #print(len(ts_no))
     
     #i_iter: i iteration timesteps per cycles for no of samples
     i_iter[0:len(ts_no),i]=ts_no
    
     #The most important part: calculating total no of timesteps
     i_sum=[]
    #values for the while loop
     for y in range(len(i_val)):
         #print(i)
         sum_ = np.sum(i_val[0:y+1])
         i_sum.append(sum_)
     #print("Values need to run the river code",i_sum) #"for the last i value,additional +1 has added")
      #print("i_values",len(i_values))
    
     #Teslin_imax values = the last item of the Tes_i_sum
     imax.append(i_sum[len(i_sum)-1])
     #print(imax)
     
                                                #   time_calculation    #
     #print("Tes_imax_values:",Tes_imax[i])
     t_ts = np.arange(0,imax[i])*rv.dt*(1/86400) # time is in days
     time_ = t_ts.tolist()
     
     #tmax 
     tmax.append(time_[len(time_)-1])
     #print(tmax)
     
     #t_iter: time for iteration steps per cycles for no of samples - 2D array
     t_iter[0:len(time_),i]=time_
     #print(t_iter[0:len(time_),i])
         
                                                #   distance_calculation    #
     dist_ts =[]
     for d in range(len(ts_val)):
         v_ts = velo[d]*rv.dt*(1/1000) # distance for dt time value: distance is in km, 1000 to convert m to km
         d_ts = ts_val[d]*v_ts
         #print(d_ts)
         dis_= d_ts.tolist()
         dist_ts.append(dis_)
     #print(dist_ts)
     #print(type(dist_ts[0]))
        
     dist_=[]
     for dd in range(len(dist_ts)):
         dist_.extend(dist_ts[dd]) 
     #print(dist_) 
     #print(len(dist_))
     dmax.append(dist_[len(dist_)-1])
     
     #d_iter: distance for iteration steps per cycles for no of samples - 2D array
     d_iter[0:len(dist_),i]=dist_
     #print(d_iter[0:len(dist_),i])
     
     
        
                         # Creating arrays for C_Conc,Tau and Prod  #
     
     C_Conc = np.ones((imax[i],rv.cmax)) 
     Tau = np.ones((imax[i],rv.cmax))
     Prod = np.ones((imax[i],rv.cmax))
     
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
               (rv.S_Phe_f_CDOM_samples[i] * C_Conc[:,14]) + \
               (rv.C_Phe_f_CDOM_samples[i] * C_Conc[:,15])
     #Phenol
     C_Conc[:,16] = C_Conc[:,13]+C_Conc[:,14]+C_Conc[:,15]      
     #DOC - Does not have Inorganic Carbon-8
     C_Conc[:,17] = (C_Conc[:,0] + C_Conc[:,1] + C_Conc[:,2] + C_Conc[:,3] +
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
     Prod[:,0] = rv.prot_Prod_samples[i]
     #Polypeptide
     #Prod[0,1] = -
     #Amino acid
     #Prod[:,2] = -
     #Heteroploycondensate(HPC)
     #Prod[:,3] = -
     #Polysaccharides
     Prod[:,4] = rv.poly_Prod_samples[i]
     #Oligosaccharides
     #Prod[:,5] = -
     #Humic acid
     #Prod[:,6] = -
     #Monosaccharides
     #Prod[:,7] = -
     #Inorganic carbon
     #Prod[:,8] = -
     #Lipids
     Prod[:,9] = rv.lip_Prod_samples[i]
     #Pigments
     Prod[:,10] = rv.pig_Prod_samples[i]
     #Porphyrin
     #Prod[:,11] = -
     #CDOM
     #Prod[:,12] = -
     #V-Phenol
     Prod[:,13] = rv.V_Phe_Prod_samples[i]
     #S-Phenol
     Prod[:,14] = rv.S_Phe_Prod_samples[i]
     #C-Phenol
     Prod[:,15] = rv.C_Phe_Prod_samples[i]
    
                                #  Entering to the river code  #
     j=0
     while j < imax[i]:
        #Protein 
        Conc0 =  f1(rv.dt,C_Conc[j,0],Tau[j,0],Prod[j,0])
        #Chem1:Peptides, two chanells
        Prod[j,1] = C_Conc[j,0]/Tau[j,0]
        Conc1 =  f1(rv.dt,C_Conc[j,1],Tau[j,1],Prod[j,1])
        #Chem2:A.Acids
        Prod[j,2] = C_Conc[j,1]/(Tau[j,1]*2)
        Conc2 =  f1(rv.dt,C_Conc[j,2],Tau[j,2],Prod[j,2])
        #Chem3:Heteropolycondensate
        Prod[j,3] = C_Conc[j,1]/(Tau[j,1]*2)+C_Conc[j,5]/(Tau[j,5]*(13/3))+C_Conc[j,9]/(Tau[j,9]*(13/3))
        Conc3 =  f1(rv.dt,C_Conc[j,3],Tau[j,3],Prod[j,3])
        #Chem4:Polysaccharides
        Conc4 =  f1(rv.dt,C_Conc[j,4],Tau[j,4],Prod[j,4])
        #Chem5:Oligosaccharide*** subject to decay by two processes, so tau is half of the value.
        Prod[j,5] = C_Conc[j,4]/Tau[j,4]
        Conc5 =  f1(rv.dt,C_Conc[j,5],Tau[j,5],Prod[j,5])
        #Chem6:Humic acid
        Prod[j,6] = C_Conc[j,3]/Tau[j,3]
        Conc6 =  f1(rv.dt,C_Conc[j,6],Tau[j,6],Prod[j,6])
        #Chem7:Monosachcharides
        Prod[j,7] = C_Conc[j,5]/(Tau[j,5]*(13/10))
        Conc7 =  f1(rv.dt,C_Conc[j,7],Tau[j,7],Prod[j,7])
        #Chem8:Inorganic Carbon
        Prod[j,8] = (C_Conc[j,7]/Tau[j,7])+(C_Conc[j,9]/(Tau[j,9])*(13/10))+(C_Conc[j,11]/Tau[j,11])
        Conc8 =  f1(rv.dt,C_Conc[j,8],Tau[j,8],Prod[j,8])
        #Chem9:Lipids
        Conc9 =  f1(rv.dt,C_Conc[j,9],Tau[j,9],Prod[j,9])
        #Chem10:Pigments
        Conc10 =  f1(rv.dt,C_Conc[j,10],Tau[j,10],Prod[j,10])
        #Chem11:Porpherin
        Prod[j,11] = C_Conc[j,10]/Tau[j,10]
        Conc11 =  f1(rv.dt,C_Conc[j,11],Tau[j,11],Prod[j,11])
        #Vphenol
        Conc13 =  f1(rv.dt,C_Conc[j,13],Tau[j,13],Prod[j,13])
        #Sphenol
        Conc14 =  f1(rv.dt,C_Conc[j,14],Tau[j,14],Prod[j,14])
        #Cphenol
        Conc15 =  f1(rv.dt,C_Conc[j,15],Tau[j,15],Prod[j,15])
        #CDOM 
        Conc12 =(rv.prot_f_CDOM_samples[i] * Conc0 ) + \
                (rv.hpc_f_CDOM_samples[i] * Conc3) + \
                (rv.pig_f_CDOM_samples[i] * Conc10) + \
                (rv.por_f_CDOM_samples[i] * Conc11) + \
                (rv.V_Phe_f_CDOM_samples[i] * Conc13) + \
                (rv.S_Phe_f_CDOM_samples[i] * Conc14) + \
                (rv.C_Phe_f_CDOM_samples[i] * Conc15)
        #Phenol
        Conc16 =  Conc13 + Conc14 + Conc15
        #DOC   - Does not have Inorganic Carbon-8
        Conc17 = (Conc0 + Conc1 + Conc2 + Conc3 +
                  Conc4 + Conc5 + Conc6 + Conc7 +
                  Conc9 + Conc10 + Conc11 + Conc13 +
                  Conc14 + Conc15)
        #
        Conc18 =  f1(rv.dt,C_Conc[i,18],Tau[i,18],Prod[i,18])
            
        if j < imax[i] - 1:
            C_Conc[j+1,0]=Conc0
            C_Conc[j+1,1]=Conc1
            C_Conc[j+1,2]=Conc2 
            C_Conc[j+1,3]=Conc3
            C_Conc[j+1,4]=Conc4
            C_Conc[j+1,5]=Conc5
            C_Conc[j+1,6]=Conc6
            C_Conc[j+1,7]=Conc7
            C_Conc[j+1,8]=Conc8
            C_Conc[j+1,9]=Conc9
            C_Conc[j+1,10]=Conc10
            C_Conc[j+1,11]=Conc11
            C_Conc[j+1,12]=Conc12
            C_Conc[j+1,13]=Conc13
            C_Conc[j+1,14]=Conc14
            C_Conc[j+1,15]=Conc15
            C_Conc[j+1,16]=Conc16
            C_Conc[j+1,17]=Conc17

        j=j+1
     
     #Inserting end values at river mouth lists for all chemical
     Conc_Che0.append(C_Conc[j-1,0])
     Conc_Che1.append(C_Conc[j-1,1])
     Conc_Che2.append(C_Conc[j-1,2])
     Conc_Che3.append(C_Conc[j-1,3])
     Conc_Che4.append(C_Conc[j-1,4])
     Conc_Che5.append(C_Conc[j-1,5])
     Conc_Che6.append(C_Conc[j-1,6])
     Conc_Che7.append(C_Conc[j-1,7])
     Conc_Che8.append(C_Conc[j-1,8])
     Conc_Che9.append(C_Conc[j-1,9])
     Conc_Che10.append(C_Conc[j-1,10])
     Conc_Che11.append(C_Conc[j-1,11])
     Conc_Che12.append(C_Conc[j-1,12])
     Conc_Che13.append(C_Conc[j-1,13])
     Conc_Che14.append(C_Conc[j-1,14])
     Conc_Che15.append(C_Conc[j-1,15])
     Conc_Che16.append(C_Conc[j-1,16])
     Conc_Che17.append(C_Conc[j-1,17])
     
     #For all chemicals: mapping all iteration values for n no of samples into 2D array
     Che0_iter[0:len(C_Conc[0:j,0]),i]=C_Conc[0:j,0]    
     Che1_iter[0:len(C_Conc[0:j,1]),i]=C_Conc[0:j,1]
     Che2_iter[0:len(C_Conc[0:j,2]),i]=C_Conc[0:j,2]
     Che3_iter[0:len(C_Conc[0:j,3]),i]=C_Conc[0:j,3]
     Che4_iter[0:len(C_Conc[0:j,4]),i]=C_Conc[0:j,4]
     Che5_iter[0:len(C_Conc[0:j,5]),i]=C_Conc[0:j,5]
     Che6_iter[0:len(C_Conc[0:j,6]),i]=C_Conc[0:j,6]
     Che7_iter[0:len(C_Conc[0:j,7]),i]=C_Conc[0:j,7]
     Che8_iter[0:len(C_Conc[0:j,8]),i]=C_Conc[0:j,8]
     Che9_iter[0:len(C_Conc[0:j,9]),i]=C_Conc[0:j,9]
     Che10_iter[0:len(C_Conc[0:j,10]),i]=C_Conc[0:j,10]
     Che11_iter[0:len(C_Conc[0:j,11]),i]=C_Conc[0:j,11]
     Che12_iter[0:len(C_Conc[0:j,12]),i]=C_Conc[0:j,12]
     Che13_iter[0:len(C_Conc[0:j,13]),i]=C_Conc[0:j,13]
     Che14_iter[0:len(C_Conc[0:j,14]),i]=C_Conc[0:j,14]
     Che15_iter[0:len(C_Conc[0:j,15]),i]=C_Conc[0:j,15]
     Che16_iter[0:len(C_Conc[0:j,16]),i]=C_Conc[0:j,16]
     Che17_iter[0:len(C_Conc[0:j,17]),i]=C_Conc[0:j,17]
     
#Print check on the end values
#print(imax)
#print(tmax)
#print(dmax)
#print(Conc_Che0)
 
#Print check on the iteration arrays values
#print(t_iter)
#print(len(t_iter))
#print(d_iter)
#print(len(d_iter))

plt.plot(dist_,C_Conc[0:j,0],'k--')
plt.show()
plt.plot(time_,Conc_Che0,'r--')
plt.show()


#                                     ***   Creating CSV files     ***
#                                 -------------------------------------
#Creating a list for no of samples
Sample_no = np.arange(0,rv.num_samples)

#Creating a CSV file for Creating lists for the end values of timesteps, time, distance and chemical concentrations
River_mouth_values = {
    'no_of_the sample':Sample_no,
    'total_timesteps':imax,
    'total_time / km':tmax,
    'total_distance/km':dmax,
    'Protein / micromolar':Conc_Che0,
    'Polypeptide / micromolar' :Conc_Che1,
    'Amino acid / micromolar':Conc_Che2,
    'Heteropolycondensate / micromolar':Conc_Che3,
    'Polysaccharide / micromolar':Conc_Che4,
    'Oligosaccharide / micromolar':Conc_Che5,
    'Humic acid / micromolar':Conc_Che6,
    'Monosaccharide / micromolar':Conc_Che7,
    'Inorganic carbon / micromolar':Conc_Che8,
    'Lipid / micromolar':Conc_Che9,
    'Pigment / micromolar':Conc_Che10,
    'Porpherin / micromolar':Conc_Che11,
    'CDOM / micromolar':Conc_Che12,
    'VPhenol / micromolar':Conc_Che13,
    'SPhenol / micromolar':Conc_Che14,
    'CPhenol / micromolar':Conc_Che15,
    'Phenol / micromolar':Conc_Che16,
    'DOC / micromolar':Conc_Che17,
          }
df = pd.DataFrame(River_mouth_values)
#print(df)
df.to_csv('output/teslin/Teslin_River_Mouth_Values.csv')

#***************
#Creating CSV file/es for all iterations for "n" no of sample runs
#i_iter
new_array = np.insert(i_iter,0,Sample_no,axis=0)
file_path = 'output/teslin/iterations/time_steps.csv'
np.savetxt(file_path, new_array, delimiter=",")

#t_iter
new_array = np.insert(t_iter,0,Sample_no,axis=0)
file_path = 'output/teslin/iterations/time.csv'
np.savetxt(file_path, new_array, delimiter=",")

#d_iter
new_array = np.insert(d_iter,0,Sample_no,axis=0)
file_path = 'output/teslin/iterations/distance.csv'
np.savetxt(file_path, new_array, delimiter=",")

#Che0_iter
new_array = np.insert(Che0_iter,0,Sample_no,axis=0)
file_path = 'output/teslin/iterations/0_protein.csv'
np.savetxt(file_path, new_array, delimiter=",")

#Che1_iter
new_array = np.insert(Che1_iter,0,Sample_no,axis=0)
file_path = 'output/teslin/iterations/1_polypeptide.csv'
np.savetxt(file_path, new_array, delimiter=",")

#Che2_iter
new_array = np.insert(Che2_iter,0,Sample_no,axis=0)
file_path = 'output/teslin/iterations/2_amino_acid.csv'
np.savetxt(file_path, new_array, delimiter=",")

#Che3_iter
new_array = np.insert(Che3_iter,0,Sample_no,axis=0)
file_path = 'output/teslin/iterations/3_heteropolycondensate.csv'
np.savetxt(file_path, new_array, delimiter=",")

#Che4_iter
new_array = np.insert(Che4_iter,0,Sample_no,axis=0)
file_path = 'output/teslin/iterations/4_polysaccharide.csv'
np.savetxt(file_path, new_array, delimiter=",")

#Che5_iter
new_array = np.insert(Che5_iter,0,Sample_no,axis=0)
file_path = 'output/teslin/iterations/5_oligosaccharide.csv'
np.savetxt(file_path, new_array, delimiter=",")

#Che6_iter
new_array = np.insert(Che6_iter,0,Sample_no,axis=0)
file_path = 'output/teslin/iterations/6_humic_acid.csv'
np.savetxt(file_path, new_array, delimiter=",")

#Che7_iter
new_array = np.insert(Che7_iter,0,Sample_no,axis=0)
file_path = 'output/teslin/iterations/7_monosaccharide.csv'
np.savetxt(file_path, new_array, delimiter=",")

#Che8_iter
new_array = np.insert(Che8_iter,0,Sample_no,axis=0)
file_path = 'output/teslin/iterations/8_inorganic_carbon.csv'
np.savetxt(file_path, new_array, delimiter=",")

#Che9_iter
new_array = np.insert(Che9_iter,0,Sample_no,axis=0)
file_path = 'output/teslin/iterations/9_lipids.csv'
np.savetxt(file_path, new_array, delimiter=",")

#Che10_iter
new_array = np.insert(Che10_iter,0,Sample_no,axis=0)
file_path = 'output/teslin/iterations/10_pigments.csv'
np.savetxt(file_path, new_array, delimiter=",")

#Che11_iter
new_array = np.insert(Che11_iter,0,Sample_no,axis=0)
file_path = 'output/teslin/iterations/11_porphyrin.csv'
np.savetxt(file_path, new_array, delimiter=",")

#Che12_iter
new_array = np.insert(Che12_iter,0,Sample_no,axis=0)
file_path = 'output/teslin/iterations/12_CDOM.csv'
np.savetxt(file_path, new_array, delimiter=",")

#Che13_iter
new_array = np.insert(Che13_iter,0,Sample_no,axis=0)
file_path = 'output/teslin/iterations/13_v-phenol.csv'
np.savetxt(file_path, new_array, delimiter=",")

#Che14_iter
new_array = np.insert(Che14_iter,0,Sample_no,axis=0)
file_path = 'output/teslin/iterations/14_s-phenol.csv'
np.savetxt(file_path, new_array, delimiter=",")

#Che15_iter
new_array = np.insert(Che15_iter,0,Sample_no,axis=0)
file_path = 'output/teslin/iterations/15_c-phenol.csv'
np.savetxt(file_path, new_array, delimiter=",")

#Che16_iter
new_array = np.insert(Che16_iter,0,Sample_no,axis=0)
file_path = 'output/teslin/iterations/16_phenol.csv'
np.savetxt(file_path, new_array, delimiter=",")

#Che17_iter
new_array = np.insert(Che17_iter,0,Sample_no,axis=0)
file_path = 'output/teslin/iterations/17_DOC.csv'
np.savetxt(file_path, new_array, delimiter=",")
