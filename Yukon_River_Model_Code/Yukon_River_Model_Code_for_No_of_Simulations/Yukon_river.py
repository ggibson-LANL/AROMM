import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import river_initial_monac_code as rv

#Importing the tributary river codes
import Teslin_river as tr
import Pelly_river as nd1rv               #node 1
import White_Donjec_river as nd2rv        #node 2
import Stewart_river as nd3rv             #node 3
import Porcupine_river as nd4rv           #node 4
import Tanana_river as nd5rv              #node 5
import Koyukuk_river as nd6rv             #node 6


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


#                        ***  lista and arrays to store the monte-carlo analysis results  ***
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

#Creating 2D arrays for store end value of all the chemicals for n no of samples
Chemicals_ = np.ones((rv.num_samples,rv.cmax))
#print(Chemicals_)

#                               ***  loop for the monte-carlo analysis   ***
#                                   -------------------------------------
for i in range(rv.num_samples):
    
                               # Code to generate the imax,tmax and dmax values #
                    
     velo = [rv.Yuk_V100_samples[i],rv.Yuk_V1_samples[i],rv.Yuk_V2_samples[i],rv.Yuk_V3_samples[i]
             ,rv.Yuk_V4_samples[i],rv.Yuk_V5_samples[i],rv.Yuk_V6_samples[i],rv.Yuk_V7_samples[i]
             ,rv.Yuk_V8_samples[i]]
     #print(velo)
     
     bounds =[]
     for j in range(len(rv.Yuk_nodes)-1):
         #print(j)
         a = rv.Yuk_nodes[j]/(velo[j]*rv.dt*(1/1000))
         b= rv.Yuk_nodes[j+1]/(velo[j]*rv.dt*(1/1000))
         bounds.append(a)
         bounds.append(b)
     #print(Yuk_bounds) 
     #print(len(Yuk_bounds))
     
                  # ----------------time step calculations------------------- #
                                        
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
     
    
     #The most important part: calculating total no of timesteps
     i_sum=[]
    #values for the while loop
     for y in range(len(i_val)):
         #print(i)
         sum_ = np.sum(i_val[0:y+1])
         i_sum.append(sum_)
     #print("Values need to run the river code",i_sum) #"for the last i value,additional +1 has added")
     #print("i_values",len(i_values))
    
     #Yukon_imax values = the last item of the Tes_i_sum
     imax.append(i_sum[len(i_sum)-1])
     #print(imax)

                        # ----------------time_calculation------------------- #     
                                        
     #print("Yuk_imax_values:",Yuk_imax[i])
     t_ts = np.arange(0,imax[i])*rv.dt*(1/86400) # time is in days
     time_ = t_ts.tolist()
     
     #tmax 
     tmax.append(time_[len(time_)-1])
     #print(tmax)
     
     
                    # ----------------distance_calculation------------------- #
         
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
     
          
             # ----------creating arrays for C_Conc,Tau and Prod------------- #        
                        
     C_Conc = np.ones((imax[i],rv.cmax)) 
     Tau = np.ones((imax[i],rv.cmax))
     Prod = np.ones((imax[i],rv.cmax))
     
        # --getting initial values for C_Conc, Tau and Prod for the chemica-- #     

     #Initial concentration of the species ia an array - units micrmolar C (micromol C / L)
     #DOC=rv.Yuk_DOC_samples[i]
     #print(tr.Conc_Che0[i])
     #Protein
     C_Conc[:,0]=tr.Conc_Che0[i]
     #Poltpeptide
     C_Conc[:,1]=tr.Conc_Che1[i]
     #Amino acid
     C_Conc[:,2]=tr.Conc_Che2[i]
     #Heteroploycondensate
     C_Conc[:,3]=tr.Conc_Che3[i]
     #Polysaccharides
     C_Conc[:,4]=tr.Conc_Che4[i]
     #Oligosaccharides
     C_Conc[:,5]=tr.Conc_Che5[i]
     #Humic acid
     C_Conc[:,6]=tr.Conc_Che6[i]
     #Monosaccharides
     C_Conc[:,7]=tr.Conc_Che7[i]
     #Inorganic carbon
     C_Conc[:,8]=tr.Conc_Che8[i]
     #Lipids
     C_Conc[:,9]=tr.Conc_Che9[i]
     #Pigments
     C_Conc[:,10]=tr.Conc_Che10[i]
     #Porpherin
     C_Conc[:,11]=tr.Conc_Che11[i]
     #VPhenol
     C_Conc[:,13]=tr.Conc_Che13[i]
     #SPhenol
     C_Conc[:,14]=tr.Conc_Che14[i]
     #CPhenol
     C_Conc[:,15]=tr.Conc_Che15[i]
     #CDOM
     C_Conc[:,12]=tr.Conc_Che12[i]
     #Phenol
     C_Conc[:,16]=tr.Conc_Che16[i]    
     #DOC - Does not have Inorganic Carbon-8
     C_Conc[:,17]=tr.Conc_Che17[i]
       
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
     #Production values are very closer to zero (both used and not used cases)
     #Protein
     Prod[:,0] = tr.Prod0[i]     # used in loop
     #Polypeptide
     Prod[:,1] = tr.Prod1[i]     # value not used in loop, calculated
     #Amino acid
     Prod[:,2] = tr.Prod2[i]     # value not used in loop, calculated
     #Heteroploycondensate(HPC)
     Prod[:,3] = tr.Prod3[i]     # value not used in loop, calculated
     #Polysaccharides
     Prod[:,4] = tr.Prod4[i]     # used in loop
     #Oligosaccharides
     Prod[:,5] = tr.Prod5[i]     # value not used in loop, calculated
     #Humic acid
     Prod[:,6] = tr.Prod6[i]     # value not used in loop, calculated
     #Monosaccharides
     Prod[:,7] = tr.Prod7[i]     # value not used in loop, calculated
     #Inorganic carbon
     Prod[:,8] = tr.Prod8[i]     # value not used in loop, calculated
     #Lipids
     Prod[:,9] = tr.Prod9[i]     # used in loop
     #Pigments
     Prod[:,10] = tr.Prod10[i]   # value not used in loop, calculated
     #Porphyrin
     Prod[:,11] = tr.Prod11[i]   # value not used in loop, calculated
     #CDOM
     #Prod[:,12] = -
     #V-Phenol
     Prod[:,13] = tr.Prod13[i]   # used in loop
     #S-Phenol
     Prod[:,14] = tr.Prod14[i]   # used in loop
     #C-Phenol
     Prod[:,15] = tr.Prod15[i]   # used in loop
     
                    # ----------------distance_calculation------------------- #     
     
    
                                #  Entering to the river code  #
     j=0
     
     while j < i_sum[1]:
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
        #Conc18 =  f1(rv.dt,C_Conc[i,18],Tau[i,18],Prod[i,18])
            
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
        
     # ---------------------------------Node_01------------------------------------------------  
     #node_1=Pelly_river
     #print("  ")
     #print("node_1:Pelly_river")
     #print("j_value",j)
     #print("node_1,stem:trib mixing")
     
     #Mixing
     C_Conc[j,:] = (1-rv.Pel_Df_samples[i])*C_Conc[j,:] + rv.Pel_Df_samples[i]*nd1rv.Chemicals_[i,:]     

     # ----------------------------------------------------------------------------------------      

     while j < i_sum[2]:
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
        #Conc18 =  f1(rv.dt,C_Conc[i,18],Tau[i,18],Prod[i,18])
            
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
        
     # ---------------------------------Node_02------------------------------------------------ 
     #node_2=White+Donjec_river
     #print("  ")
     #print("node_2:White+Donjec_river")
     #print("j_value",j)
     #print("node_2,stem:trib mixing")
 
     #Mixing
     C_Conc[j,:] = (1-rv.W_D_Df_samples[i])*C_Conc[j,:] + rv.W_D_Df_samples[i]*nd2rv.Chemicals_[i,:] 

     # ----------------------------------------------------------------------------------------  
         
     while j < i_sum[3]:
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
        #Conc18 =  f1(rv.dt,C_Conc[i,18],Tau[i,18],Prod[i,18])
            
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
        
     # ---------------------------------Node_03------------------------------------------------ 
     #node_3=Stewart_river
     #print("  ")
     #print("node_3:Stewart_river")
     #print("j_value",j)
     #print("node_3,stem:trib mixing")
     
     #Mixing
     C_Conc[j,:] = (1-rv.Ste_Df_samples[i])*C_Conc[j,:] + rv.Ste_Df_samples[i]*nd3rv.Chemicals_[i,:]  

     # ----------------------------------------------------------------------------------------      
         
     while j < i_sum[4]:
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
        #Conc18 =  f1(rv.dt,C_Conc[i,18],Tau[i,18],Prod[i,18])
            
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
     # ---------------------------------Node_04------------------------------------------------  
     #node_4=Porcupine_river
     #print("  ")
     #print("node_4:Porcupine_river")
     #print("j_value",j)
     #print("node_4,stem:trib mixing")    

     #Mixing
     C_Conc[j,:] = (1-rv.Por_Df_samples[i])*C_Conc[j,:] + rv.Por_Df_samples[i]*nd4rv.Chemicals_[i,:] 

     # ----------------------------------------------------------------------------------------     
         
     while j < i_sum[5]:
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
        #Conc18 =  f1(rv.dt,C_Conc[i,18],Tau[i,18],Prod[i,18])
            
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
        
     # ---------------------------------Node_05------------------------------------------------ 
     #node_5=Tanana_river
     #print("  ")
     #print("node_5:Tanana_river")
     #print("j_value",j)
     #print("node_5,stem:trib mixing")

     #Mixing
     C_Conc[j,:] = (1-rv.Tan_Df_samples[i])*C_Conc[j,:] + rv.Tan_Df_samples[i]*nd5rv.Chemicals_[i,:] 

     # ----------------------------------------------------------------------------------------          
         
     while j < i_sum[6]:  
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
        #Conc18 =  f1(rv.dt,C_Conc[i,18],Tau[i,18],Prod[i,18])
            
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
     
     # ---------------------------------Node_06------------------------------------------------ 
     #node_6=Koyukuk_river
     #print("  ")
     #print("node_6:Pelly_river")
     #print("j_value",j)
     #print("node_6,stem:trib mixing")
     
     #Mixing
     C_Conc[j,:] = (1-rv.Koy_Df_samples[i])*C_Conc[j,:] + rv.Koy_Df_samples[i]*nd6rv.Chemicals_[i,:] 

     # ----------------------------------------------------------------------------------------     
                        
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
        #Conc18 =  f1(rv.dt,C_Conc[i,18],Tau[i,18],Prod[i,18])
            
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
        
# ---------------------------------Filling the results into lists and arrays------------------------------------------------        
           
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
     
     #Creating 2D arrays for store end value of all the chemicals for n no of samples
     Chemicals_[i,0] = C_Conc[j-1,0]
     Chemicals_[i,1] = C_Conc[j-1,1]
     Chemicals_[i,2] = C_Conc[j-1,2]
     Chemicals_[i,3] = C_Conc[j-1,3]
     Chemicals_[i,4] = C_Conc[j-1,4]
     Chemicals_[i,5] = C_Conc[j-1,5]
     Chemicals_[i,6] = C_Conc[j-1,6]
     Chemicals_[i,7] = C_Conc[j-1,7]
     Chemicals_[i,8] = C_Conc[j-1,8]
     Chemicals_[i,9] = C_Conc[j-1,9]
     Chemicals_[i,10] = C_Conc[j-1,10]
     Chemicals_[i,11] = C_Conc[j-1,11]
     Chemicals_[i,12] = C_Conc[j-1,12]
     Chemicals_[i,13] = C_Conc[j-1,13]
     Chemicals_[i,14] = C_Conc[j-1,14]
     Chemicals_[i,15] = C_Conc[j-1,15]
     Chemicals_[i,16] = C_Conc[j-1,16]
     Chemicals_[i,17] = C_Conc[j-1,17]              
     #print(Chemicals_[i,17])
          
#Print check on the end values
#print(imax)
#print(tmax)
#print(dmax)
#print(Conc_Che17)
#print(type(Conc_Che17))
print(len(Conc_Che17))
 
#Print check on the iteration arrays values
#print(t_iter)
#print(len(t_iter))
#print(d_iter)
#print(len(d_iter))


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
df.to_csv('output/yukon/Yukon_River_Mouth_Values.csv')
df.to_csv('output/analysis/Yukon_River_Mouth_Values.csv')


#Creating input CSV file for Monte Carlo Analysis
InputData = {
    'Teslin_V1_ms_1':rv.Tes_V1_samples,
    'Pelly_V1_ms_1':rv.Pel_V1_samples,
    'White_Donjec_V1_ms_1':rv.W_D_V1_samples,
    'Stewart_V1_ms_1':rv.Ste_V1_samples,
    'Porcupine_V1_ms_1':rv.Por_V1_samples,
    'Tanana_V1_ms_1':rv.Tan_V1_samples,
    'Koyukuk_V1_ms_1':rv.Koy_V1_samples,
    'Yukon_V8_ms_1':rv.Yuk_V8_samples, 
    'Yukon_Flowtime_days':tmax,
    'Pelly_Df':rv.Pel_Df_samples,
    'White+Donjec_Df':rv.W_D_Df_samples,
    'Stewart_Df':rv.Ste_Df_samples,
    'Porcupine_Df':rv.Por_Df_samples,
    'Tanana_Df':rv.Tan_Df_samples,
    'Koyukuk_Df':rv.Koy_Df_samples,
    'Protein_turnovertime_seconds':rv.prot_Tau_samples,
    'Polypeptide_turnovertime_seconds':rv.polypep_Tau_samples,
    'Amino_acid_turnovertime_seconds':rv.amino_acid_Tau_samples,
    'Heteroploycondensate(HPC)_turnovertime_seconds':rv.hpc_Tau_samples,
    'Polysaccharides_turnovertime_seconds':rv.poly_Tau_samples,
    'Oligosaccharides_turnovertime_seconds':rv.olig_Tau_samples,
    'Humic_acid_turnovertime_seconds':rv.humic_acid_Tau_samples,
    'Monosaccharides_turnovertime_seconds':rv.mono_Tau_samples,
    'Inorganic_carbon_turnovertime_seconds':rv.ino_C_Tau_samples,
    'Lipids_turnovertime_seconds':rv.lip_Tau_samples,
    'Pigments_turnovertime_seconds':rv.pig_Tau_samples,
    'Porphyrin_turnovertime_seconds':rv.por_Tau_samples,
    'V-Phenol_turnovertime_seconds':rv.V_Phe_Tau_samples,    
    'S-Phenol_turnovertime_seconds':rv.S_Phe_Tau_samples,   
    'C-Phenol_turnovertime_seconds':rv.C_Phe_Tau_samples, 
    'Protein_fraction':rv.prot_f_samples,
    'Heteroploycondensate(HPC)_fraction':rv.hpc_f_samples,
    'Polysaccharides_fraction':rv.poly_f_samples,
    'Lipids_fraction':rv.lip_f_samples,
    'Pigments_fraction':rv.pig_f_samples,
    'V-Phenol_fraction':rv.V_Phe_f_samples,    
    'S-Phenol_fraction':rv.S_Phe_f_samples,   
    'C-Phenol_fraction':rv.C_Phe_f_samples,   
    'Teslin_DOC_micromolar':rv.Tes_DOC_samples,
    'Pelly_DOC_micromolar':rv.Pel_DOC_samples,
    'White_Donjec_DOC_micromolar':rv.W_D_DOC_samples,
    'Stewart_DOC_micromolar':rv.Ste_DOC_samples,
    'Porcupine_DOC_micromolar':rv.Por_DOC_samples,
    'Tanana_DOC_micromolar':rv.Tan_DOC_samples,
    'Koyukuk_DOC_micromolar':rv.Koy_DOC_samples       
          }

df = pd.DataFrame(InputData)
df.to_csv('output/analysis/InputData.csv',index=False)

#Creating output CSV file for Monte Carlo Analysis
OutputData = {
    'Protein_micromolar':Conc_Che0,
    'Polypeptide_micromolar' :Conc_Che1,
    'Amino acid_micromolar':Conc_Che2,
    'Heteropolycondensate_micromolar':Conc_Che3,
    'Polysaccharide_micromolar':Conc_Che4,
    'Oligosaccharide_micromolar':Conc_Che5,
    'Humic acid_micromolar':Conc_Che6,
    'Monosaccharide_micromolar':Conc_Che7,
    'Inorganic carbon_micromolar':Conc_Che8,
    'Lipid_micromolar':Conc_Che9,
    'Pigment_micromolar':Conc_Che10,
    'Porpherin_micromolar':Conc_Che11,
    'CDOM_micromolar':Conc_Che12,
    'VPhenol_micromolar':Conc_Che13,
    'SPhenol_micromolar':Conc_Che14,
    'CPhenol_micromolar':Conc_Che15,
    'Phenol_micromolar':Conc_Che16,
    'DOC_micromolar':Conc_Che17,
          }

df = pd.DataFrame(OutputData)
df.to_csv('output/analysis/OutputData.csv',index=False)
