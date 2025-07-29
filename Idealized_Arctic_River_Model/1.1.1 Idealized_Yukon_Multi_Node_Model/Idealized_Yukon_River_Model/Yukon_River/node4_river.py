import node4 as rv
import math
import numpy as np
import matplotlib.pyplot as plt

# Equation
# C =C0*exp(-dt/Tau) + P*Tau*(1-exp(-dt/Tau))

def f1(dt,C0,Tau,P):
    term1 = C0 * math.exp((-1*dt)/Tau)
    term2 = P * Tau * (1-math.exp((-1*dt)/Tau))
    #print("Decay value", term1)
    #print("Production term", term2)
    Conc = term1 + term2
    # print("The concentration is " , Conc)
    return Conc

#f1(1,1,1,1)
#Coding for calculating the concentration
#print(rv.imax)
#print(rv.cmax)
#print("Tau value for first chemical",rv.Tau[0,0])
#print("Prod value for first chemical", rv.Prod[0,0])
#print("Initial concentration value for first chemical", rv.C_Conc[0,0])
i = 0
#The time step
dt = 10000

#print("Real calculations starts here")
while i < rv.imax:
   #print("               ")
   #print("The step no:",i)
   #Concentraion of chemicals
   Conc0 =  f1(dt,rv.C_Conc[i,0],rv.Tau[i,0],rv.Prod[i,0])
   #Chem1:Peptides, two chanells
   rv.Prod[i,1] = rv.C_Conc[i,0]/rv.Tau[i,0]
   Conc1 =  f1(dt,rv.C_Conc[i,1],rv.Tau[i,1],rv.Prod[i,1])
   #Chem2:A.Acids
   rv.Prod[i,2] = rv.C_Conc[i,1]/(rv.Tau[i,1]*2)
   Conc2 =  f1(dt,rv.C_Conc[i,2],rv.Tau[i,2],rv.Prod[i,2])
   #Chem3:Heteropolycondensate
   rv.Prod[i,3] = rv.C_Conc[i,1]/(rv.Tau[i,1]*2)+rv.C_Conc[i,5]/(rv.Tau[i,5]*(13/3))+rv.C_Conc[i,9]/(rv.Tau[i,9]*(13/3))
   Conc3 =  f1(dt,rv.C_Conc[i,3],rv.Tau[i,3],rv.Prod[i,3])
   #Chem4:Starch
   Conc4 =  f1(dt,rv.C_Conc[i,4],rv.Tau[i,4],rv.Prod[i,4])
   #Chem5:Oligosaccharide*** subject to decay by two processes, so tau is half of the value.
   rv.Prod[i,5] = rv.C_Conc[i,4]/rv.Tau[i,4]
   Conc5 =  f1(dt,rv.C_Conc[i,5],rv.Tau[i,5],rv.Prod[i,5])
   #Chem6:Humic acid
   rv.Prod[i,6] = rv.C_Conc[i,3]/rv.Tau[i,3]
   Conc6 =  f1(dt,rv.C_Conc[i,6],rv.Tau[i,6],rv.Prod[i,6])
   #Chem7:Monosachcharides
   rv.Prod[i,7] = rv.C_Conc[i,5]/(rv.Tau[i,5]*(13/10))
   Conc7 =  f1(dt,rv.C_Conc[i,7],rv.Tau[i,7],rv.Prod[i,7])
   #Chem8:Inorganic Carbon
   rv.Prod[i,8] = (rv.C_Conc[i,7]/rv.Tau[i,7])+(rv.C_Conc[i,9]/(rv.Tau[i,9])*(13/10))+(rv.C_Conc[i,11]/rv.Tau[i,11])
   Conc8 =  f1(dt,rv.C_Conc[i,8],rv.Tau[i,8],rv.Prod[i,8])
   #Chem9:Lipids
   Conc9 =  f1(dt,rv.C_Conc[i,9],rv.Tau[i,9],rv.Prod[i,9])
   #Chem10:Pigments
   Conc10 =  f1(dt,rv.C_Conc[i,10],rv.Tau[i,10],rv.Prod[i,10])
   #Chem11:Porpherin
   rv.Prod[i,11] = rv.C_Conc[i,10]/rv.Tau[i,10]
   Conc11 =  f1(dt,rv.C_Conc[i,11],rv.Tau[i,11],rv.Prod[i,11])
   #Vphenol
   Conc13 =  f1(dt,rv.C_Conc[i,13],rv.Tau[i,13],rv.Prod[i,13])
   #Sphenol
   Conc14 =  f1(dt,rv.C_Conc[i,14],rv.Tau[i,14],rv.Prod[i,14])
   #CDOM 
   Conc12 =  0.1*Conc0 + (0.1*Conc3)+(Conc10)+((1/3)*Conc13)+((1/3)*(Conc14))
   Conc15 =  f1(dt,rv.C_Conc[i,15],rv.Tau[i,15],rv.Prod[i,15])
   Conc16 =  f1(dt,rv.C_Conc[i,16],rv.Tau[i,16],rv.Prod[i,16])
   Conc17 =  f1(dt,rv.C_Conc[i,17],rv.Tau[i,17],rv.Prod[i,17])
   Conc18 =  f1(dt,rv.C_Conc[i,18],rv.Tau[i,18],rv.Prod[i,18])
   Conc19 =  f1(dt,rv.C_Conc[i,19],rv.Tau[i,19],rv.Prod[i,19])
   Conc20 =  f1(dt,rv.C_Conc[i,20],rv.Tau[i,20],rv.Prod[i,20])
   Conc21 =  f1(dt,rv.C_Conc[i,21],rv.Tau[i,21],rv.Prod[i,21])
   Conc22 =  f1(dt,rv.C_Conc[i,22],rv.Tau[i,22],rv.Prod[i,22])
   Conc23 =  f1(dt,rv.C_Conc[i,23],rv.Tau[i,23],rv.Prod[i,23])
   Conc24 =  f1(dt,rv.C_Conc[i,24],rv.Tau[i,24],rv.Prod[i,24])
   Conc25 =  f1(dt,rv.C_Conc[i,25],rv.Tau[i,25],rv.Prod[i,25])
   Conc26 =  f1(dt,rv.C_Conc[i,26],rv.Tau[i,26],rv.Prod[i,26])
   Conc27 =  f1(dt,rv.C_Conc[i,27],rv.Tau[i,27],rv.Prod[i,27])
   Conc28 =  f1(dt,rv.C_Conc[i,28],rv.Tau[i,28],rv.Prod[i,28])
   Conc29 =  f1(dt,rv.C_Conc[i,29],rv.Tau[i,29],rv.Prod[i,29])
   
   if i < rv.imax - 1:
        rv.C_Conc[i+1,0]=Conc0
        rv.C_Conc[i+1,1]=Conc1
        rv.C_Conc[i+1,2]=Conc2 
        rv.C_Conc[i+1,3]=Conc3
        rv.C_Conc[i+1,4]=Conc4
        rv.C_Conc[i+1,5]=Conc5
        rv.C_Conc[i+1,6]=Conc6
        rv.C_Conc[i+1,7]=Conc7
        rv.C_Conc[i+1,8]=Conc8
        rv.C_Conc[i+1,9]=Conc9
        rv.C_Conc[i+1,10]=Conc10
        rv.C_Conc[i+1,11]=Conc11
        rv.C_Conc[i+1,12]=Conc12
        rv.C_Conc[i+1,13]=Conc13
        rv.C_Conc[i+1,14]=Conc14
        rv.C_Conc[i+1,15]=Conc15
        rv.C_Conc[i+1,16]=Conc16
        rv.C_Conc[i+1,17]=Conc17
        rv.C_Conc[i+1,18]=Conc18
        rv.C_Conc[i+1,19]=Conc19
        rv.C_Conc[i+1,20]=Conc20
        rv.C_Conc[i+1,21]=Conc21
        rv.C_Conc[i+1,22]=Conc22
        rv.C_Conc[i+1,23]=Conc23
        rv.C_Conc[i+1,24]=Conc24
        rv.C_Conc[i+1,25]=Conc25
        rv.C_Conc[i+1,26]=Conc26
        rv.C_Conc[i+1,27]=Conc27
        rv.C_Conc[i+1,28]=Conc28
        rv.C_Conc[i+1,29]=Conc29
        
   #print("Array of current concentration values",rv.C_Conc)
   i=i+1

TDOC =rv.C_Conc[:,0]+rv.C_Conc[:,1]+rv.C_Conc[:,2]+rv.C_Conc[:,3]+rv.C_Conc[:,4]+rv.C_Conc[:,5]+rv.C_Conc[:,6]+rv.C_Conc[:,7]+rv.C_Conc[:,9]+rv.C_Conc[:,10]+rv.C_Conc[:,11]+rv.C_Conc[:,13]+rv.C_Conc[:,14] 
Conc_Che30 = np.log10(TDOC)
Conc_Che30 = Conc_Che30.tolist()
#print(Conc_Che30[0],"+++logvalue")  
#print("node_4_[TDOC]_before_connecting_Primary_river_",Conc_Che30[i-1])

print("node_4_ending_i_value",i)