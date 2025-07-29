import real2_Primary_river as rv
import real2_node1_river as nd1rv
import real2_node2_river as nd2rv
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
dt = 10

#print("Real calculations starts here")
while i < 190000:
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

print(i)
#the node1 1:1 mixing 
print(rv.C_Conc[i,0])
print(nd1rv.rv.C_Conc[i,0])
#rv.C_Conc[i,:] =((1/1)*rv.C_Conc[i,:])+((1/1)*nd1rv.rv.C_Conc[i,:])
rv.C_Conc[i,:] =((1/1)*rv.C_Conc[i,:])+((1/1)*nd2rv.rv.C_Conc[i,:])

while i <= 200000:
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
   rv.Prod[i,8] = (rv.C_Conc[i,7]/rv.Tau[i,7])+(rv.C_Conc[i,9]/   (rv.Tau[i,9])*(13/10))+(rv.C_Conc[i,11]/rv.Tau[i,11])
   Conc8 =  f1(dt,rv.C_Conc[i,8],rv.Tau[i,8],rv.Prod[i,8])
   #Chem9:Lipids
   Conc9 =  f1(dt,rv.C_Conc[i,9],rv.Tau[i,9],rv.Prod[i,9])
   #Chem10:Pigments
   Conc10 =  f1(dt,rv.C_Conc[i,10],rv.Tau[i,10],rv.Prod[i,10])
   #Chem11:Porpherin
   rv.Prod[i,11] = rv.C_Conc[i,10]/rv.Tau[i,10]
   Conc11 =  f1(dt,rv.C_Conc[i,11],rv.Tau[i,11],rv.Prod[i,11])
   #CDOM 
   Conc12 =  0.1*Conc0 + (0.1*Conc3)+(Conc10)
   #Vphenol-comes with tiaga
   Conc13 =  f1(dt,rv.C_Conc[i,13],rv.Tau[i,13],rv.Prod[i,13])
   #Sphenol-Comes with tundra
   Conc14 =  f1(dt,rv.C_Conc[i,14],rv.Tau[i,14],rv.Prod[i,14])
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

   i=i+1

#the node2 1:1 mixing   
rv.C_Conc[i,:] =((1/1)*rv.C_Conc[i,:])+((1/1)*nd2rv.rv.C_Conc[i,:])

while i < 300000:
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
   rv.Prod[i,8] = (rv.C_Conc[i,7]/rv.Tau[i,7])+(rv.C_Conc[i,9]/   (rv.Tau[i,9])*(13/10))+(rv.C_Conc[i,11]/rv.Tau[i,11])
   Conc8 =  f1(dt,rv.C_Conc[i,8],rv.Tau[i,8],rv.Prod[i,8])
   #Chem9:Lipids
   Conc9 =  f1(dt,rv.C_Conc[i,9],rv.Tau[i,9],rv.Prod[i,9])
   #Chem10:Pigments
   Conc10 =  f1(dt,rv.C_Conc[i,10],rv.Tau[i,10],rv.Prod[i,10])
   #Chem11:Porpherin
   rv.Prod[i,11] = rv.C_Conc[i,10]/rv.Tau[i,10]
   Conc11 =  f1(dt,rv.C_Conc[i,11],rv.Tau[i,11],rv.Prod[i,11])
   #CDOM 
   Conc12 =  0.1*Conc0 + (0.1*Conc3)+(Conc10)
   #Vphenol-comes with tiaga
   Conc13 =  f1(dt,rv.C_Conc[i,13],rv.Tau[i,13],rv.Prod[i,13])
   #Sphenol-Comes with tundra
   Conc14 =  f1(dt,rv.C_Conc[i,14],rv.Tau[i,14],rv.Prod[i,14])
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

   i=i+1


#print("                  ")
#list of no of steps
time = np.arange(0,rv.imax)*dt
time = time.tolist()
#print("* Time list for the chemical",Chem_ls[j])
#print(time)

#Velocity for first hundred km
velo = 0.1
#print("The velocity is" , velo , "m/s")
#list of no of steps, in distance(km)
dist1 = np.arange(0,100000)*dt*velo*(1/1000)
dist1 = dist1.tolist()
#print("* Distance")
#print(dist1)
#print(rv.C_Conc[99999,9])
#print(np.log10(rv.C_Conc[99999,9]))

#Velocity for first hundred km
velo = 1
#print("The velocity is" , velo , "m/s")
#list of no of steps, in distance(km)
dist2 = np.arange(10000,290000)*dt*velo*(1/1000)
dist2 = dist2.tolist()
#print("* Distance")
#print(dist2)
#print(rv.C_Conc[100000,9])
#print(np.log10(rv.C_Conc[100000,9]))

#Velocity for first hundred km
velo = 0.1
#print("The velocity is" , velo , "m/s")
#list of no of steps, in distance(km)
dist3 = np.arange(2900000,3000000)*dt*velo*(1/1000)
dist3 = dist3.tolist()
#print("* Distance")
#print(dist3)
#print(dist2)
dist = dist1+dist2+dist3
#print(dist)


#velo = 1
#print("The velocity is" , velo , "m/s")
#list of no of steps, in distance(km)
#dist = np.arange(0,rv.imax)*rv1.dt*velo*(1/1000)
#dist = dist.tolist()
#print("* Distance")
#print(dist)



#list of concentration of the species
#0
#print("* Concentration list of","Chem0","in micromolar")
#Conc_Che = rv.C_Conc[:,0]
Conc_Che = np.log10(rv.C_Conc[:,0])
Conc_Che = Conc_Che.tolist()
#print(Conc_Che)
#1
Conc_Che1 = np.log10(rv.C_Conc[:,1])
Conc_Che1 = Conc_Che1.tolist()
#2
Conc_Che2 = np.log10(rv.C_Conc[:,2])
Conc_Che2 = Conc_Che2.tolist()
#3
Conc_Che3 = np.log10(rv.C_Conc[:,3])
Conc_Che3 = Conc_Che3.tolist()
#4
Conc_Che4 = np.log10(rv.C_Conc[:,4])
Conc_Che4 = Conc_Che4.tolist()
#5
Conc_Che5 = np.log10(rv.C_Conc[:,5])
Conc_Che5 = Conc_Che5.tolist()
#6
Conc_Che6 = np.log10(rv.C_Conc[:,6])
Conc_Che6 = Conc_Che6.tolist()
#7
Conc_Che7 = np.log10(rv.C_Conc[:,7])
Conc_Che7 = Conc_Che7.tolist()
#8

Conc_Che8 = np.log10(rv.C_Conc[:,8])
Conc_Che8 = Conc_Che8.tolist()
#9
Conc_Che9 = np.log10(rv.C_Conc[:,9])
Conc_Che9 = Conc_Che9.tolist()

#10
Conc_Che10 = np.log10(rv.C_Conc[:,10])
Conc_Che10 = Conc_Che10.tolist()
#11
Conc_Che11 = np.log10(rv.C_Conc[:,11])
Conc_Che11 = Conc_Che11.tolist()
#12
Conc_Che12 = np.log10(rv.C_Conc[:,12])
Conc_Che12 = Conc_Che12.tolist()
#13
Conc_Che13 = np.log10(rv.C_Conc[:,13])
Conc_Che13 = Conc_Che13.tolist()
#14
Conc_Che14 = np.log10(rv.C_Conc[:,14])
Conc_Che14 = Conc_Che14.tolist()
#15
Conc_Che15 = np.log10(rv.C_Conc[:,15])
Conc_Che15 = Conc_Che15.tolist()
#16
Conc_Che16 = np.log10(rv.C_Conc[:,16])
Conc_Che16 = Conc_Che16.tolist()
#17
Conc_Che17 = np.log10(rv.C_Conc[:,17])
Conc_Che17 = Conc_Che17.tolist()
#18
Conc_Che18 = np.log10(rv.C_Conc[:,18])
Conc_Che18 = Conc_Che18.tolist()
#19
Conc_Che19 = np.log10(rv.C_Conc[:,19])
Conc_Che19 = Conc_Che19.tolist()
#20
Conc_Che20 = np.log10(rv.C_Conc[:,20])
Conc_Che20 = Conc_Che20.tolist()
#21
Conc_Che21 = np.log10(rv.C_Conc[:,21])
Conc_Che21 = Conc_Che21.tolist()
#22
Conc_Che22 = np.log10(rv.C_Conc[:,22])
Conc_Che22 = Conc_Che22.tolist()
#23
Conc_Che23 = np.log10(rv.C_Conc[:,23])
Conc_Che23 = Conc_Che23.tolist()
#24
Conc_Che24 = np.log10(rv.C_Conc[:,24])
Conc_Che24 = Conc_Che24.tolist()
#25
Conc_Che25 = np.log10(rv.C_Conc[:,25])
Conc_Che25 = Conc_Che25.tolist()
#26
Conc_Che26 = np.log10(rv.C_Conc[:,26])
Conc_Che26 = Conc_Che26.tolist()
#27
Conc_Che27 = np.log10(rv.C_Conc[:,27])
Conc_Che27 = Conc_Che27.tolist()
#28
Conc_Che28 = np.log10(rv.C_Conc[:,28])
Conc_Che28 = Conc_Che28.tolist()
#29
Conc_Che29 = np.log10(rv.C_Conc[:,29])
Conc_Che29 = Conc_Che29.tolist()

#Sum of chemical group
#Total DOC**HPC Duplicating ????
TDOC =rv.C_Conc[:,0]+rv.C_Conc[:,1]+rv.C_Conc[:,2]+rv.C_Conc[:,3]+rv.C_Conc[:,4]+rv.C_Conc[:,5]+rv.C_Conc[:,6]+rv.C_Conc[:,7]+rv.C_Conc[:,9]+rv.C_Conc[:,10]+rv.C_Conc[:,11]+rv.C_Conc[:,13]+rv.C_Conc[:,14] 
Conc_Che30 = np.log10(TDOC)
Conc_Che30 = Conc_Che30.tolist()
#Total protein+polypep+A.A
TProtein=rv.C_Conc[:,0]+rv.C_Conc[:,1]+rv.C_Conc[:,2]
Conc_Che31 = np.log10(TProtein)
Conc_Che31 = Conc_Che31.tolist()
#Total Carbohydrates:Starch+Oligo+Mono
TCarb=rv.C_Conc[:,4]+rv.C_Conc[:,5]+rv.C_Conc[:,7]
Conc_Che32 = np.log10(TCarb)
Conc_Che32 = Conc_Che32.tolist()
#Total Lipids
TLip=rv.C_Conc[:,9]
Conc_Che33 = np.log10(TLip)
Conc_Che33 = Conc_Che33.tolist()
#Total Colored materials:CDOM
TCOM=rv.C_Conc[:,12]
Conc_Che34 = np.log10(TCOM)
Conc_Che34 = Conc_Che34.tolist()

#plt.plot(dist,Conc_Che,'k-',label='Chem0')
#plt.plot(dist,Conc_Che1,'r--',label='Chem1')
#plt.plot(dist,Conc_Che2,'b-.',label='Chem2')
#plt.plot(dist,Conc_Che3,'g-',label='Chem3')
#plt.plot(dist,Conc_Che4,'c--',label='Chem4')
#plt.plot(dist,Conc_Che5,'y-.',label='Chem5')
#plt.plot(dist,Conc_Che6,'k-.',label='Chem6')
#plt.plot(dist,Conc_Che7,'r-.',label='Chem7')
#plt.plot(dist,Conc_Che8,'b--',label='Chem8')
#plt.plot(dist,Conc_Che9,'c-.',label='Chem9')
#plt.plot(dist,Conc_Che10,'g-.',label='Chem10')
#plt.plot(dist,Conc_Che11,'b--',label='Chem11')
#plt.plot(dist,Conc_Che12,'k:',label='Chem12')
#Protein,Carbohydrates,Heteroplo(Humic,CDOM),Pigments
#plt.axis([0,1000,-1,4])
#plt.grid(True)
#plt.plot(time,Conc_Che)
#plt.ylabel("Log of concentration of chemicals in micromolar")
#plt.xlabel("Distance/km")
#plt.title("The plot of Log of chemicals concentrations vs distance")
#plt.legend(loc=2)
#plt.show()


#fig,axes= plt.subplots(nrows=2, ncols=2)
fig = plt.figure()
fig.patch.set_facecolor('grey')
ax  = fig.add_subplot(2,2,1)
ax1 = fig.add_subplot(2,2,2)
ax2 = fig.add_subplot(2,2,3)
ax3 = fig.add_subplot(2,2,4)
#ax.set_facecolor('xkcd:salmon')
ax.axis([0,3000,-1,4])
ax.plot(dist,Conc_Che3,color='grey',label='(HPC)')
ax.plot(dist,Conc_Che,'r-',label='Protein')
ax.plot(dist,Conc_Che1,'k-',label='Polypeptide')
ax.plot(dist,Conc_Che2,'k--',label='Amino acids')
ax.set_title('Protein group')
ax.set_ylabel("Log[C]/micromolarL-1")
#ax.set_xlabel("Distance/km")
ax.grid(True)
ax.legend(loc=9, fontsize='x-small')

ax1.axis([0,3000,-1,4])
ax1.plot(dist,Conc_Che3,color='grey',label='(HPC)')
ax1.plot(dist,Conc_Che4,'r-',label='Carbohydrates')
ax1.plot(dist,Conc_Che5,'k-',label='Oligosaccharides')
ax1.plot(dist,Conc_Che2,'k--',label='Monosaccharides')
ax1.set_title('Carbohydrates group')
#ax1.set_ylabel("Log[C]/micromolarL-1")
#ax1.set_xlabel("Distance/km")
ax1.grid(True)
ax1.legend(fontsize='x-small')

ax2.axis([0,3000,-5,4])
ax2.plot(dist,Conc_Che3,color='grey',label='(HPC)')
ax2.plot(dist,Conc_Che9,'r-',label='Lipid')
ax2.plot(dist,Conc_Che13,'k-',label='VPhenols')
ax2.plot(dist,Conc_Che14,'k--',label='SPhenols')
ax2.set_title('Lipids and phenols')
ax2.set_ylabel("Log[C]/micromolarL-1")
ax2.set_xlabel("Distance/km")
ax2.grid(True)
ax2.legend(fontsize='x-small')

ax3.axis([0,3000,-5,4])
ax3.plot(dist,Conc_Che3,color='grey',label='(HPC)')
ax3.plot(dist,Conc_Che12,'r-',label='CDOM')
ax3.plot(dist,Conc_Che10,'k-',label='Pigments')
ax3.plot(dist,Conc_Che11,'k--',label='Porphyrins')
ax3.set_title('Colored materials group')
#ax3.set_ylabel("Log[C]/micromolarL-1")
ax3.set_xlabel("Distance/km")
ax3.grid(True)
ax3.legend(fontsize='x-small')

plt.tight_layout()
plt.show()

####Plots for summations

#fig,axes= plt.subplots(nrows=2, ncols=2)
fig = plt.figure()
fig.patch.set_facecolor('grey')
ax  = fig.add_subplot(2,2,1)
ax1 = fig.add_subplot(2,2,2)
ax2 = fig.add_subplot(2,2,3)
ax3 = fig.add_subplot(2,2,4)
#ax.set_facecolor('xkcd:salmon')
ax.axis([0,3000,-1,4])
ax.plot(dist,Conc_Che30,color='grey',label='(DOC)')
ax.plot(dist,Conc_Che31,'k-',label='Proteins')
ax.set_title('Surface Active Proteins')
ax.set_ylabel("Log[C]/micromolarL-1")
#ax.set_xlabel("Distance/km")
ax.grid(True)
ax.legend(fontsize='x-small')

ax1.axis([0,3000,-1,4])
ax1.plot(dist,Conc_Che30,color='grey',label='(DOC)')
#ax1.plot(dist,Conc_Che32,'k-',label='Carbohydrates')
ax1.set_title('Carbohydrates not active')
#ax1.set_ylabel("Log[C]/micromolarL-1")
#ax1.set_xlabel("Distance/km")
ax1.grid(True)
ax1.legend(fontsize='x-small')

ax2.axis([0,3000,-5,4])
ax2.plot(dist,Conc_Che30,color='grey',label='(DOC)')
ax2.plot(dist,Conc_Che33,'k-',label='Lipids')
ax2.set_title('Surface active Lipids')
ax2.set_ylabel("Log[C]/micromolarL-1")
ax2.set_xlabel("Distance/km")
ax2.grid(True)
ax2.legend(fontsize='x-small')

ax3.axis([0,3000,-5,4])
ax3.plot(dist,Conc_Che30,color='grey',label='(DOC)')
ax3.plot(dist,Conc_Che34,'k-',label='CDOM')
ax3.set_title('Colored materials')
#ax3.set_ylabel("Log[C]/micromolarL-1")
ax3.set_xlabel("Distance/km")
ax3.grid(True)
ax3.legend(fontsize='x-small')

plt.tight_layout()
plt.show()
