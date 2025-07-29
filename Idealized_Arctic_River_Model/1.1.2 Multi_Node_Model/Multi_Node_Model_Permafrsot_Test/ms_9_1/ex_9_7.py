#HighSolute(x3)+Tundra(x3)
import math
import numpy as np
import matplotlib.pyplot as plt

#No of distnce ponits: start+nodes+end
imax=13
#print("imax",imax)

#cmax=no of chemical species
cmax= 30
#print("cmax",cmax)

#Nodes ratio
#Node1
#V0 primary river ratio
V0=0.5
#tributary river ratio
V1=0.5

#Node2
#V20 primary ratio
V20=0.5
#tributary ratio
V2=0.5

#Node3
V30=0.5
V3=0.5

#Node4
V40=0.5
V4=0.5

#Node5
V50=0.5
V5=0.5

#Node6
V60=0.5
V6=0.5

#Node7
V70=0.5
V7=0.5

#Node8
V80=0.5
V8=0.5

#Node9
V90=0.5
V9=0.5

#Node10
V100=0.5
V10=0.5

#Node11
V110=0.5
V11=0.5


#Maximum no of time steps
n=1
#print("imax",imax)

#cmax=no of chemical species
cmax= 30
#print("cmax",cmax)

#Initial concentration of the species ia an array
#Concentrations in the micromolar
#Mountain_node
n1C_Conc = np.ones((n,cmax))
#High_solute_node
n2C_Conc = np.ones((n,cmax))
#Tundra_node
n3C_Conc = np.ones((n,cmax))

#Mountain node
#Protein
n1C_Conc[:,0]=5
#Poltpep
n1C_Conc[:,1]=0
#Amino_acid
n1C_Conc[:,2]=0
#Heteroploycondensate
n1C_Conc[:,3]=35
#Starch_poly
n1C_Conc[:,4]=4
#Oligosaccharides
n1C_Conc[:,5]=0
#Humic acid
n1C_Conc[:,6]=50
#Monosaccharides
n1C_Conc[:,7]=0
#Inorganic carbon
n1C_Conc[:,8]=0
#Lipids
n1C_Conc[:,9]=3
#Pigments
n1C_Conc[:,10]=0.7
#Porpherin
n1C_Conc[:,11]=0
#CDOM
#C_Conc[:,12]=0
#VPhenol
n1C_Conc[:,13]=1
#SPhenol
n1C_Conc[:,14]=0.000001

#High_solute_node
#Protein
n2C_Conc[:,0]=50*3
#Poltpep
n2C_Conc[:,1]=0
#Amino_acid
n2C_Conc[:,2]=0
#Heteroploycondensate
n2C_Conc[:,3]=350*3
#Starch_poly
n2C_Conc[:,4]=40*3
#Oligosaccharides
n2C_Conc[:,5]=0
#Humic acid
n2C_Conc[:,6]=500*3
#Monosaccharides
n2C_Conc[:,7]=0
#Inorganic carbon
n2C_Conc[:,8]=0
#Lipids
n2C_Conc[:,9]=30*3
#Pigments
n2C_Conc[:,10]=7*3
#Porpherin
n2C_Conc[:,11]=0
#CDOM
#C_Conc[:,12]=0
#VPhenol
n2C_Conc[:,13]=10*3
#SPhenol
n2C_Conc[:,14]=0.000001*3


#Tundra_node
#Protein
n3C_Conc[:,0]=15*3
#Poltpep
n3C_Conc[:,1]=0
#Amino_acid
n3C_Conc[:,2]=0
#Heteroploycondensate
n3C_Conc[:,3]=105*3
#Starch_poly
n3C_Conc[:,4]=12*3
#Oligosaccharides
n3C_Conc[:,5]=0
#Humic acid
n3C_Conc[:,6]=150*3
#Monosaccharides
n3C_Conc[:,7]=0
#Inorganic carbon
n3C_Conc[:,8]=0
#Lipids
n3C_Conc[:,9]=9*3
#Pigments
n3C_Conc[:,10]=7*3
#Porpherin
n3C_Conc[:,11]=0
#CDOM
#C_Conc[:,12]=0
#VPhenol
n3C_Conc[:,13]=0
#SPhenol
n3C_Conc[:,14]=3*3


#Initial concentration of the species ia an array
#Concentrations in the micromolar
C_Conc = np.ones((imax,cmax))
#C_Conc[:,:]=0
#Protein
C_Conc[0,0]=5
#Poltpep
C_Conc[0,1]=0
#Amino_acid
C_Conc[0,2]=0
#Heteroploycondensate
C_Conc[0,3]=35
print(C_Conc[0,3])
#Starch_poly
C_Conc[0,4]=4
#Oligosaccharides
C_Conc[0,5]=0
#Humic acid
C_Conc[0,6]=50
#Monosaccharides
C_Conc[0,7]=0
#Inorganic carbon
C_Conc[0,8]=0
#Lipids
C_Conc[0,9]=3
#Pigments
C_Conc[0,10]=0.7
#Porpherin
C_Conc[0,11]=0
#CDOM
#C_Conc[0,12]=0
#VPhenol
C_Conc[0,13]=1
#SPhenol
C_Conc[0,14]=0.000001

#mixing_node_1
#Protein
C_Conc[1,0]=(V0*C_Conc[0,0])+(V1*n1C_Conc[:,0])
#Poltpep
C_Conc[1,1]=(V0*C_Conc[0,1])+(V1*n1C_Conc[:,1])
#Amino_acid
C_Conc[1,2]=(V0*C_Conc[0,2])+(V1*n1C_Conc[:,2])
#Heteroploycondensate
C_Conc[1,3]=(V0*C_Conc[0,3])+(V1*n1C_Conc[:,3])
#Starch_poly
C_Conc[1,4]=(V0*C_Conc[0,4])+(V1*n1C_Conc[:,4])
#Oligosaccharides
C_Conc[1,5]=(V0*C_Conc[0,5])+(V1*n1C_Conc[:,5])
#Humic acid
C_Conc[1,6]=(V0*C_Conc[0,6])+(V1*n1C_Conc[:,6])
#Monosaccharides
C_Conc[1,7]=(V0*C_Conc[0,7])+(V1*n1C_Conc[:,7])
#Inorganic carbon
C_Conc[1,8]=(V0*C_Conc[0,8])+(V1*n1C_Conc[:,8])
#Lipids
C_Conc[1,9]=(V0*C_Conc[0,9])+(V1*n1C_Conc[:,9])
#Pigments
C_Conc[1,10]=(V0*C_Conc[0,10])+(V1*n1C_Conc[:,10])
#Porpherin
C_Conc[1,11]=(V0*C_Conc[0,11])+(V1*n1C_Conc[:,11])
#CDOM
#C_Conc[:,12]=0
#VPhenol
C_Conc[1,13]=(V0*C_Conc[0,13])+(V1*n1C_Conc[:,13])
#SPhenol
C_Conc[1,14]=(V0*C_Conc[0,14])+(V1*n1C_Conc[:,14])


#mixing_node_2
#Protein
C_Conc[2,0]=(V20*C_Conc[1,0])+(V2*n1C_Conc[:,0])
#Poltpep
C_Conc[2,1]=(V20*C_Conc[1,1])+(V2*n1C_Conc[:,1])
#Amino_acid
C_Conc[2,2]=(V20*C_Conc[1,2])+(V2*n1C_Conc[:,2])
#Heteroploycondensate
C_Conc[2,3]=(V20*C_Conc[1,3])+(V2*n1C_Conc[:,3])
#Starch_poly
C_Conc[2,4]=(V20*C_Conc[1,4])+(V2*n1C_Conc[:,4])
#Oligosaccharides
C_Conc[2,5]=(V20*C_Conc[1,5])+(V2*n1C_Conc[:,5])
#Humic acid
C_Conc[2,6]=(V20*C_Conc[1,6])+(V2*n1C_Conc[:,6])
#Monosaccharides
C_Conc[2,7]=(V20*C_Conc[1,7])+(V2*n1C_Conc[:,7])
#Inorganic carbon
C_Conc[2,8]=(V20*C_Conc[1,8])+(V2*n1C_Conc[:,8])
#Lipids
C_Conc[2,9]=(V20*C_Conc[1,9])+(V2*n1C_Conc[:,9])
#Pigments
C_Conc[2,10]=(V20*C_Conc[1,10])+(V2*n1C_Conc[:,10])
#Porpherin
C_Conc[2,11]=(V20*C_Conc[1,11])+(V2*n1C_Conc[:,11])
#CDOM
#C_Conc[:,12]=0
#VPhenol
C_Conc[2,13]=(V20*C_Conc[1,13])+(V2*n1C_Conc[:,13])
#SPhenol
C_Conc[2,14]=(V20*C_Conc[1,14])+(V2*n1C_Conc[:,14])

#mixing_node_3
#Protein
C_Conc[3,0]=(V30*C_Conc[2,0])+(V3*n1C_Conc[:,0])
#Poltpep
C_Conc[3,1]=(V30*C_Conc[2,1])+(V3*n1C_Conc[:,1])
#Amino_acid
C_Conc[3,2]=(V30*C_Conc[2,2])+(V3*n1C_Conc[:,2])
#Heteroploycondensate
C_Conc[3,3]=(V30*C_Conc[2,3])+(V3*n1C_Conc[:,3])
#Starch_poly
C_Conc[3,4]=(V30*C_Conc[2,4])+(V3*n1C_Conc[:,4])
#Oligosaccharides
C_Conc[3,5]=(V30*C_Conc[2,5])+(V3*n1C_Conc[:,5])
#Humic acid
C_Conc[3,6]=(V30*C_Conc[2,6])+(V3*n1C_Conc[:,6])
#Monosaccharides
C_Conc[3,7]=(V30*C_Conc[2,7])+(V3*n1C_Conc[:,7])
#Inorganic carbon
C_Conc[3,8]=(V30*C_Conc[2,8])+(V3*n1C_Conc[:,8])
#Lipids
C_Conc[3,9]=(V30*C_Conc[2,9])+(V3*n1C_Conc[:,9])
#Pigments
C_Conc[3,10]=(V30*C_Conc[2,10])+(V3*n1C_Conc[:,10])
#Porpherin
C_Conc[3,11]=(V30*C_Conc[2,11])+(V3*n1C_Conc[:,11])
#CDOM
#C_Conc[:,12]=0
#VPhenol
C_Conc[3,13]=(V30*C_Conc[2,13])+(V3*n1C_Conc[:,13])
#SPhenol
C_Conc[3,14]=(V30*C_Conc[2,14])+(V3*n1C_Conc[:,14])


#mixing_node_4
#Protein
C_Conc[4,0]=(V40*C_Conc[3,0])+(V4*n2C_Conc[:,0])
#Poltpep
C_Conc[4,1]=(V40*C_Conc[3,1])+(V4*n2C_Conc[:,1])
#Amino_acid
C_Conc[4,2]=(V40*C_Conc[3,2])+(V4*n2C_Conc[:,2])
#Heteroploycondensate
C_Conc[4,3]=(V40*C_Conc[3,3])+(V4*n2C_Conc[:,3])
#Starch_poly
C_Conc[4,4]=(V40*C_Conc[3,4])+(V4*n2C_Conc[:,4])
#Oligosaccharides
C_Conc[4,5]=(V40*C_Conc[3,5])+(V4*n2C_Conc[:,5])
#Humic acid
C_Conc[4,6]=(V40*C_Conc[3,6])+(V4*n2C_Conc[:,6])
#Monosaccharides
C_Conc[4,7]=(V40*C_Conc[3,7])+(V4*n2C_Conc[:,7])
#Inorganic carbon
C_Conc[4,8]=(V40*C_Conc[3,8])+(V4*n2C_Conc[:,8])
#Lipids
C_Conc[4,9]=(V40*C_Conc[3,9])+(V4*n2C_Conc[:,9])
#Pigments
C_Conc[4,10]=(V40*C_Conc[3,10])+(V4*n2C_Conc[:,10])
#Porpherin
C_Conc[4,11]=(V40*C_Conc[3,11])+(V4*n2C_Conc[:,11])
#CDOM
#C_Conc[:,12]=0
#VPhenol
C_Conc[4,13]=(V40*C_Conc[3,13])+(V4*n2C_Conc[:,13])
#SPhenol
C_Conc[4,14]=(V40*C_Conc[3,14])+(V4*n2C_Conc[:,14])

#mixing_node_5
#Protein
C_Conc[5,0]=(V50*C_Conc[4,0])+(V5*n2C_Conc[:,0])
#Poltpep
C_Conc[5,1]=(V50*C_Conc[4,1])+(V5*n2C_Conc[:,1])
#Amino_acid
C_Conc[5,2]=(V50*C_Conc[4,2])+(V5*n2C_Conc[:,2])
#Heteroploycondensate
C_Conc[5,3]=(V50*C_Conc[4,3])+(V5*n2C_Conc[:,3])
print(C_Conc[5,3])
#Starch_poly
C_Conc[5,4]=(V50*C_Conc[4,4])+(V5*n2C_Conc[:,4])
#Oligosaccharides
C_Conc[5,5]=(V50*C_Conc[4,5])+(V5*n2C_Conc[:,5])
#Humic acid
C_Conc[5,6]=(V50*C_Conc[4,6])+(V5*n2C_Conc[:,6])
#Monosaccharides
C_Conc[5,7]=(V50*C_Conc[4,7])+(V5*n2C_Conc[:,7])
#Inorganic carbon
C_Conc[5,8]=(V50*C_Conc[4,8])+(V5*n2C_Conc[:,8])
#Lipids
C_Conc[5,9]=(V50*C_Conc[4,9])+(V5*n2C_Conc[:,9])
#Pigments
C_Conc[5,10]=(V50*C_Conc[4,10])+(V5*n2C_Conc[:,10])
#Porpherin
C_Conc[5,11]=(V50*C_Conc[4,11])+(V5*n2C_Conc[:,11])
#CDOM
#C_Conc[:,12]=0
#VPhenol
C_Conc[5,13]=(V50*C_Conc[4,13])+(V5*n2C_Conc[:,13])
#SPhenol
C_Conc[5,14]=(V50*C_Conc[4,14])+(V5*n2C_Conc[:,14])


#mixing_node_6
#Protein
C_Conc[6,0]=(V60*C_Conc[5,0])+(V6*n2C_Conc[:,0])
#Poltpep
C_Conc[6,1]=(V60*C_Conc[5,1])+(V6*n2C_Conc[:,1])
#Amino_acid
C_Conc[6,2]=(V60*C_Conc[5,2])+(V6*n2C_Conc[:,2])
#Heteroploycondensate
C_Conc[6,3]=(V60*C_Conc[5,3])+(V6*n2C_Conc[:,3])
print(C_Conc[6,3])
#Starch_poly
C_Conc[6,4]=(V60*C_Conc[5,4])+(V6*n2C_Conc[:,4])
#Oligosaccharides
C_Conc[6,5]=(V60*C_Conc[5,5])+(V6*n2C_Conc[:,5])
#Humic acid
C_Conc[6,6]=(V60*C_Conc[5,6])+(V6*n2C_Conc[:,6])
#Monosaccharides
C_Conc[6,7]=(V60*C_Conc[5,7])+(V6*n2C_Conc[:,7])
#Inorganic carbon
C_Conc[6,8]=(V60*C_Conc[5,8])+(V6*n2C_Conc[:,8])
#Lipids
C_Conc[6,9]=(V60*C_Conc[5,9])+(V6*n2C_Conc[:,9])
#Pigments
C_Conc[6,10]=(V60*C_Conc[5,10])+(V6*n2C_Conc[:,10])
#Porpherin
C_Conc[6,11]=(V60*C_Conc[5,11])+(V6*n2C_Conc[:,11])
#CDOM
#C_Conc[:,12]=0
#VPhenol
C_Conc[6,13]=(V60*C_Conc[5,13])+(V6*n2C_Conc[:,13])
#SPhenol
C_Conc[6,14]=(V60*C_Conc[5,14])+(V6*n2C_Conc[:,14])


#mixing_node_7
#Protein
C_Conc[7,0]=(V70*C_Conc[6,0])+(V7*n2C_Conc[:,0])
#Poltpep
C_Conc[7,1]=(V70*C_Conc[6,1])+(V7*n2C_Conc[:,1])
#Amino_acid
C_Conc[7,2]=(V70*C_Conc[6,2])+(V7*n2C_Conc[:,2])
#Heteroploycondensate
C_Conc[7,3]=(V70*C_Conc[6,3])+(V7*n2C_Conc[:,3])
#Starch_poly
C_Conc[7,4]=(V70*C_Conc[6,4])+(V7*n2C_Conc[:,4])
#Oligosaccharides
C_Conc[7,5]=(V70*C_Conc[6,5])+(V7*n2C_Conc[:,5])
#Humic acid
C_Conc[7,6]=(V70*C_Conc[6,6])+(V7*n2C_Conc[:,6])
#Monosaccharides
C_Conc[7,7]=(V70*C_Conc[6,7])+(V7*n2C_Conc[:,7])
#Inorganic carbon
C_Conc[7,8]=(V70*C_Conc[6,8])+(V7*n2C_Conc[:,8])
#Lipids
C_Conc[7,9]=(V70*C_Conc[6,9])+(V7*n2C_Conc[:,9])
#Pigments
C_Conc[7,10]=(V70*C_Conc[6,10])+(V7*n2C_Conc[:,10])
#Porpherin
C_Conc[7,11]=(V70*C_Conc[6,11])+(V7*n2C_Conc[:,11])
#CDOM
#C_Conc[:,12]=0
#VPhenol
C_Conc[7,13]=(V70*C_Conc[6,13])+(V7*n2C_Conc[:,13])
#SPhenol
C_Conc[7,14]=(V70*C_Conc[6,14])+(V7*n2C_Conc[:,14])

#mixing_node_8
#Protein
C_Conc[8,0]=(V80*C_Conc[7,0])+(V8*n3C_Conc[:,0])
#Poltpep
C_Conc[8,1]=(V80*C_Conc[7,1])+(V8*n3C_Conc[:,1])
#Amino_acid
C_Conc[8,2]=(V80*C_Conc[7,2])+(V8*n3C_Conc[:,2])
#Heteroploycondensate
C_Conc[8,3]=(V80*C_Conc[7,3])+(V8*n3C_Conc[:,3])
print(C_Conc[8,3])
#Starch_poly
C_Conc[8,4]=(V80*C_Conc[7,4])+(V8*n3C_Conc[:,4])
#Oligosaccharides
C_Conc[8,5]=(V80*C_Conc[7,5])+(V8*n3C_Conc[:,5])
#Humic acid
C_Conc[8,6]=(V80*C_Conc[7,6])+(V8*n3C_Conc[:,6])
#Monosaccharides
C_Conc[8,7]=(V80*C_Conc[7,7])+(V8*n3C_Conc[:,7])
#Inorganic carbon
C_Conc[8,8]=(V80*C_Conc[7,8])+(V8*n3C_Conc[:,8])
#Lipids
C_Conc[8,9]=(V80*C_Conc[7,9])+(V8*n3C_Conc[:,9])
#Pigments
C_Conc[8,10]=(V80*C_Conc[7,10])+(V8*n3C_Conc[:,10])
#Porpherin
C_Conc[8,11]=(V80*C_Conc[7,11])+(V8*n3C_Conc[:,11])
#CDOM
#C_Conc[:,12]=0
#VPhenol
C_Conc[8,13]=(V80*C_Conc[7,13])+(V8*n3C_Conc[:,13])
#SPhenol
C_Conc[8,14]=(V80*C_Conc[7,14])+(V8*n3C_Conc[:,14])

#mixing_node_9
#Protein
C_Conc[9,0]=(V90*C_Conc[8,0])+(V9*n3C_Conc[:,0])
#Poltpep
C_Conc[9,1]=(V90*C_Conc[8,1])+(V9*n3C_Conc[:,1])
#Amino_acid
C_Conc[9,2]=(V90*C_Conc[8,2])+(V9*n3C_Conc[:,2])
#Heteroploycondensate
C_Conc[9,3]=(V90*C_Conc[8,3])+(V9*n3C_Conc[:,3])
#Starch_poly
C_Conc[9,4]=(V90*C_Conc[8,4])+(V9*n3C_Conc[:,4])
#Oligosaccharides
C_Conc[9,5]=(V90*C_Conc[8,5])+(V9*n3C_Conc[:,5])
#Humic acid
C_Conc[9,6]=(V90*C_Conc[8,6])+(V9*n3C_Conc[:,6])
#Monosaccharides
C_Conc[9,7]=(V90*C_Conc[8,7])+(V9*n3C_Conc[:,7])
#Inorganic carbon
C_Conc[9,8]=(V90*C_Conc[8,8])+(V9*n3C_Conc[:,8])
#Lipids
C_Conc[9,9]=(V90*C_Conc[8,9])+(V9*n3C_Conc[:,9])
#Pigments
C_Conc[9,10]=(V90*C_Conc[8,10])+(V9*n3C_Conc[:,10])
#Porpherin
C_Conc[9,11]=(V90*C_Conc[8,11])+(V9*n3C_Conc[:,11])
#CDOM
#C_Conc[:,12]=0
#VPhenol
C_Conc[9,13]=(V90*C_Conc[8,13])+(V9*n3C_Conc[:,13])
#SPhenol
C_Conc[9,14]=(V90*C_Conc[8,14])+(V9*n3C_Conc[:,14])

#mixing_node_10
#Protein
C_Conc[10,0]=(V100*C_Conc[9,0])+(V10*n3C_Conc[:,0])
#Poltpep
C_Conc[10,1]=(V100*C_Conc[9,1])+(V10*n3C_Conc[:,1])
#Amino_acid
C_Conc[10,2]=(V100*C_Conc[9,2])+(V10*n3C_Conc[:,2])
#Heteroploycondensate
C_Conc[10,3]=(V100*C_Conc[9,3])+(V10*n3C_Conc[:,3])
#Starch_poly
C_Conc[10,4]=(V100*C_Conc[9,4])+(V10*n3C_Conc[:,4])
#Oligosaccharides
C_Conc[10,5]=(V100*C_Conc[9,5])+(V10*n3C_Conc[:,5])
#Humic acid
C_Conc[10,6]=(V100*C_Conc[9,6])+(V10*n3C_Conc[:,6])
#Monosaccharides
C_Conc[10,7]=(V100*C_Conc[9,7])+(V10*n3C_Conc[:,7])
#Inorganic carbon
C_Conc[10,8]=(V100*C_Conc[9,8])+(V10*n3C_Conc[:,8])
#Lipids
C_Conc[10,9]=(V100*C_Conc[9,9])+(V10*n3C_Conc[:,9])
#Pigments
C_Conc[10,10]=(V100*C_Conc[9,10])+(V10*n3C_Conc[:,10])
#Porpherin
C_Conc[10,11]=(V100*C_Conc[9,11])+(V10*n3C_Conc[:,11])
#CDOM
#C_Conc[:,12]=0
#VPhenol
C_Conc[10,13]=(V100*C_Conc[9,13])+(V10*n3C_Conc[:,13])
#SPhenol
C_Conc[10,14]=(V100*C_Conc[9,14])+(V10*n3C_Conc[:,14])

#mixing_node_11
#Protein
C_Conc[11,0]=(V110*C_Conc[10,0])+(V11*n3C_Conc[:,0])
#Poltpep
C_Conc[11,1]=(V110*C_Conc[10,1])+(V11*n3C_Conc[:,1])
#Amino_acid
C_Conc[11,2]=(V110*C_Conc[10,2])+(V11*n3C_Conc[:,2])
#Heteroploycondensate
C_Conc[11,3]=(V110*C_Conc[10,3])+(V11*n3C_Conc[:,3])
#Starch_poly
C_Conc[11,4]=(V110*C_Conc[10,4])+(V11*n3C_Conc[:,4])
#Oligosaccharides
C_Conc[11,5]=(V110*C_Conc[10,5])+(V11*n3C_Conc[:,5])
#Humic acid
C_Conc[11,6]=(V110*C_Conc[10,6])+(V11*n3C_Conc[:,6])
#Monosaccharides
C_Conc[11,7]=(V110*C_Conc[10,7])+(V11*n3C_Conc[:,7])
#Inorganic carbon
C_Conc[11,8]=(V110*C_Conc[10,8])+(V11*n3C_Conc[:,8])
#Lipids
C_Conc[11,9]=(V110*C_Conc[10,9])+(V11*n3C_Conc[:,9])
#Pigments
C_Conc[11,10]=(V110*C_Conc[10,10])+(V11*n3C_Conc[:,10])
#Porpherin
C_Conc[11,11]=(V110*C_Conc[10,11])+(V11*n3C_Conc[:,11])
#CDOM
#C_Conc[:,12]=0
#VPhenol
C_Conc[11,13]=(V110*C_Conc[10,13])+(V11*n3C_Conc[:,13])
#SPhenol
C_Conc[11,14]=(V110*C_Conc[10,14])+(V11*n3C_Conc[:,14])



#mixing_no#Assuming the chemistry is verl low , in that case the mixing node concentration will equal to rivermouth concentration
#Protein
C_Conc[12,0]=C_Conc[11,0]
#Poltpep
C_Conc[12,1]=C_Conc[11,1]
#Amino_acid
C_Conc[12,2]=C_Conc[11,2]
#Heteroploycondensate
C_Conc[12,3]=C_Conc[11,3]
print(C_Conc[10,3])
#Starch_poly
C_Conc[12,4]=C_Conc[11,4]
#Oligosaccharides
C_Conc[12,5]=C_Conc[11,5]
#Humic acid
C_Conc[12,6]=C_Conc[11,6]
#Monosaccharides
C_Conc[12,7]=C_Conc[11,7]
#Inorganic carbon
C_Conc[12,8]=C_Conc[11,8]
#Lipids
C_Conc[12,9]=C_Conc[11,9]
#Pigments
C_Conc[12,10]=C_Conc[11,10]
#Porpherin
C_Conc[12,11]=C_Conc[11,11]
#CDOM
#C_Conc[:,12]=0
#VPhenol
C_Conc[12,13]=C_Conc[11,13]
#SPhenol
C_Conc[12,14]=C_Conc[11,14]

#list of concentration of the species

#print("* Concentration list of","Chem0","in micromolar")
#Conc_Che = rv.C_Conc[:,0]
print(C_Conc[:,0])
Conc_Che = np.log10(C_Conc[:,0])
Conc_Che = Conc_Che.tolist()
#print(Conc_Che)
#1
Conc_Che1 = np.log10(C_Conc[:,1])
Conc_Che1 = Conc_Che1.tolist()
#2
Conc_Che2 = np.log10(C_Conc[:,2])
Conc_Che2 = Conc_Che2.tolist()
#3
Conc_Che3 = np.log10(C_Conc[:,3])
Conc_Che3 = Conc_Che3.tolist()
#4
Conc_Che4 = np.log10(C_Conc[:,4])
Conc_Che4 = Conc_Che4.tolist()
#5
Conc_Che5 = np.log10(C_Conc[:,5])
Conc_Che5 = Conc_Che5.tolist()
#6
Conc_Che6 = np.log10(C_Conc[:,6])
Conc_Che6 = Conc_Che6.tolist()
#7
Conc_Che7 = np.log10(C_Conc[:,7])
Conc_Che7 = Conc_Che7.tolist()
#8
Conc_Che8 = np.log10(C_Conc[:,8])
Conc_Che8 = Conc_Che8.tolist()
#9
Conc_Che9 = np.log10(C_Conc[:,9])
Conc_Che9 = Conc_Che9.tolist()
#10
Conc_Che10 = np.log10(C_Conc[:,10])
Conc_Che10 = Conc_Che10.tolist()
#11
Conc_Che11 = np.log10(C_Conc[:,11])
Conc_Che11 = Conc_Che11.tolist()
#13
Conc_Che13 = np.log10(C_Conc[:,13])
Conc_Che13 = Conc_Che13.tolist()
#14
Conc_Che14 = np.log10(C_Conc[:,14])
Conc_Che14 = Conc_Che14.tolist()
#Sum of chemical group
#Total DOC**HPC Duplicating ????
TDOC =C_Conc[:,0]+C_Conc[:,1]+C_Conc[:,2]+C_Conc[:,3]+C_Conc[:,4]+C_Conc[:,5]+C_Conc[:,6]+C_Conc[:,7]+C_Conc[:,9]+C_Conc[:,10]+C_Conc[:,11]+C_Conc[:,13]+C_Conc[:,14] 
Conc_Che30 = np.log10(TDOC)
Conc_Che30 = Conc_Che30.tolist()
#Humic+ HPC
Hu = C_Conc[:,3]+C_Conc[:,6]
Conc_Che31 = np.log10(Hu)
Conc_Che31 = Conc_Che31.tolist()
#CDOM
CDOM =(0.1*C_Conc[:,0])+(0.1*C_Conc[:,3])+C_Conc[:,10]+(1/3)*C_Conc[:,13]+(1/3)*C_Conc[:,14] 
Conc_Che32 = np.log10(CDOM)
Conc_Che32 = Conc_Che32.tolist()


dist = [0,250,500,750,1000,1250,1500,1750,2000,2250,2500,2700,3000]

fig = plt.figure()
fig.patch.set_facecolor('grey')
fig.suptitle("Dilution/Mixing_HighSolute(x3)+Tundra(x3)")
plt.axis([0,3000,-7,4])
plt.plot(dist,Conc_Che30,color='orchid',label="TDOC")
#plt.plot(dist,Conc_Che31,color='black',label="Humic+HPC")
plt.plot(dist,Conc_Che6,color='lightcoral',label="Humic acid")
plt.plot(dist,Conc_Che3,color='blue',label="HPC")
plt.plot(dist,Conc_Che,color='grey',label="Protein")
#plt.plot(dist,Conc_Che1,color='red',label="polypep")
#plt.plot(dist,Conc_Che2,color='black',label="amino_acid")
plt.plot(dist,Conc_Che4,color='red',label="Polysaccharide")
#plt.plot(dist,Conc_Che5,color='black',label="Oligosaccharide")
#plt.plot(dist,Conc_Che7,color='lightcyan',label="Monosaccharide")
#plt.plot(dist,Conc_Che8,color='olive',label="Inorganic_carbon")
plt.plot(dist,Conc_Che9,color='olive',label="Lipids")
plt.plot(dist,Conc_Che10,color='plum',label="Pigments")
#plt.plot(dist,Conc_Che11,color='orchid',label="Porpherin")
plt.plot(dist,Conc_Che13,color='orange',label="VPhenol")
#plt.plot(dist,Conc_Che14,color='moccasin',label="SPhenol")
plt.plot(dist,Conc_Che14,color='darkcyan',label="SPhenol")
plt.legend(fontsize='x-small')
plt.ylabel("Log[C]/micromolar")
plt.xlabel("Distance/km")
plt.grid(True)
plt.tight_layout()
plt.show()

        
#fig,axes= plt.subplots(nrows=2, ncols=2)
fig = plt.figure()
fig.patch.set_facecolor('grey')
fig.suptitle("Dilution/Mixing_HighSolute(x3)+Tundra(x3)")
ax  = fig.add_subplot(2,2,1)
ax1 = fig.add_subplot(2,2,2)
ax2 = fig.add_subplot(2,2,3)
ax3 = fig.add_subplot(2,2,4)

#ax.set_facecolor('xkcd:salmon')
ax.axis([0,3000,-4,4])
ax.plot(dist,Conc_Che30,'k-',label='(DOC)')
ax.plot(dist,Conc_Che6,color='grey',label='(Humic)')
ax.plot(dist,Conc_Che3,color='y',label='(HPC)')
ax.plot(dist,Conc_Che,'r-',label='Protein')
ax.set_title('Protein group')
ax.set_ylabel("Log[C]/micromolar")
#ax.set_xlabel("Distance/km")
ax.grid(True)
#ax.legend(loc=9, fontsize='x-small')
ax.legend(fontsize='x-small')

ax1.axis([0,3000,-4,4])
ax1.plot(dist,Conc_Che30,'k-',label='(DOC)')
ax1.plot(dist,Conc_Che6,color='grey',label='(Humic)')
ax1.plot(dist,Conc_Che3,color='y',label='(HPC)')
ax1.plot(dist,Conc_Che4,'r-',label='Polysaccharide')
#ax1.plot(dist,Conc_Che5,'k-',label='Oligosaccharides')
#ax1.plot(dist,Conc_Che2,'k--',label='Monosaccharides')
ax1.set_title('Polysachcharide group')
#ax1.set_ylabel("Log[C]/micromolar")
#ax1.set_xlabel("Distance/km")
ax1.grid(True)
ax1.legend(fontsize='x-small')

ax2.axis([0,3000,-7,4])
ax2.plot(dist,Conc_Che30,'k-',label='(DOC)')
ax2.plot(dist,Conc_Che6,color='grey',label='(Humic)')
ax2.plot(dist,Conc_Che3,color='y',label='(HPC)')
ax2.plot(dist,Conc_Che9,'r-',label='Lipid')
ax2.plot(dist,Conc_Che13,'b-',label='VPhenols')
ax2.plot(dist,Conc_Che14,'c-',label='SPhenols')
ax2.set_title('Lipids and phenols')
ax2.set_ylabel("Log[C]/micromolar")
ax2.set_xlabel("Distance/km")
ax2.grid(True)
ax2.legend(fontsize='x-small')

ax3.axis([0,3000,-5,4])
ax3.plot(dist,Conc_Che30,'k-',label='(DOC)')
ax3.plot(dist,Conc_Che6,color='grey',label='(Humic)')
ax3.plot(dist,Conc_Che3,color='y',label='(HPC)')
ax3.plot(dist,Conc_Che32,'r-',label='CDOM')
ax3.plot(dist,Conc_Che10,'b-',label='Pigments')
#ax3.plot(dist,Conc_Che11,'k--',label='Porphyrins')
ax3.set_title('Colored materials group')
#ax3.set_ylabel("Log[C]/micromolar")
ax3.set_xlabel("Distance/km")
ax3.grid(True)
ax3.legend(fontsize='x-small')

plt.tight_layout()
plt.show()
