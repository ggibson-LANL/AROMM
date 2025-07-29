import numpy as np
import math

#Maximum no of time steps
imax=480
#print("imax",imax)

#cmax=no of chemical species
cmax= 30
#print("cmax",cmax)

#chemical species in an array
Tau = np.ones((imax,cmax))
#Tau[:,:] =0
Tau[:,0] =864000
Tau[:,1] =864000*(3/13)
Tau[:,2] =864000*(1/10)
Tau[:,3] =864000
Tau[:,4] =864000*3
Tau[:,5] =(10/11)*864000
Tau[:,6] =864000*3
Tau[:,7] =864000*3
Tau[:,8] =86400*3
Tau[:,9] =864000*(3/13)
Tau[:,10] =864000
Tau[:,11] =864000
#Tau[:,12] =864000
Tau[:,13] =86400*3
Tau[:,14] =86400*3
#Tau[:,15] =0
#Tau[:,16] =0
#Tau[:,17] =0
#Tau[:,18] =0
#Tau[:,19] =0
#Tau[:,20] =0
#Tau[:,21] =0
#Tau[:,22] =0
#Tau[:,23] =0
#Tau[:,24] =0
#Tau[:,25] =0
#Tau[:,26] =0
#Tau[:,27] =0
#Tau[:,28] =0
#Tau[:,29] =0


#print("Tau = ",Tau)


#Initial concentration of the species ia an array
#Concentrations in the micromolar
C_Conc = np.ones((imax,cmax))
#C_Conc[:,:]=0
#Protein
C_Conc[:,0]=15
#Poltpep
C_Conc[:,1]=0
#Amino acid
C_Conc[:,2]=0
#Heteroploycondensate
C_Conc[:,3]=255
#Starch/poly
C_Conc[:,4]=15
#Oligosaccharides
C_Conc[:,5]=0
#Humic acid
C_Conc[:,6]=0
#Monosaccharides
C_Conc[:,7]=0
#Inorganic carbon
C_Conc[:,8]=0
#Lipids
C_Conc[:,9]=15
#Pigments
C_Conc[:,10]=3
#Porpherin
C_Conc[:,11]=0
#CDOM
#C_Conc[:,12]=0
#VPhenol
C_Conc[:,13]=3
#SPhenol
C_Conc[:,14]=0
#TDOC
TDOC =C_Conc[:,0]+C_Conc[:,1]+C_Conc[:,2]+C_Conc[:,3]+C_Conc[:,4]+C_Conc[:,5]+C_Conc[:,6]+C_Conc[:,7]+C_Conc[:,9]+C_Conc[:,10]+C_Conc[:,11]+C_Conc[:,13]+C_Conc[:,14]
#Scaling
#Protein
C_Conc[:,0]=C_Conc[:,0]*(300/TDOC)
#Poltpep
C_Conc[:,1]=C_Conc[:,1]*(300/TDOC)
#Amino acid
C_Conc[:,2]=C_Conc[:,2]*(300/TDOC)
#Heteroploycondensate
C_Conc[:,3]=C_Conc[:,3]*(300/TDOC)
print(TDOC[0],"Primary")
print(C_Conc[0,3])
#Starch
C_Conc[:,4]=C_Conc[:,4]*(300/TDOC)
#Oligosaccharides
C_Conc[:,5]=C_Conc[:,5]*(300/TDOC)
#Humic acid
C_Conc[:,6]=C_Conc[:,6]*(300/TDOC)
#Monosaccharides
C_Conc[:,7]=C_Conc[:,7]*(300/TDOC)
#Inorganic carbon
C_Conc[:,8]=C_Conc[:,8]*(300/TDOC)
#Lipids
C_Conc[:,9]=C_Conc[:,9]*(300/TDOC)
#Pigments
C_Conc[:,10]=C_Conc[:,10]*(300/TDOC)
#Porpherin
C_Conc[:,11]=C_Conc[:,11]*(300/TDOC)
#CDOM
#C_Conc[:,12]=0
#VPhenol
C_Conc[:,13]=C_Conc[:,13]*(300/TDOC)
#SPhenol
C_Conc[:,14]=C_Conc[:,14]*(300/TDOC)
#
#C_Conc[:,15]=0
#C_Conc[:,16]=0
#C_Conc[:,17]=0
#C_Conc[:,18]=0
#C_Conc[:,19]=0
#C_Conc[:,20]=0
#C_Conc[:,21]=0
#C_Conc[:,22]=0
#C_Conc[:,23]=0
#C_Conc[:,24]=0
#C_Conc[:,25]=0
#C_Conc[:,26]=0
#C_Conc[:,27]=0
#C_Conc[:,28]=0
#C_Conc[:,29]=0

#print("C_Conc in micro molar",C_Conc)

#Production of the chemical species in an array
Prod = np.ones((imax,cmax))
#Prod[:,:] =0
Prod[:,0] =0
#Prod[0,1] =C_Conc[0,0]/Tau[0,0]
#Prod[:,2] =1000/100000
#Prod[:,3] =0
Prod[:,4] =0
#Prod[:,5] =0
#Prod[:,6] =0
#Prod[:,7] =0
#Prod[:,8] =0
Prod[:,9] =0
Prod[:,10] =0
#Prod[:,11] =0
#Prod[:,12] =0
Prod[:,13] =0
Prod[:,14] =0
#Prod[:,15] =0
#Prod[:,16] =0
#Prod[:,17] =0
#Prod[:,18] =0
#Prod[:,19] =0
#Prod[:,20] =0
#Prod[:,21] =0
#Prod[:,22] =0
#Prod[:,23] =0
#Prod[:,24] =0
#Prod[:,25] =0
#Prod[:,26] =0
#Prod[:,27] =0
#Prod[:,28] =0
#Prod[:,29] =0


#print("Prod = ",Prod)