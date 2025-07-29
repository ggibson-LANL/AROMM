import math
import numpy as np
import matplotlib.pyplot as plt
import ex_9_1 as e1
import ex_9_2 as e2
import ex_9_3 as e3
import ex_9_4 as e4
import ex_9_5 as e5
import ex_9_6 as e6
import ex_9_7 as e7

#Protein
e1.Conc_Che
e4.Conc_Che
e7.Conc_Che
#Poltpep
e1.Conc_Che1
e4.Conc_Che1
e7.Conc_Che1
#Amino_acid
e1.Conc_Che2
e4.Conc_Che2
e7.Conc_Che2
#Heteroploycondensate
e1.Conc_Che3
e4.Conc_Che3
e7.Conc_Che3
#Starch_poly
e1.Conc_Che4
e4.Conc_Che4
e7.Conc_Che4
#Oligosaccharides
e1.Conc_Che5
e4.Conc_Che5
e7.Conc_Che5
#Humic acid
e1.Conc_Che6
e4.Conc_Che6
e7.Conc_Che6
#Monosaccharides
e1.Conc_Che7
e4.Conc_Che7
e7.Conc_Che7
#Inorganic carbon
e1.Conc_Che8
e4.Conc_Che8
e7.Conc_Che8
#Lipids
e1.Conc_Che9
e4.Conc_Che9
e7.Conc_Che9
#Pigments
e1.Conc_Che10
e4.Conc_Che10
e7.Conc_Che10
#Porpherin
e1.Conc_Che11
e4.Conc_Che11
e7.Conc_Che11
#CDOM
#C_Conc[:,12]=0
#VPhenol
e1.Conc_Che13
e4.Conc_Che13
e7.Conc_Che13
#SPhenol
e1.Conc_Che14
e4.Conc_Che14
e7.Conc_Che14
#Total DOC
e1.Conc_Che30
print('**',e1.Conc_Che30)
e4.Conc_Che30
e7.Conc_Che30 
#Humic+ HPC
e1.Conc_Che31
e4.Conc_Che31
e7.Conc_Che31

dist = [0,250,500,750,1000,1250,1500,1750,2000,2250,2500,2700,3000]

fig = plt.figure()
fig.patch.set_facecolor('grey')
fig.suptitle("Dilution/Mixing")
plt.axis([0,3000,1,4])
plt.plot(dist,e1.Conc_Che30,color='k',label="e1.TDOC")
plt.plot(dist,e2.Conc_Che30,color='orange',label="e2.TDOC")
plt.plot(dist,e3.Conc_Che30,color='orchid',label="e3.TDOC")
plt.plot(dist,e4.Conc_Che30,color='r',label="e4.TDOC")
plt.plot(dist,e5.Conc_Che30,color='darkcyan',label="e5.TDOC")
plt.plot(dist,e6.Conc_Che30,color='blue',label="e6.TDOC")
plt.plot(dist,e7.Conc_Che30,color='grey',label="e7.TDOC")
#plt.plot(dist,Conc_Che31,color='black',label="Humic+HPC")
#plt.plot(dist,Conc_Che6,color='lightcoral',label="Humic acid")
#plt.plot(dist,Conc_Che3,color='blue',label="HPC")
#plt.plot(dist,Conc_Che,color='grey',label="Protein")
#plt.plot(dist,Conc_Che1,color='red',label="polypep")
#plt.plot(dist,Conc_Che2,color='black',label="amino_acid")
#plt.plot(dist,Conc_Che4,color='red',label="Polysaccharide")
#plt.plot(dist,Conc_Che5,color='black',label="Oligosaccharide")
#plt.plot(dist,Conc_Che7,color='lightcyan',label="Monosaccharide")
#plt.plot(dist,Conc_Che8,color='olive',label="Inorganic_carbon")
#plt.plot(dist,Conc_Che9,color='olive',label="Lipids")
#plt.plot(dist,Conc_Che10,color='plum',label="Pigments")
#plt.plot(dist,Conc_Che11,color='orchid',label="Porpherin")
#plt.plot(dist,Conc_Che13,color='orange',label="VPhenol")
#plt.plot(dist,Conc_Che14,color='moccasin',label="SPhenol")
#plt.plot(dist,Conc_Che14,color='darkcyan',label="SPhenol")
plt.legend(fontsize='x-small')
plt.ylabel("Log[C]/micromolar")
plt.xlabel("Distance/km")
plt.grid(True)
plt.tight_layout()
plt.show()

#fig,axes= plt.subplots(nrows=2, ncols=2)
fig = plt.figure()
fig.patch.set_facecolor('grey')
fig.suptitle("Dilution/Mixing")
ax  = fig.add_subplot(1,2,1)
ax9 = fig.add_subplot(1,2,2)
#ax.set_facecolor('xkcd:salmon')
ax.axis([0,3000,-1,3])
ax.plot(dist,e1.Conc_Che,'k-',label='Ref')
ax.plot(dist,e2.Conc_Che,'g:',label='Tundra(1/3)')
ax.plot(dist,e5.Conc_Che,'g--',label='Tundra(3)')
ax.plot(dist,e3.Conc_Che,'b:',label='HighSolute(1/3)')
ax.plot(dist,e6.Conc_Che,'b--',label='HighSolute(3)')
ax.plot(dist,e4.Conc_Che,'r:',label='Tundra(1/3)+HighSolute(1/3)')
ax.plot(dist,e7.Conc_Che,'r--',label='Tundra(3)+HighSolute(3)')
ax.set_title('Protein')
ax.set_ylabel("Log[C]/micromolar")
ax.set_xlabel("Distance/km")
ax.grid(True)
#ax.legend(loc=9, fontsize='x-small')
ax.legend(fontsize='x-small')

ax9.axis([0,3000,-7,3])
ax9.plot(dist,e1.Conc_Che14,'k-',label='Ref')
ax9.plot(dist,e2.Conc_Che14,'g:',label='Tundra(1/3)')
ax9.plot(dist,e5.Conc_Che14,'g--',label='Tundra(3)')
ax9.plot(dist,e3.Conc_Che14,'b:',label='HighSolute(1/3)')
ax9.plot(dist,e6.Conc_Che14,'b--',label='HighSolute(3)')
ax9.plot(dist,e4.Conc_Che14,'r:',label='Tundra(1/3)+HighSolute(1/3)')
ax9.plot(dist,e7.Conc_Che14,'r--',label='Tundra(3)+HighSolute(3)')
ax9.set_title('S-Phenol')
#ax3.set_ylabel("Log[C]/micromolar")
ax9.set_xlabel("Distance/km")
ax9.grid(True)
ax9.legend(fontsize='x-small')


plt.tight_layout()
plt.show()
