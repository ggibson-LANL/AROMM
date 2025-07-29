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
yr.Conc_Che17        #DOC



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

Conc = [yr.Conc_Che0,yr.Conc_Che1, yr.Conc_Che2, yr.Conc_Che3,yr.Conc_Che4, yr.Conc_Che5, yr.Conc_Che6, yr.Conc_Che7,
        yr.Conc_Che8,yr.Conc_Che9, yr.Conc_Che10, yr.Conc_Che11, yr.Conc_Che12, yr.Conc_Che13, yr.Conc_Che14, yr.Conc_Che15, yr.Conc_Che16,yr.Conc_Che17]

C_Name = ['Protein', 'Polypeptide', 'Amino acid', 'Heteropolycondensate', 'Polysaccharide', 'Oligosaccharide','Humic acid',
          'Monosaccharides','Inorganic carbon','Lipids','Pigments','Porphyrin','CDOM','V-Phenol','S-Phenol','C-Phenol',
          'Phenol', 'DOC']
Names = ['Tributary Velocity', 'Yukon Stem Velocity', 'Dilution Factor','DOC','Tau' , 'Chemical Fraction']          
          
#                                     ***  Create a figure and subplots    ***


#print(len(rv_variables[0]))
#print(len(yr.Conc_Che17))
# Loop through each subplot
for j in range(len(Conc)):
    ###Tributary Velocity
    # List of rv variables
    rv_variables = [
        rv.Tes_V1_samples,
        rv.Pel_V1_samples,
        rv.W_D_V1_samples,
        rv.Ste_V1_samples,
        rv.Por_V1_samples,
        rv.Tan_V1_samples,
        rv.Koy_V1_samples
    ]
    titles = [
        'Teslin', 'Pelly', 'White+Donjec', 'Stewart', 'Porcupine', 'Tanana', 'Koyukuk'
    ]
    fig, axs = plt.subplots(3, 3, figsize=(15, 15)) 
    fig.suptitle(f'{C_Name[j]} Vs Tributary Velocity',fontsize = 20,fontweight= "bold")
    fig.subplots_adjust(top=0.95)
    for i, ax in enumerate(axs.flat):

        # Check if the index is within the range of rv_variables
        if i < len(rv_variables):
            # Plot the scatter plot
            ax.scatter(rv_variables[i], Conc[j], color = 'black')
            #ax.set_title(f'{titles[i]}',fontsize = 25)
            ax.set_xlabel(f'{titles[i]}',fontsize = 15)
            #ax.set_ylabel('Y-axis')
            ax.tick_params(axis='x', labelsize=15)  # Change 15 to your desired font size
            ax.tick_params(axis='y', labelsize=15)  # Change 12 to your desired font size
           
        else:
            # Hide the subplot if there's no corresponding rv variable
            ax.axis('off')    
    # Adjust layout to prevent overlapping
    plt.tight_layout()
    # Show the plots
    plt.show()
    fig.savefig(f"output/analysis/Figures/ {C_Name[j]} Vs Tributary Velocity.png")
    
    ###Yukon Stem Velocity
    # List of rv variables
    rv_variables = [
        rv.Yuk_V100_samples,
        rv.Yuk_V1_samples,
        rv.Yuk_V2_samples,
        rv.Yuk_V3_samples,
        rv.Yuk_V4_samples,
        rv.Yuk_V5_samples,
        rv.Yuk_V6_samples,
        rv.Yuk_V7_samples,
        rv.Yuk_V8_samples,
    ]

    # List of titles
    titles = [
        'V100', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8'
    ]   
    fig, axs = plt.subplots(3, 3, figsize=(15, 15)) 
    fig.suptitle(f'{C_Name[j]} Vs Yukon Stem Velocity ',fontsize = 20,fontweight= "bold")
    fig.subplots_adjust(top=0.95)
    for i, ax in enumerate(axs.flat):

        # Check if the index is within the range of rv_variables
        if i < len(rv_variables):
            # Plot the scatter plot
            ax.scatter(rv_variables[i], Conc[j])
            #ax.set_title(f'{titles[i]}',fontsize = 25)
            ax.set_xlabel(f'{titles[i]}',fontsize = 15)
            #ax.set_ylabel('Y-axis')
            ax.tick_params(axis='x', labelsize=15)  # Change 15 to your desired font size
            ax.tick_params(axis='y', labelsize=15)  # Change 12 to your desired font size
           
        else:
            # Hide the subplot if there's no corresponding rv variable
            ax.axis('off')    
    # Adjust layout to prevent overlapping
    plt.tight_layout()
    # Show the plots
    plt.show()
    fig.savefig(f"output/analysis/Figures/ {C_Name[j]} Vs Yukon Stem Velocity.png")

    ###Dilution factor         
    # List of rv variables
    rv_variables = [
        rv.Tes_Df_samples,
        rv.Pel_Df_samples,
        rv.W_D_Df_samples,
        rv.Ste_Df_samples,
        rv.Por_Df_samples,
        rv.Tan_Df_samples,
        rv.Koy_Df_samples
    ]
    # List of titles
    titles = [
        'Teslin', 'Pelly', 'White+Donjec', 'Stewart', 'Porcupine', 'Tanana', 'Koyukuk'
    ]
    # Create a figure and subplots
    fig, axs = plt.subplots(3, 3, figsize=(15, 15))
    fig.suptitle(f'{C_Name[j]} Vs Dilution Fractions',fontsize = 20,fontweight= "bold")
    fig.subplots_adjust(top=0.95)
    # Loop through each subplot
    for i, ax in enumerate(axs.flat):
        # Check if the index is within the range of rv_variables
        if i < len(rv_variables):
            # Plot the scatter plot
            ax.scatter(rv_variables[i], Conc[j], color = 'blue')
            #ax.set_title(f'Scatter Plot Example ({titles[i]})')
            ax.set_xlabel(f'{titles[i]}', fontsize = 15)
            #ax.set_ylabel('Y-axis')
            ax.tick_params(axis='x', labelsize=15)  # Change 15 to your desired font size
            ax.tick_params(axis='y', labelsize=15)  # Change 12 to your desired font size
        else:
            # Hide the subplot if there's no corresponding rv variable
            ax.axis('off')

    # Adjust layout to prevent overlapping
    plt.tight_layout()

    # Show the plots
    plt.show()
    fig.savefig(f"output/analysis/Figures/ {C_Name[j]} Vs Dilution Fractions.png")
    
    ###DOC
    # List of rv variables
    rv_variables = [
        rv.Tes_DOC_samples,
        rv.Pel_DOC_samples,
        rv.W_D_DOC_samples,
        rv.Ste_DOC_samples,
        rv.Por_DOC_samples,
        rv.Tan_DOC_samples,
        rv.Koy_DOC_samples
    ]

    # List of titles
    titles = [
        'Teslin', 'Pelly', 'White+Donjec', 'Stewart', 'Porcupine', 'Tanana', 'Koyukuk'
    ]
    # Create a figure and subplots
    fig, axs = plt.subplots(3, 3, figsize=(15, 15))
    fig.suptitle(f'{C_Name[j]} Vs Input DOC',fontsize = 20,fontweight= "bold")
    fig.subplots_adjust(top=0.95)
    # Loop through each subplot
    for i, ax in enumerate(axs.flat):
        # Check if the index is within the range of rv_variables
        if i < len(rv_variables):
            # Plot the scatter plot
            ax.scatter(rv_variables[i], Conc[j], color = 'orange')
            #ax.set_title(f'Scatter Plot Example ({titles[i]})')
            ax.set_xlabel(f'{titles[i]}',fontsize = 15)
            #ax.set_ylabel('Y-axis')
            ax.tick_params(axis='x', labelsize=15)  # Change 15 to your desired font size
            ax.tick_params(axis='y', labelsize=15)  # Change 12 to your desired font size
        else:
            # Hide the subplot if there's no corresponding rv variable
            ax.axis('off')

    # Adjust layout to prevent overlapping
    plt.tight_layout()

    # Show the plots
    plt.show()
    fig.savefig(f"output/analysis/Figures/ {C_Name[j]} Vs Input DOC.png")
    
    ###Tau
    # List of rv variables
    rv_variables = [
        rv.prot_Tau_samples,
        rv.polypep_Tau_samples,
        rv.amino_acid_Tau_samples,
        rv.hpc_Tau_samples,
        rv.poly_Tau_samples,
        rv.olig_Tau_samples,
        rv.humic_acid_Tau_samples,
        rv.mono_Tau_samples,
        rv.ino_C_Tau_samples,
        rv.lip_Tau_samples,
        rv.pig_Tau_samples,
        rv.por_Tau_samples,
        rv.V_Phe_Tau_samples,
        rv.S_Phe_Tau_samples,
        rv.C_Phe_Tau_samples
    ]

    # List of titles
    titles = [
       'Protein', 'Polypeptide', 'Amino Acid', 'Heteroplycondenstate', 'Polysaccharide', 
        'Oligosaccharide', 'Humic Acid','Monosaccharide','Inorganic C','Lipids',
        'Pigments','Porphyrin','V-Phenol','S-Phenol','C-Phenol'
    ]
    # Create a figure and subplots
    fig, axs = plt.subplots(4, 4, figsize=(15, 15))
    fig.suptitle(f'{C_Name[j]} DOC Vs Tau',fontsize = 20,fontweight= "bold")
    fig.subplots_adjust(top=0.95)
    # Loop through each subplot
    for i, ax in enumerate(axs.flat):
        # Check if the index is within the range of rv_variables
        if i < len(rv_variables):
            # Plot the scatter plot
            ax.scatter(rv_variables[i], Conc[j], color = 'red')
            #ax.set_title(f'Scatter Plot Example ({titles[i]})')
            ax.set_xlabel(f'{titles[i]}',fontsize = 15)
            #ax.set_ylabel('Y-axis')
            ax.tick_params(axis='x', labelsize=15)  # Change 15 to your desired font size
            ax.tick_params(axis='y', labelsize=15)  # Change 12 to your desired font size
        else:
            # Hide the subplot if there's no corresponding rv variable
            ax.axis('off')

    # Adjust layout to prevent overlapping
    plt.tight_layout()

    # Show the plots
    plt.show()
    fig.savefig(f"output/analysis/Figures/ {C_Name[j]} Vs Tau.png")
    
    #Chemical fraction
    # Create a figure and subplots
    # List of rv variables
    rv_variables = [
        rv.prot_f_samples,
        rv.hpc_f_samples,
        rv.poly_f_samples,
        rv.lip_f_samples,
        rv.pig_f_samples,
        rv.V_Phe_f_samples,
        rv.S_Phe_f_samples,
        rv.C_Phe_f_samples
    ]

    # List of titles
    titles = [
       'Protein', 'Heteroplycondenstate', 'Polysaccharide', 'Lipids',
        'Pigments','V-Phenol','S-Phenol','C-Phenol'
    ]
    fig, axs = plt.subplots(4, 2, figsize=(15, 15))
    fig.suptitle(f'{C_Name[j]} Vs Chemical fraction',fontsize = 20,fontweight= "bold")
    fig.subplots_adjust(top=0.95)
    # Loop through each subplot
    for i, ax in enumerate(axs.flat):
        # Check if the index is within the range of rv_variables
        if i < len(rv_variables):
            # Plot the scatter plot
            ax.scatter(rv_variables[i], Conc[j], color = 'green')
            #ax.set_title(f'Scatter Plot Example ({titles[i]})')
            ax.set_xlabel(f'{titles[i]}',fontsize = 15)
            #ax.set_ylabel('Y-axis')
            ax.tick_params(axis='x', labelsize=15)  # Change 15 to your desired font size
            ax.tick_params(axis='y', labelsize=15)  # Change 12 to your desired font size
        else:
            # Hide the subplot if there's no corresponding rv variable
            ax.axis('off')

    # Adjust layout to prevent overlapping
    plt.tight_layout()

    # Show the plots
    plt.show()
    fig.savefig(f"output/analysis/Figures/ {C_Name[j]} Vs Chemical fraction.png")
    
    
    
        
    
