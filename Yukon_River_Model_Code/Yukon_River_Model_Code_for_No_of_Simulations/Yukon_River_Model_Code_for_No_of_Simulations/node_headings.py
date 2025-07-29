

     #node_1=Pelly_river
     print("  ")
     print("node_1:Pelly_river")
     print("j_value",j)
     print("node_1,stem:trib mixing")
     
     #Mixing
     C_Conc[j,:] = (1-rv.Pel_Df_samples[i])*C_Conc[j,:] + rv.Pel_Df_samples[i]*nd1rv.Chemicals_[i,:]  
     
     #node_2=White+Donjec_river
     print("  ")
     print("node_2:White+Donjec_river")
     print("j_value",j)
     print("node_2,stem:trib mixing")
 
     #Mixing
     C_Conc[j,:] = (1-rv.W_D_Df_samples[i])*C_Conc[j,:] + rv.W_D_Df_samples[i]*nd2rv.Chemicals_[i,:]
     
     #node_3=Stewart_river
     print("  ")
     print("node_3:Stewart_river")
     print("j_value",j)
     print("node_3,stem:trib mixing")
     
     #Mixing
     C_Conc[j,:] = (1-rv.Ste_Df_samples[i])*C_Conc[j,:] + rv.Ste_Df_samples[i]*nd3rv.Chemicals_[i,:]  

     #node_4=Porcupine_river
     print("  ")
     print("node_4:Porcupine_river")
     print("j_value",j)
     print("node_4,stem:trib mixing")    

     #Mixing
     C_Conc[j,:] = (1-rv.Por_Df_samples[i])*C_Conc[j,:] + rv.Por_Df_samples[i]*nd4rv.Chemicals_[i,:] 
     
     #node_5=Tanana_river
     print("  ")
     print("node_5:Tanana_river")
     print("j_value",j)
     print("node_5,stem:trib mixing")

     #Mixing
     C_Conc[j,:] = (1-rv.Tan_Df_samples[i])*C_Conc[j,:] + rv.Tan_Df_samples[i]*nd5rv.Chemicals_[i,:] 
     
     #node_6=Koyukuk_river
     print("  ")
     print("node_6:Pelly_river")
     print("j_value",j)
     print("node_6,stem:trib mixing")
     
     #Mixing
     C_Conc[j,:] = (1-rv.Koy_Df_samples[i])*C_Conc[j,:] + rv.Koy_Df_samples[i]*nd6rv.Chemicals_[i,:] 
        














#node_1=Pelly_river
     print("  ")
     print("node_1:Pelly_river")
     print("j_value",j)
     print("node_1,stem:trib mixing")
     
     #Mixing
     C_Conc[j,0] = (1-rv.Pel_Df_samples[i])*C_Conc[j,0] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che0[i]
     C_Conc[j,1] = (1-rv.Pel_Df_samples[i])*C_Conc[j,1] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che1[i]
     C_Conc[j,2] = (1-rv.Pel_Df_samples[i])*C_Conc[j,2] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che2[i]
     C_Conc[j,3] = (1-rv.Pel_Df_samples[i])*C_Conc[j,3] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che3[i]
     C_Conc[j,4] = (1-rv.Pel_Df_samples[i])*C_Conc[j,4] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che4[i]  
     C_Conc[j,5] = (1-rv.Pel_Df_samples[i])*C_Conc[j,5] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che5[i]
     C_Conc[j,6] = (1-rv.Pel_Df_samples[i])*C_Conc[j,6] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che6[i]
     C_Conc[j,7] = (1-rv.Pel_Df_samples[i])*C_Conc[j,7] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che7[i]
     C_Conc[j,8] = (1-rv.Pel_Df_samples[i])*C_Conc[j,8] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che8[i]
     C_Conc[j,9] = (1-rv.Pel_Df_samples[i])*C_Conc[j,9] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che9[i]
     C_Conc[j,10] = (1-rv.Pel_Df_samples[i])*C_Conc[j,10] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che10[i]
     C_Conc[j,11] = (1-rv.Pel_Df_samples[i])*C_Conc[j,11] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che11[i]
     C_Conc[j,12] = (1-rv.Pel_Df_samples[i])*C_Conc[j,12] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che12[i]
     C_Conc[j,13] = (1-rv.Pel_Df_samples[i])*C_Conc[j,13] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che13[i]
     C_Conc[j,14] = (1-rv.Pel_Df_samples[i])*C_Conc[j,14] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che14[i]
     C_Conc[j,15] = (1-rv.Pel_Df_samples[i])*C_Conc[j,15] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che15[i]
     C_Conc[j,16] = (1-rv.Pel_Df_samples[i])*C_Conc[j,16] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che16[i]
     C_Conc[j,17] = (1-rv.Pel_Df_samples[i])*C_Conc[j,17] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che17[i]

     #node_2=White+Donjec_river
     print("  ")
     print("node_2:White+Donjec_river")
     print("j_value",j)
     print("node_2,stem:trib mixing")
     
     #Mixing
     C_Conc[j,0] = (1-rv.Pel_Df_samples[i])*C_Conc[j,0] + rv.Pel_Df_samples[i]*nd2rv.Conc_Che0[i]
     C_Conc[j,1] = (1-rv.Pel_Df_samples[i])*C_Conc[j,1] + rv.Pel_Df_samples[i]*nd2rv.Conc_Che1[i]
     C_Conc[j,2] = (1-rv.Pel_Df_samples[i])*C_Conc[j,2] + rv.Pel_Df_samples[i]*nd2rv.Conc_Che2[i]
     C_Conc[j,3] = (1-rv.Pel_Df_samples[i])*C_Conc[j,3] + rv.Pel_Df_samples[i]*nd2rv.Conc_Che3[i]
     C_Conc[j,4] = (1-rv.Pel_Df_samples[i])*C_Conc[j,4] + rv.Pel_Df_samples[i]*nd2rv.Conc_Che4[i]  
     C_Conc[j,5] = (1-rv.Pel_Df_samples[i])*C_Conc[j,5] + rv.Pel_Df_samples[i]*nd2rv.Conc_Che5[i]
     C_Conc[j,6] = (1-rv.Pel_Df_samples[i])*C_Conc[j,6] + rv.Pel_Df_samples[i]*nd2rv.Conc_Che6[i]
     C_Conc[j,7] = (1-rv.Pel_Df_samples[i])*C_Conc[j,7] + rv.Pel_Df_samples[i]*nd2rv.Conc_Che7[i]
     C_Conc[j,8] = (1-rv.Pel_Df_samples[i])*C_Conc[j,8] + rv.Pel_Df_samples[i]*nd2rv.Conc_Che8[i]
     C_Conc[j,9] = (1-rv.Pel_Df_samples[i])*C_Conc[j,9] + rv.Pel_Df_samples[i]*nd2rv.Conc_Che9[i]
     C_Conc[j,10] = (1-rv.Pel_Df_samples[i])*C_Conc[j,10] + rv.Pel_Df_samples[i]*nd2rv.Conc_Che10[i]
     C_Conc[j,11] = (1-rv.Pel_Df_samples[i])*C_Conc[j,11] + rv.Pel_Df_samples[i]*nd2rv.Conc_Che11[i]
     C_Conc[j,12] = (1-rv.Pel_Df_samples[i])*C_Conc[j,12] + rv.Pel_Df_samples[i]*nd2rv.Conc_Che12[i]
     C_Conc[j,13] = (1-rv.Pel_Df_samples[i])*C_Conc[j,13] + rv.Pel_Df_samples[i]*nd2rv.Conc_Che13[i]
     C_Conc[j,14] = (1-rv.Pel_Df_samples[i])*C_Conc[j,14] + rv.Pel_Df_samples[i]*nd2rv.Conc_Che14[i]
     C_Conc[j,15] = (1-rv.Pel_Df_samples[i])*C_Conc[j,15] + rv.Pel_Df_samples[i]*nd2rv.Conc_Che15[i]
     C_Conc[j,16] = (1-rv.Pel_Df_samples[i])*C_Conc[j,16] + rv.Pel_Df_samples[i]*nd2rv.Conc_Che16[i]
     C_Conc[j,17] = (1-rv.Pel_Df_samples[i])*C_Conc[j,17] + rv.Pel_Df_samples[i]*nd2rv.Conc_Che17[i]
     
     
     #node_3=Stewart_river
     print("  ")
     print("node_3:Stewart_river")
     print("j_value",j)
     print("node_3,stem:trib mixing")
     
     #Mixing
     C_Conc[j,0] = (1-rv.Pel_Df_samples[i])*C_Conc[j,0] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che0[i]
     C_Conc[j,1] = (1-rv.Pel_Df_samples[i])*C_Conc[j,1] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che1[i]
     C_Conc[j,2] = (1-rv.Pel_Df_samples[i])*C_Conc[j,2] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che2[i]
     C_Conc[j,3] = (1-rv.Pel_Df_samples[i])*C_Conc[j,3] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che3[i]
     C_Conc[j,4] = (1-rv.Pel_Df_samples[i])*C_Conc[j,4] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che4[i]  
     C_Conc[j,5] = (1-rv.Pel_Df_samples[i])*C_Conc[j,5] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che5[i]
     C_Conc[j,6] = (1-rv.Pel_Df_samples[i])*C_Conc[j,6] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che6[i]
     C_Conc[j,7] = (1-rv.Pel_Df_samples[i])*C_Conc[j,7] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che7[i]
     C_Conc[j,8] = (1-rv.Pel_Df_samples[i])*C_Conc[j,8] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che8[i]
     C_Conc[j,9] = (1-rv.Pel_Df_samples[i])*C_Conc[j,9] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che9[i]
     C_Conc[j,10] = (1-rv.Pel_Df_samples[i])*C_Conc[j,10] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che10[i]
     C_Conc[j,11] = (1-rv.Pel_Df_samples[i])*C_Conc[j,11] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che11[i]
     C_Conc[j,12] = (1-rv.Pel_Df_samples[i])*C_Conc[j,12] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che12[i]
     C_Conc[j,13] = (1-rv.Pel_Df_samples[i])*C_Conc[j,13] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che13[i]
     C_Conc[j,14] = (1-rv.Pel_Df_samples[i])*C_Conc[j,14] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che14[i]
     C_Conc[j,15] = (1-rv.Pel_Df_samples[i])*C_Conc[j,15] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che15[i]
     C_Conc[j,16] = (1-rv.Pel_Df_samples[i])*C_Conc[j,16] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che16[i]
     C_Conc[j,17] = (1-rv.Pel_Df_samples[i])*C_Conc[j,17] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che17[i]
     
     
     #node_4=Porcupine_river
     print("  ")
     print("node_4:Porcupine_river")
     print("j_value",j)
     print("node_4,stem:trib mixing")
     
     #Mixing
     C_Conc[j,0] = (1-rv.Pel_Df_samples[i])*C_Conc[j,0] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che0[i]
     C_Conc[j,1] = (1-rv.Pel_Df_samples[i])*C_Conc[j,1] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che1[i]
     C_Conc[j,2] = (1-rv.Pel_Df_samples[i])*C_Conc[j,2] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che2[i]
     C_Conc[j,3] = (1-rv.Pel_Df_samples[i])*C_Conc[j,3] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che3[i]
     C_Conc[j,4] = (1-rv.Pel_Df_samples[i])*C_Conc[j,4] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che4[i]  
     C_Conc[j,5] = (1-rv.Pel_Df_samples[i])*C_Conc[j,5] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che5[i]
     C_Conc[j,6] = (1-rv.Pel_Df_samples[i])*C_Conc[j,6] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che6[i]
     C_Conc[j,7] = (1-rv.Pel_Df_samples[i])*C_Conc[j,7] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che7[i]
     C_Conc[j,8] = (1-rv.Pel_Df_samples[i])*C_Conc[j,8] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che8[i]
     C_Conc[j,9] = (1-rv.Pel_Df_samples[i])*C_Conc[j,9] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che9[i]
     C_Conc[j,10] = (1-rv.Pel_Df_samples[i])*C_Conc[j,10] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che10[i]
     C_Conc[j,11] = (1-rv.Pel_Df_samples[i])*C_Conc[j,11] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che11[i]
     C_Conc[j,12] = (1-rv.Pel_Df_samples[i])*C_Conc[j,12] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che12[i]
     C_Conc[j,13] = (1-rv.Pel_Df_samples[i])*C_Conc[j,13] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che13[i]
     C_Conc[j,14] = (1-rv.Pel_Df_samples[i])*C_Conc[j,14] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che14[i]
     C_Conc[j,15] = (1-rv.Pel_Df_samples[i])*C_Conc[j,15] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che15[i]
     C_Conc[j,16] = (1-rv.Pel_Df_samples[i])*C_Conc[j,16] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che16[i]
     C_Conc[j,17] = (1-rv.Pel_Df_samples[i])*C_Conc[j,17] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che17[i]
     

     #node_5=Tanana_river
     print("  ")
     print("node_5:Tanana_river")
     print("j_value",j)
     print("node_5,stem:trib mixing")
     
     #Mixing
     C_Conc[j,0] = (1-rv.Pel_Df_samples[i])*C_Conc[j,0] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che0[i]
     C_Conc[j,1] = (1-rv.Pel_Df_samples[i])*C_Conc[j,1] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che1[i]
     C_Conc[j,2] = (1-rv.Pel_Df_samples[i])*C_Conc[j,2] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che2[i]
     C_Conc[j,3] = (1-rv.Pel_Df_samples[i])*C_Conc[j,3] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che3[i]
     C_Conc[j,4] = (1-rv.Pel_Df_samples[i])*C_Conc[j,4] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che4[i]  
     C_Conc[j,5] = (1-rv.Pel_Df_samples[i])*C_Conc[j,5] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che5[i]
     C_Conc[j,6] = (1-rv.Pel_Df_samples[i])*C_Conc[j,6] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che6[i]
     C_Conc[j,7] = (1-rv.Pel_Df_samples[i])*C_Conc[j,7] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che7[i]
     C_Conc[j,8] = (1-rv.Pel_Df_samples[i])*C_Conc[j,8] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che8[i]
     C_Conc[j,9] = (1-rv.Pel_Df_samples[i])*C_Conc[j,9] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che9[i]
     C_Conc[j,10] = (1-rv.Pel_Df_samples[i])*C_Conc[j,10] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che10[i]
     C_Conc[j,11] = (1-rv.Pel_Df_samples[i])*C_Conc[j,11] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che11[i]
     C_Conc[j,12] = (1-rv.Pel_Df_samples[i])*C_Conc[j,12] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che12[i]
     C_Conc[j,13] = (1-rv.Pel_Df_samples[i])*C_Conc[j,13] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che13[i]
     C_Conc[j,14] = (1-rv.Pel_Df_samples[i])*C_Conc[j,14] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che14[i]
     C_Conc[j,15] = (1-rv.Pel_Df_samples[i])*C_Conc[j,15] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che15[i]
     C_Conc[j,16] = (1-rv.Pel_Df_samples[i])*C_Conc[j,16] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che16[i]
     C_Conc[j,17] = (1-rv.Pel_Df_samples[i])*C_Conc[j,17] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che17[i]


     #node_6=Pelly_river
     print("  ")
     print("node_6:Pelly_river")
     print("j_value",j)
     print("node_6,stem:trib mixing")
     
     #Mixing
     C_Conc[j,0] = (1-rv.Pel_Df_samples[i])*C_Conc[j,0] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che0[i]
     C_Conc[j,1] = (1-rv.Pel_Df_samples[i])*C_Conc[j,1] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che1[i]
     C_Conc[j,2] = (1-rv.Pel_Df_samples[i])*C_Conc[j,2] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che2[i]
     C_Conc[j,3] = (1-rv.Pel_Df_samples[i])*C_Conc[j,3] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che3[i]
     C_Conc[j,4] = (1-rv.Pel_Df_samples[i])*C_Conc[j,4] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che4[i]  
     C_Conc[j,5] = (1-rv.Pel_Df_samples[i])*C_Conc[j,5] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che5[i]
     C_Conc[j,6] = (1-rv.Pel_Df_samples[i])*C_Conc[j,6] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che6[i]
     C_Conc[j,7] = (1-rv.Pel_Df_samples[i])*C_Conc[j,7] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che7[i]
     C_Conc[j,8] = (1-rv.Pel_Df_samples[i])*C_Conc[j,8] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che8[i]
     C_Conc[j,9] = (1-rv.Pel_Df_samples[i])*C_Conc[j,9] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che9[i]
     C_Conc[j,10] = (1-rv.Pel_Df_samples[i])*C_Conc[j,10] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che10[i]
     C_Conc[j,11] = (1-rv.Pel_Df_samples[i])*C_Conc[j,11] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che11[i]
     C_Conc[j,12] = (1-rv.Pel_Df_samples[i])*C_Conc[j,12] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che12[i]
     C_Conc[j,13] = (1-rv.Pel_Df_samples[i])*C_Conc[j,13] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che13[i]
     C_Conc[j,14] = (1-rv.Pel_Df_samples[i])*C_Conc[j,14] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che14[i]
     C_Conc[j,15] = (1-rv.Pel_Df_samples[i])*C_Conc[j,15] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che15[i]
     C_Conc[j,16] = (1-rv.Pel_Df_samples[i])*C_Conc[j,16] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che16[i]
     C_Conc[j,17] = (1-rv.Pel_Df_samples[i])*C_Conc[j,17] + rv.Pel_Df_samples[i]*nd1rv.Conc_Che17[i]
     