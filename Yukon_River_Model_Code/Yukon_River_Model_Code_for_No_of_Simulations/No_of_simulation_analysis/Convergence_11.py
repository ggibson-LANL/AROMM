import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties

# Updated data
x = [
    500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000,
    5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500, 10000,
    10500, 11000, 11500, 12000, 12500, 13000, 13500, 14000, 14500, 15000,
    15500, 16000, 16500, 17000, 17500, 18000, 18500, 19000, 19500, 20000,
    20500, 21000, 21500, 22000, 22500, 23000, 23500, 24000, 24500, 25000
    ]
#print(len(x))
y = [
    82.97, 90.91, 96.27, 99.62, 95.10, 93.73, 95.16, 95.74, 94.24, 95.89,
    91.93, 94.60, 95.04, 94.28, 93.05, 94.07, 94.40, 92.34, 91.16, 94.29,
    94.33, 93.74, 93.00, 92.76, 92.62, 94.27, 93.98, 92.98, 93.76, 94.79,
    93.92, 93.29, 94.15, 93.68, 92.02, 94.74, 93.48, 93.82, 93.38, 93.24,
    92.82, 94.76, 94.89, 94.02, 92.57, 94.72, 93.99, 92.97, 92.47, 93.59
]
print(len(y))

# Convert to pandas DataFrame
data = pd.DataFrame({'x': x, 'y': y})
print()

grad = []
grad1 = []
for i in range(len(x)-1):
    #print(i)
    g =  (data["y"][i+1]-data["y"][i]) /  (data["x"][i+1]-data["x"][i])
    g1 = abs(g)
    grad.append(g)
    grad1.append(g1)
    #print(g)
    
print(grad1)    
print(len(grad1))
print(sum(grad1))
print(sum(grad1)/49)

grad1.insert(0,float('nan'))
#my_list.insert(0, value_to_add)

data1 = pd.DataFrame({'x1':data['x'], 'y1': grad1})
print(data1)
print(data1.shape)
print(data1['x1'][49])
print(data1['x1'][5])

# Determine the convergence point
threshold_grad = 0.0022# You can adjust this threshold based on your data
convergence_point = data[data1['y1'] < threshold_grad].iloc[0]['x'] if not data[data1['y1'] < threshold_grad].empty else None

# Create FontProperties object for the label
font_properties = FontProperties(family='serif', style='italic', weight='normal', size=12)


# Plotting the data in grayscale
fig = plt.figure(figsize=(10, 6))
#plt.subplots(2,1)


#plt.axvline(x=2000, color='darkgrey', linestyle='--', label= 'x = 2000')  

# Vertical line for convergence point
if convergence_point:
    plt.axvline(x=convergence_point, color='k', linestyle='--', label=f'Convergence Point: {convergence_point}')

        
# Plot original data
#plt.plot(data['x'], data['y'], marker='o', color='black', label='Original Data')
plt.plot(data['x'], data['y'], marker='o', color='black')
#plt.plot(data['x'][1:],grad,  linestyle='-', color='darkgrey', label='Original Data')
# Updated axis labels
#plt.title("DOC_micromolar vs. X with Moving Average and Std Dev (Grayscale)")
plt.xlabel("No of Simulations",fontstyle = "italic",fontweight = 'bold',fontsize = 14,family = 'Times New Roman')  # Changed X-axis label
plt.xticks(np.arange(0, 25001, 2000),fontstyle = "normal",fontweight = 'bold',fontsize = 14,family = 'Times New Roman')  # Changed X-axis label
#plt.xlabel("Concentration (µM)")  # Changed X-axis label
plt.ylabel(r"$ \sigma^{2}$ of Output [DOC]", fontstyle = 'italic', fontweight='bold',fontsize ='14',family ='Times New Roman')  # Changed Y-axis label
plt.yticks(fontstyle = "normal",fontweight = 'bold',fontsize = 14,family = 'Times New Roman')  # Changed X-axis label
# Add legend with custom font properties
plt.legend(prop=font_properties,edgecolor='none')
plt.show()

# Plotting the data in grayscale
fig = plt.figure(figsize=(10, 6))

#plt.axvline(x=2000, color='darkgrey', linestyle='--', label= 'x = 2000')  
# Vertical line for convergence point
if convergence_point:
    plt.axvline(x=convergence_point, color='k', linestyle='--', label=f'Convergence Point: {convergence_point}')
         

# Add a horizontal line at y = 0.5
plt.plot([x[4],x[49]],[0.0022,0.0022], color='grey', linestyle='-.', label='Average = 0.0022')     
    
plt.plot(data['x'],grad1,  marker='o', color='k')
# Updated axis labels
#                                                                                plt.title("DOC_micromolar vs. X with Moving Average and Std Dev (Grayscale)")
plt.xlabel("No of Simulations",fontstyle = "italic",fontweight = 'bold',fontsize = 14,family = 'Times New Roman')  # Changed X-axis label
plt.xticks(np.arange(0, 25001, 2000),fontstyle = "normal",fontweight = 'bold',fontsize = 14,family = 'Times New Roman')  # Changed X-axis label
plt.ylabel("Absolute Value of the Gradient", fontstyle = 'italic', fontweight='bold',fontsize ='14',family ='Times New Roman')  # Changed Y-axis label
plt.yticks(fontstyle = "normal",fontweight = 'bold',fontsize = 14,family = 'Times New Roman')  # Changed X-axis label
# Add legend with custom font properties
plt.legend(prop=font_properties,edgecolor='none')
#plt.grid(True)
plt.show()

