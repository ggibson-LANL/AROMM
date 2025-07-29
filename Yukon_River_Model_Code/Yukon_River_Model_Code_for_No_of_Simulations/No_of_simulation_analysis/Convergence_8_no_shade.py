import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

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
#print(len(y))
# Convert to pandas DataFrame
data = pd.DataFrame({'x': x, 'y': y})

# Calculate moving average and moving standard deviation with a window size of 10
window_size = 12
data['moving_average'] = data['y'].rolling(window=window_size).mean()
data['moving_std'] = data['y'].rolling(window=window_size).std()

# Determine the convergence point
threshold_std = 2  # You can adjust this threshold based on your data
convergence_point = data[data['moving_std'] < threshold_std].iloc[0]['x'] if not data[data['moving_std'] < threshold_std].empty else None
#print(convergence_point)
# Plotting the data in grayscale
fig = plt.figure(figsize=(10, 6))

# Plot original data
plt.plot(data['x'], data['y'], marker='o', color='black', label='Original Data')

# Plot moving average
plt.plot(data['x'], data['moving_average'], color='grey', label='Moving Average')

# Fill between for moving standard deviation
#plt.fill_between(data['x'], data['moving_average'] - data['moving_std'], 
#                 data['moving_average'] + data['moving_std'], color='lightgrey', alpha=0.5, label='Moving Std Dev')

# Vertical line for convergence point
if convergence_point:
    plt.axvline(x=convergence_point, color='darkgrey', linestyle='--', label=f'Convergence Point: {convergence_point}')

# Updated axis labels
plt.title("DOC_micromolar vs. X with Moving Average and Std Dev (Grayscale)")
plt.xlabel("Concentration (µM)")  # Changed X-axis label
plt.ylabel("Measurement")  # Changed Y-axis label
plt.legend()
plt.grid(True)
plt.show()


fig.savefig(f"Figures/Figure_{window_size}_{convergence_point}_no_shade.jpeg")

