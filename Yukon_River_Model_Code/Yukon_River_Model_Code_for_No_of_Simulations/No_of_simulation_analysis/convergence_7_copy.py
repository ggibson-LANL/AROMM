import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Given data
x = [
    100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
    1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000,
    2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000,
    3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000,
    4500, 5000,
    
    5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500, 10000,
    10500, 11000, 11500, 12000, 12500, 13000, 13500, 14000, 14500, 15000,
    15500, 16000, 16500, 17000, 17500, 18000, 18500, 19000, 19500, 20000,
    20500, 21000, 21500, 22000, 22500, 23000, 23500, 24000, 24500, 25000
]
y = [
    118.050, 107.534, 103.918, 102.719, 82.966, 93.305, 104.037, 90.703, 94.095, 90.908,
    93.377, 99.134, 90.140, 92.473, 96.273, 95.431, 94.526, 92.977, 93.372, 99.617,
    90.930, 92.208, 93.323, 95.463, 95.098, 94.918, 90.865, 95.678, 93.032, 93.729,
    94.363, 93.544, 95.244, 91.585, 95.161, 96.813, 93.913, 95.366, 95.110, 95.735,
    94.24, 95.89,91.93, 94.60, 95.04, 94.28, 93.05, 94.07, 94.40, 92.34, 91.16, 94.29,
    94.33, 93.74, 93.00, 92.76, 92.62, 94.27, 93.98, 92.98, 93.76, 94.79,
    93.92, 93.29, 94.15, 93.68, 92.02, 94.74, 93.48, 93.82, 93.38, 93.24,
    92.82, 94.76, 94.89, 94.02, 92.57, 94.72, 93.99, 92.97, 92.47, 93.59
    
]

# Convert to pandas DataFrame
data = pd.DataFrame({'x': x, 'y': y})

# Calculate moving average and moving standard deviation with a window size of 10
window_size = 10
data['moving_average'] = data['y'].rolling(window=window_size).mean()
data['moving_std'] = data['y'].rolling(window=window_size).std()

# Determine the convergence point
threshold_std = 2  # You can adjust this threshold based on your data
convergence_point = data[data['moving_std'] < threshold_std].iloc[0]['x'] if not data[data['moving_std'] < threshold_std].empty else None

# Plotting the data
plt.figure(figsize=(10, 6))

# Plot original data with smaller marker size
plt.plot(data['x'], data['y'], marker='o', color='black', markersize=5, label='Original Data')  # Reduced marker size

# Plot moving average
plt.plot(data['x'], data['moving_average'], color='grey', label='Moving Average')

# Fill between for moving standard deviation
plt.fill_between(data['x'], data['moving_average'] - data['moving_std'], 
                 data['moving_average'] + data['moving_std'], color='lightgrey', alpha=0.5, label='Moving Std Dev')

# Vertical line for convergence point
if convergence_point:
    plt.axvline(x=convergence_point, color='darkgrey', linestyle='--', label=f'Convergence Point: {convergence_point}')

plt.title("DOC_micromolar vs. X with Moving Average and Std Dev (Grayscale)")
plt.xlabel("X")
plt.ylabel("DOC_micromolar")
plt.legend()
plt.grid(True)
plt.show()

convergence_point

