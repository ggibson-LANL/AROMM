import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Given data
x = [
    100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
    1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000,
    2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000,
    3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000,
    4500, 5000, 5500, 6000, 7000, 8000, 9000, 10000, 11000, 12000,
    13000, 14000, 15000
]
y = [
    118.050, 107.534, 103.918, 102.719, 82.966, 93.305, 104.037, 90.703, 94.095, 90.908,
    93.377, 99.134, 90.140, 92.473, 96.273, 95.431, 94.526, 92.977, 93.372, 99.617,
    90.930, 92.208, 93.323, 95.463, 95.098, 94.918, 90.865, 95.678, 93.032, 93.729,
    94.363, 93.544, 95.244, 91.585, 95.161, 96.813, 93.913, 95.366, 95.110, 95.735,
    94.242, 95.887, 91.930, 94.604, 94.284, 94.074, 92.336, 94.288, 93.744, 92.756,
    94.274, 92.984, 94.787
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
plt.plot(data['x'], data['y'], marker='o', label='Original Data')
plt.plot(data['x'], data['moving_average'], color='red', label='Moving Average')
plt.fill_between(data['x'], data['moving_average'] - data['moving_std'], 
                 data['moving_average'] + data['moving_std'], color='red', alpha=0.2, label='Moving Std Dev')
plt.axvline(x=convergence_point, color='green', linestyle='--', label=f'Convergence Point: {convergence_point}')
plt.title("DOC_micromolar vs. X with Moving Average and Std Dev")
plt.xlabel("X")
plt.ylabel("DOC_micromolar")
plt.legend()
plt.grid(True)
plt.show()

convergence_point

