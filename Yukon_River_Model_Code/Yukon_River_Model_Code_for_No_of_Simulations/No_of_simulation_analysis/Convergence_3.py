import pandas as pd
import matplotlib.pyplot as plt

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
print("convergence_point",convergence_point)

# Plotting the data
plt.figure(figsize=(10, 6))
#plt.plot(data['x'], data['y'], marker='o', color='black', label='Original Data')
plt.plot(data['x'], data['y'], color='black', label='Original Data')
plt.ylim(bottom=0)
#plt.plot(data['x'], data['moving_average'], color='grey', label='Moving Average')
#plt.fill_between(data['x'], data['moving_average'] - data[moving_std'], 
#                 data['moving_average'] + data['moving_std'], color='grey#', alpha=0.3, label='Moving Std Dev')
#if convergence_point:
#    plt.axvline(x=convergence_point, color='black', linestyle='--', label=f'Convergence Point: {convergence_point}')
#plt.title("DOC_micromolar vs. X with Moving Average and Std Dev", color='black')
plt.xlabel("No of simulations", color='black',fontsize=16, fontstyle='italic', fontweight='bold', family='Times New Roman')
plt.ylabel("Variance of [DOC]", color='black',fontsize=16, fontstyle='italic', fontweight='bold', family='Times New Roman')
#plt.xlabel("No of simulations", color='black',fontsize=16, fontweight='bold', family='Times New Roman')
#plt.ylabel("Variance of [DOC]", color='black',fontsize=16, fontweight='bold', family='Times New Roman')
plt.xticks(fontsize=14, fontweight='bold', family='Times New Roman')
plt.yticks(fontsize=14, fontweight='bold', family='Times New Roman')
#plt.legend()
#plt.grid(True)
plt.show()

convergence_point

x = [
    100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
    1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000,
    2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000,
    3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000,
    4500, 5000, 5500, 6000, 7000, 8000, 9000, 10000, 11000, 12000,
    13000, 14000, 15000
]
y = [
    8.520, 8.310, 8.501, 8.311, 8.398, 8.358, 8.141, 8.301, 8.267, 8.446, 8.293, 
    8.311, 8.304, 8.325, 8.289, 8.365, 8.368, 8.384, 8.320, 8.251, 8.378, 8.364, 
    8.409, 8.349, 8.339, 8.299, 8.401, 8.349, 8.366, 8.308, 8.303, 8.302, 8.307, 
    8.301, 8.354, 8.322, 8.350, 8.353, 8.297, 8.294, 8.314, 8.312, 8.358, 8.337, 
    8.267, 8.331, 8.340, 8.310, 8.289, 8.324, 8.284, 8.336, 8.318
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
print("convergence_point",convergence_point)

# Plotting the data
plt.figure(figsize=(10, 6))
#plt.plot(data['x'], data['y'], marker='o', color='black', label='Original Data')
plt.plot(data['x'], data['y'], color='black', label='Original Data')
plt.ylim(bottom=0,top=10)
#plt.plot(data['x'], data['moving_average'], color='grey', label='Moving Average')
#plt.fill_between(data['x'], data['moving_average'] - data['moving_std'], 
#                 data['moving_average'] + data['moving_std'], color='grey', alpha=0.3, label='Moving Std Dev')
if convergence_point:
    plt.axvline(x=convergence_point, color='black', linestyle='--', label=f'Convergence Point: {convergence_point}')
#plt.title("DOC_micromolar vs. X with Moving Average and Std Dev", color='black')
plt.xlabel("No of simulations", color='black')
plt.ylabel("Variance of [Protein]", color='black')
#plt.legend()
#plt.grid(True)
plt.show()

convergence_point

x = [
    100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
    1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000,
    2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000,
    3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000,
    4500, 5000, 5500, 6000, 7000, 8000, 9000, 10000, 11000, 12000,
    13000, 14000, 15000
]
y = [
    302.848, 286.034, 267.893, 272.768, 251.422, 258.161, 272.215, 257.592, 265.264, 
    251.858, 265.368, 274.740, 256.664, 259.226, 268.543, 260.243, 263.720, 258.702, 
    265.822, 274.071, 256.875, 257.353, 259.495, 260.787, 256.807, 262.880, 253.470, 
    267.553, 263.640, 258.875, 264.642, 259.596, 264.427, 261.189, 264.946, 265.158, 
    261.446, 263.462, 265.651, 270.079, 262.009, 266.007, 260.588, 262.830, 262.660, 
    262.136, 260.150, 262.139, 259.855, 262.587, 264.841, 260.637, 264.572
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
print("convergence_point",convergence_point)

# Plotting the data
plt.figure(figsize=(10, 6))
#plt.plot(data['x'], data['y'], marker='o', color='black', label='Original Data')
plt.plot(data['x'], data['y'], color='black', label='Original Data')
plt.ylim(bottom=0,top=350)
#plt.plot(data['x'], data['moving_average'], color='grey', label='Moving Average')
#plt.fill_between(data['x'], data['moving_average'] - data['moving_std'], 
#                 data['moving_average'] + data['moving_std'], color='grey', alpha=0.3, label='Moving Std Dev')
if convergence_point:
    plt.axvline(x=convergence_point, color='black', linestyle='--', label=f'Convergence Point: {convergence_point}')
#plt.title("DOC_micromolar vs. X with Moving Average and Std Dev", color='black')
plt.xlabel("No of simulations", color='black')
plt.ylabel("Variance of [Heteropolycondenstate]", color='black')
#plt.legend()
#plt.grid(True)
plt.show()

convergence_point



x = [
    100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
    1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000,
    2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000,
    3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000,
    4500, 5000, 5500, 6000, 7000, 8000, 9000, 10000, 11000, 12000,
    13000, 14000, 15000
]
y = [
    16.217, 16.206, 16.433, 16.335, 16.329, 16.401, 16.278, 16.342, 16.194, 16.558, 
    16.167, 16.154, 16.275, 16.343, 16.278, 16.329, 16.241, 16.336, 16.231, 16.169, 
    16.314, 16.329, 16.402, 16.366, 16.340, 16.326, 16.446, 16.230, 16.263, 16.365, 
    16.231, 16.291, 16.258, 16.303, 16.326, 16.309, 16.352, 16.291, 16.262, 16.153, 
    16.341, 16.282, 16.293, 16.269, 16.284, 16.253, 16.288, 16.294, 16.319, 16.277, 
    16.227, 16.292, 16.257
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
print("convergence_point",convergence_point)

# Plotting the data
plt.figure(figsize=(10, 6))
#plt.plot(data['x'], data['y'], marker='o', color='black', label='Original Data')
plt.plot(data['x'], data['y'], color='black', label='Original Data')
plt.ylim(bottom=0,top=20)
#plt.plot(data['x'], data['moving_average'], color='grey', label='Moving Average')
#plt.fill_between(data['x'], data['moving_average'] - data['moving_std'], 
#                 data['moving_average'] + data['moving_std'], color='grey', alpha=0.3, label='Moving Std Dev')
if convergence_point:
    plt.axvline(x=convergence_point, color='black', linestyle='--', label=f'Convergence Point: {convergence_point}')
#plt.title("DOC_micromolar vs. X with Moving Average and Std Dev", color='black')
plt.xlabel("No of simulations", color='black')
plt.ylabel("Variance of [Polysaccharide]", color='black')
#plt.legend()
#plt.grid(True)
plt.show()

convergence_point

x = [
    100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
    1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000,
    2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000,
    3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000,
    4500, 5000, 5500, 6000, 7000, 8000, 9000, 10000, 11000, 12000,
    13000, 14000, 15000
]
y = [
    0.368, 0.377, 0.361, 0.374, 0.355, 0.377, 0.367, 0.365, 0.371,
    0.367, 0.367, 0.369, 0.376, 0.369, 0.379, 0.364, 0.381, 0.374,
    0.366, 0.365, 0.372, 0.375, 0.397, 0.359, 0.369, 0.369, 0.375,
    0.370, 0.374, 0.370, 0.374, 0.373, 0.367, 0.368, 0.373, 0.366,
    0.377, 0.376, 0.373, 0.371, 0.380, 0.375, 0.369, 0.368, 0.364,
    0.370, 0.373, 0.371, 0.371, 0.372, 0.371, 0.374, 0.368
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
print("convergence_point",convergence_point)

# Plotting the data
plt.figure(figsize=(10, 6))
#plt.plot(data['x'], data['y'], marker='o', color='black', label='Original Data')
plt.plot(data['x'], data['y'], color='black', label='Original Data')
plt.ylim(bottom=0,top=.5)
#plt.plot(data['x'], data['moving_average'], color='grey', label='Moving Average')
#plt.fill_between(data['x'], data['moving_average'] - data['moving_std'], 
#                 data['moving_average'] + data['moving_std'], color='grey', alpha=0.3, label='Moving Std Dev')
if convergence_point:
    plt.axvline(x=convergence_point, color='black', linestyle='--', label=f'Convergence Point: {convergence_point}')
#plt.title("DOC_micromolar vs. X with Moving Average and Std Dev", color='black')
plt.xlabel("No of simulations", color='black')
plt.ylabel("Variance of [Lipid]", color='black')
#plt.legend()
#plt.grid(True)
plt.show()

convergence_point


x = [
    100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
    1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000,
    2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000,
    3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000,
    4500, 5000, 5500, 6000, 7000, 8000, 9000, 10000, 11000, 12000,
    13000, 14000, 15000
]
y = [
    1.057, 0.966, 0.948, 0.934, 0.738, 0.844, 0.950, 0.814, 0.844,
    0.833, 0.835, 0.882, 0.812, 0.836, 0.863, 0.869, 0.849, 0.841,
    0.831, 0.888, 0.819, 0.836, 0.848, 0.866, 0.871, 0.855, 0.830,
    0.855, 0.831, 0.851, 0.846, 0.846, 0.857, 0.821, 0.856, 0.876,
    0.849, 0.861, 0.852, 0.850, 0.849, 0.862, 0.825, 0.853, 0.848,
    0.846, 0.831, 0.850, 0.848, 0.831, 0.844, 0.838, 0.851
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
print("convergence_point",convergence_point)

# Plotting the data
plt.figure(figsize=(10, 6))
#plt.plot(data['x'], data['y'], marker='o', color='black', label='Original Data')
plt.plot(data['x'], data['y'], color='black', label='Original Data')
plt.ylim(bottom=0)
#plt.plot(data['x'], data['moving_average'], color='grey', label='Moving Average')
#plt.fill_between(data['x'], data['moving_average'] - data['moving_std'], 
#                 data['moving_average'] + data['moving_std'], color='grey', alpha=0.3, label='Moving Std Dev')
if convergence_point:
    plt.axvline(x=convergence_point, color='black', linestyle='--', label=f'Convergence Point: {convergence_point}')
#plt.title("DOC_micromolar vs. X with Moving Average and Std Dev", color='black')
plt.xlabel("No of simulations", color='black')
plt.ylabel("Variance of [CDOM]", color='black')
#plt.legend()
#plt.grid(True)
plt.show()

convergence_point