import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Create a DataFrame with your data
data = {
    "x": [
        100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
        1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000,
        2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000,
        3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000,
        4500, 5000, 5500, 6000, 7000, 8000, 9000, 10000, 11000, 12000,
        13000, 14000, 15000
    ],
    "Protein_micromolar": [
        8.520, 8.310, 8.501, 8.311, 8.398, 8.358, 8.141, 8.301, 8.267, 8.446,
        8.293, 8.311, 8.304, 8.325, 8.289, 8.365, 8.368, 8.384, 8.320, 8.251,
        8.378, 8.364, 8.409, 8.349, 8.339, 8.299, 8.401, 8.349, 8.366, 8.308,
        8.303, 8.302, 8.307, 8.301, 8.354, 8.322, 8.350, 8.353, 8.297, 8.294,
        8.314, 8.312, 8.358, 8.337, 8.267, 8.331, 8.340, 8.310, 8.289, 8.324,
        8.284, 8.336, 8.318
    ],
    # Add other y-values here in the same format
    "DOC_micromolar": [
        118.050, 107.534, 103.918, 102.719, 82.966, 93.305, 104.037, 90.703, 94.095, 90.908,
        93.377, 99.134, 90.140, 92.473, 96.273, 95.431, 94.526, 92.977, 93.372, 99.617,
        90.930, 92.208, 93.323, 95.463, 95.098, 94.918, 90.865, 95.678, 93.032, 93.729,
        94.363, 93.544, 95.244, 91.585, 95.161, 96.813, 93.913, 95.366, 95.110, 95.735,
        94.242, 95.887, 91.930, 94.604, 94.284, 94.074, 92.336, 94.288, 93.744, 92.756,
        94.274, 92.984, 94.787
    ]
    # Repeat for each different y-value set
}

# Convert the data to a DataFrame
df = pd.DataFrame(data)

# Define the window size for moving average
window_size = 5
threshold = 0.5  # Define a threshold for standard deviation

# Function to calculate moving average and plot
def plot_and_find_convergence(df, column_name):
    df['Moving_Avg'] = df[column_name].rolling(window=window_size).mean()
    df['Std_Dev'] = df[column_name].rolling(window=window_size).std()

    # Determine the convergence point
    convergence_point = np.argmax(df['Std_Dev'].dropna() < threshold) * window_size

    # Plot the data
    plt.figure(figsize=(12, 6))
    plt.plot(df['x'], df[column_name], label=column_name, color='blue')
    plt.plot(df['x'], df['Moving_Avg'], label='Moving Average', color='orange')
    if convergence_point > 0:
        plt.axvline(x=df['x'][convergence_point], color='green', linestyle='--', label='Convergence Point')
    plt.xlabel('X')
    plt.ylabel(column_name)
    plt.legend()
    plt.title(f'{column_name} Convergence Plot')
    plt.show()

    return convergence_point

# Iterate over each column (each y-value set)
for column in df.columns[1:]:  # Skip the first column which is 'x'
    convergence_point = plot_and_find_convergence(df, column)
    print(f'Convergence point for {column}: {convergence_point}')



