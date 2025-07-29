import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV file into a DataFrame
df1 = pd.read_csv('17_DOC.csv')
df2 = pd.read_csv('distance.csv')

# Assuming your CSV file has columns named 'x' and 'y'
x = df1['1']
y = df2['1']

# Plot the data
plt.plot(x, y)
plt.xlabel('X-axis Label')
plt.ylabel('Y-axis Label')
plt.title('Title of the Plot')
plt.grid(True)
plt.show()

