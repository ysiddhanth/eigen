import numpy as np
import matplotlib.pyplot as plt

# Read the data from the file
data = np.loadtxt('value.txt')
data1 = np.loadtxt('value1.txt')

# Extract n values (array sizes) and corresponding time taken
n_values = data[:, 0]
time_values = data[:, 1]
n1_values = data1[:, 0]
time1_values = data1[:, 1]

# Plot the original data and the fitted curve
plt.figure(figsize=(10, 6))
plt.plot(n_values, time_values, color='blue', label='Aggressive')
plt.plot(n1_values, time1_values, color='red', label='Non Aggressive')
plt.xlabel('Array Size (n)')
plt.ylabel('Time Taken (seconds)')
plt.title('Aggressive vs Non Agressive Deflation in Random Complex Matrices')
plt.legend()
plt.grid(True)
plt.show()

