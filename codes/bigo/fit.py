import numpy as np
import matplotlib.pyplot as plt

# Read the data from the file
data = np.loadtxt('value.txt')

# Extract n values (array sizes) and corresponding time taken
n_values = data[:, 0]
time_values = data[:, 1]

# Plot the original data and the fitted curve
plt.figure(figsize=(10, 6))
plt.scatter(n_values, time_values, color='blue', label='Observed Data')
plt.xlabel('Array Size (n)')
plt.ylabel('Time Taken (seconds)')
plt.title('Rea Rand Sym')
plt.legend()
plt.grid(True)
plt.show()

