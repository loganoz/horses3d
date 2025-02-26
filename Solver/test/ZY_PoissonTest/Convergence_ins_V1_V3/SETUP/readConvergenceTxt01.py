import matplotlib.pyplot as plt
import re
import numpy as np


def read_error_data(filename):
    iterations = []
    errors = []

    with open(filename, 'r') as file:
        for line in file:
            # Split each line by spaces and extract iteration and error value
            parts = line.split()
            iteration = int(parts[0])  # First element is iteration number
            error = float(parts[1])    # Second element is error value
            iterations.append(iteration)
            errors.append(error)
    
    return iterations, errors

# Reading data from file
filename = 'testConvergence.txt'
filenameSin = 'testConvergence_sin2PiXsin2PiY_unStr.txt'

fileNameSin2D = 'testConvergence_sin2PiXsin2PiY_Str.txt'
iterations, errors = read_error_data(filename)
iterationsSin, errorsSin = read_error_data(filenameSin)
iterationsSin02, errorsSin02 = read_error_data(fileNameSin2D)

labelList = [
    r"$\partial_x^2 P = 300$",
    r"$\partial_x^2 P = 100/(2.5 \pi)^2 \sin(2.5 \pi x)$",
    r"$\nabla^2 P = \sin(\pi x) \sin(\pi y)$",
]
# labelList = [
#     r"$\nabla^2 = 300$",
#     r"$\nabla^2 = 100/(2.5 \pi)^2 \sin(2.5 \pi x)$",
# ]

print(iterations)
print(errors)

iterations = np.array(iterations)

errExp = 1e-3*np.exp(-np.array(iterations)/2)
# err2 = 1e-4/(np.array(iterations))**(1)
err2 = (1./iterations)**(iterations)
print(errExp)

# Plotting the error data
plt.figure(figsize=(10, 6))
# plt.plot(iterations, errors, marker='o', linestyle='-.', color='b', label=labelList[0])
plt.plot(iterationsSin, errorsSin, marker='^', linestyle='--', color='r', label=labelList[1])
plt.plot(iterationsSin02, errorsSin02, marker='H', linestyle='--', color='orange', label=labelList[2])

plt.plot(iterations, errExp, marker=' ', linestyle='-', color='green', label='covExp',linewidth = 2)
plt.plot(iterations, err2, marker=' ', linestyle='-', color='cyan', label='err1',linewidth = 2)
# plt.rcParams['font.family']


# Setting ticks to point inward and adjusting font size for ticks
plt.tick_params(axis='both', which='both', direction='in', labelsize=18)

# Adding title and labels
plt.title('Error vs p',fontsize =18)
plt.xlabel('k',fontsize =18)
plt.ylabel('Error',fontsize =18)

# Setting logarithmic scale for Y-axis
plt.yscale('log')
# plt.xscale('log')

# Displaying the grid and legend
plt.grid(True)
plt.legend(fontsize = 10)

# Save the plot to a file
plt.savefig('error_plot.png')

# Show the plot
plt.show()
