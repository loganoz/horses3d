import pandas as pd
import re
import matplotlib.pyplot as plt


# Specify the file path
file_path = 'convergenceTestOmegaPower2_order3.txt'


def readData001(file_path):
    # Read and parse the file
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Use regular expressions to extract all numbers, then keep only the last five from each line
    data = []
    for line in lines:
        # Extract all numbers, including those in scientific notation
        values = re.findall(r"[-+]?\d*\.\d+(?:[eE][-+]?\d+)?|\d+", line)
        if values:
            # Convert extracted values to floats
            float_values = [float(value) for value in values]
            # Keep only the last five values
            data.append(float_values[-5:])

    # Convert the list of last five values from each line to a Pandas DataFrame
    df = pd.DataFrame(data, columns=['col1', 'col2', 'col3', 'col4', 'col5'])

    return df

 
# Output the result
print("Data for the last five columns of each line:")
# print(df)
file_pathPre = 'convergenceTestOmegaPower2_order'



 # Plot 'mesh' against 'col4'
plt.figure(figsize=(10, 8))


# Define markers for different plots
markers = ['o', 's', 'D', '^', 'v', '>', '<', 'p', '*', 'H']
# Show plot
plt.show()
# Process each file from order 1 to 10
for i in range(1, 11):
    file_path = f"{file_pathPre}{i}.txt"
    try:
        df = readData001(file_path)
        
        # Calculate the division of the second column by the cube of the third column
        df['result'] = df['col2'] / ((df['col3']+1) ** 2)

        df['mesh'] = df['result']**(2/3)

        df['dofsqrt'] = df['mesh']**0.5*(df['col3']+1)

        # Filter out rows where col4 is less than 1e-20
        df = df[df['col4'] >= 1e-20]
        
        print(f"Data with calculated results for file {file_path}:")
        print(df[['col2', 'col3', 'result', 'mesh']])  # Print only the relevant columns for clarity
        # Use the marker corresponding to the current index
        marker0 = markers[i - 1] if i - 1 < len(markers) else markers[-1]
        plt.plot(df['dofsqrt'], df['col4'], 
                 marker=marker0,
                 linestyle='-.', 
                 label=f'p= {i}',
                 markerfacecolor='none',  # Set face color to none for hollow effect
                 )
        
        res = df['dofsqrt'] ** (-df['col3'] - 1)
        
        if(i <=2):
            # Plot the additional function with a reduced line width
            plt.plot(df['dofsqrt'], res/res[2]*df['col4'][2],
                 linestyle='--', 
                 color=plt.gca().lines[-1].get_color(),  # Get the same color as the last plot
                 linewidth=0.8,  # Reduce the line width by a factor of three
                 alpha=1,  # Optional: set transparency for better visibility
                 label=f'slope= -({i}+1)',
                 )
        if(i>2):
            # Plot the additional function with a reduced line width
            plt.plot(df['dofsqrt'], res/res[0]*df['col4'][0],
                 linestyle='--', 
                 color=plt.gca().lines[-1].get_color(),  # Get the same color as the last plot
                 linewidth=0.8,  # Reduce the line width by a factor of three
                 alpha=1,  # Optional: set transparency for better visibility
                 label=f'slope= -({i}+1)',
                 )
        # plt.xlabel('col4')
        # plt.ylabel('mesh')

        # # Plot vertical lines for unique col3 values
        # for col3_value in df['col3'].unique():
        #     plt.axvline(x=col3_value ** 0.5, color='gray', linestyle='--', alpha=0.5)  # Draw vertical lines at col3 values
        
    except FileNotFoundError:
        print(f"File {file_path} not found.")
    except Exception as e:
        print(f"An error occurred while processing {file_path}: {e}")


# Limit the y-axis to the specified range
plt.ylim(1e-8, 7e-1)
plt.xscale('log')  # Set x-axis to logarithmic scale
plt.yscale('log')  # Set y-axis to logarithmic scale
plt.xlabel(r'$\sqrt{dof}$')  # Using LaTeX for the x-label
plt.ylabel('absolute error')
plt.title(r'$-\nabla^2 = sin(2 \pi x) sin(2 \pi z)$')
plt.legend()
plt.grid(True)
plt.savefig('convergence001.png')