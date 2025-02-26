import numpy as np
import matplotlib.pyplot as plt
import matplotlib
# matplotlib.use('TkAgg')  # 或者 'Qt5Agg' 或 'GTK3Agg'

def read_csr_matrix(filename):
    row_ptr = []
    col_indices = []
    values = []
    num_rows = 0
    num_nonzeros = 0

    with open(filename, 'r') as file:
        # Read first line for matrix dimensions
        first_line = file.readline().strip()
        num_rows, num_nonzeros = map(int, first_line.split())

        # Read row pointer
        for _ in range(num_rows + 1):
            row_ptr.append(int(file.readline().strip()))

        # Read column indices
        for _ in range(num_nonzeros):
            col_indices.append(int(file.readline().strip()))

        # Read values
        for _ in range(num_nonzeros):
            values.append(float(file.readline().strip()))  # Assuming values are floats

    return num_rows, row_ptr, col_indices, values

def csr_to_dense(num_rows, row_ptr, col_indices, values):
    # Initialize dense matrix with zeros
    dense_matrix = [[0.0] * num_rows for _ in range(num_rows)]

    # Fill dense matrix with values from CSR format
    for i in range(num_rows):
        for j in range(row_ptr[i]-1, row_ptr[i+1]-1):
            col_idx = col_indices[j]-1
            # print(col_idx, i,j, col_idx)
            dense_matrix[i][col_idx] = values[j]

    return dense_matrix

# Example usage
filename = 'Jacobian.txt'
num_rows, row_ptr, col_indices, values = read_csr_matrix(filename)

# print(num_rows,row_ptr,col_indices,values)
dense_matrix = csr_to_dense(num_rows, row_ptr, col_indices, values)

# # Print the dense matrix
# print("Dense Matrix:")
# for row in dense_matrix:
#     print(row)



dense_matrix = np.array(dense_matrix)

# Define x and y coordinates
x = np.arange(num_rows)
y = np.arange(num_rows)
X, Y = np.meshgrid(x, y)



# Create a figure and two subplots (1 row, 2 columns)
fig, axs = plt.subplots(1, 2, figsize=(24, 13))

contour1 = axs[0].pcolormesh(X, Y, dense_matrix, shading='auto', cmap='seismic')
axs[0].set_title('original')
axs[0].set_xlabel('x')
axs[0].set_ylabel('y')
ax = plt.gca()  # Get current axis
ax.set_aspect('equal', 'box')  # Set the aspect of the plot to be equal

dense_matrixMax = np.max(abs(dense_matrix))
contour1.set_clim(-np.abs(dense_matrixMax), np.abs(dense_matrixMax))

axs[0].legend()
axs[0].set_aspect('equal', 'box')  # Set the aspect of the plot to be equal


axs[0].set_xticks(np.arange(0, num_rows+2, 9)+0.0)
axs[0].set_yticks(np.arange(0, num_rows+2, 9)+0.0)
fig.colorbar(contour1, ax=axs[0])
axs[0].grid(True)

# ==============================================================================
dense_matrix = dense_matrix - dense_matrix.transpose()
dense_matrixMax = np.max(abs(dense_matrix))
# Second plot (cos)
contour2 = axs[1].pcolormesh(X, Y, dense_matrix, shading='auto', cmap='seismic')
axs[1].set_title('minus ---')
axs[1].set_xlabel('x')
axs[1].set_ylabel('y')
axs[1].legend()
axs[1].set_xticks(np.arange(0, num_rows+2, 9)+0.0)
axs[1].set_yticks(np.arange(0, num_rows+2, 9)+0.0)
contour2.set_clim(-np.abs(dense_matrixMax), np.abs(dense_matrixMax))

axs[1].grid(True)
axs[1].set_aspect('equal', 'box')  # Set the aspect of the plot to be equal
fig.colorbar(contour2, ax=axs[1])

plt.savefig('contour_plot.png')  # 保存图像到文件
print(" plt.savefig('contour_plot.png')  # 保存图像到文件 ")
plt.show()