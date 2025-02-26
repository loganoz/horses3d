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
# dense_matrix = dense_matrix - dense_matrix.transpose()

# Define x and y coordinates
x = np.arange(num_rows)
y = np.arange(num_rows)
X, Y = np.meshgrid(x, y)

# Plotting using Matplotlib
plt.figure(figsize=(14, 12))
ax = plt.gca()  # Get current axis

contour = plt.pcolormesh(X, Y, dense_matrix, shading='auto', cmap='seismic')
# contour = plt.contourf(X, Y, dense_matrix, cmap='seismic')
# contour = plt.contourf(X, Y, dense_matrix, cmap='Greys')

# plt.contourf(dense_matrix, cmap='viridis')
plt.colorbar(contour)
plt.title('Contour Plot of Dense Matrix')
plt.xlabel('Column Index')
plt.ylabel('Row Index')
ax.set_aspect('equal', 'box')  # Set the aspect of the plot to be equal
dense_matrixMax = np.max(abs(dense_matrix))
plt.clim(-np.abs(dense_matrixMax), np.abs(dense_matrixMax))
# plt.clim(-np.abs(dense_matrix.max()), np.abs(dense_matrix.max()))
# 设置网格线间隔
plt.xticks(np.arange(0, num_rows+2, 9)+0.0)
plt.yticks(np.arange(0, num_rows+2, 9)+0.0)
# plt.xticks(np.arange(-0.5, num_rows+2, 9)+0.0)
# plt.yticks(np.arange(-0.5, num_rows+2, 9)+0.0)
# plt.xticks(np.array([0,8,9,17,18,26,27,35,36,44,45,53]))
# plt.yticks(np.array([0,8,9,17,18,26,27,35,36,44,45,53]))
# plt.xticks(np.array([0,8,9,17,18,26,27,35,36,44,45,53]))
# plt.yticks(np.array([0,8,9,17,18,26,27,35,36,44,45,53]))
# plt.xticks(np.arange(0, num_rows+2, 9))
# plt.yticks(np.arange(0, num_rows+2, 9))
plt.grid(True)

plt.savefig('contour_plot.png')  # 保存图像到文件
print(" plt.savefig('contour_plot.png')  # 保存图像到文件 ")
plt.show()