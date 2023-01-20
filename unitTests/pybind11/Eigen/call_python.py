from example import *
import numpy as np

# Create a Random Eigen Matrix
a = Matrix(2,2)

# Print the Matrix
print(a.view_matrix())

# Set the index 0,0 to 1
a.get_matrix()[0,0] = 1

# Print the Matrix again
print(a.view_matrix())

