import numpy as np
import scipy.sparse

sparse_matrix = scipy.sparse.csc_matrix(np.array([[1.1, 2, 0, 3, 10, 10], [4, 0, 0, 6, 5, 0], [4, 0, 0, 6, 0, 0], [4, 0, 0, 6, 0, 0]]))

print(sparse_matrix.todense())

scipy.sparse.save_npz('forCnpy.npz', sparse_matrix, compressed=False)

