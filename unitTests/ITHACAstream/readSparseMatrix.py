import numpy as np
import scipy.sparse
import os

## Generate a sparse matrix by running StreamTest.C and then run this script to read it back in.
## Sanity check: check that the shape is correct (RowMajor mapping between C++ and Python) and that the data is correct.
## TODO: automate this check
npyFile = np.load("testElement_mat.npy")
print(f"Read sparse matrix from npy file: {npyFile}")
os.remove("testElement_mat.npy")

npyFile = np.load('forNpy.npz')

np.savez('forNpy.npz',indices=npyFile.f.indices,format=np.array('csc',dtype='|S3'),data=npyFile.f.data,indptr = npyFile.f.indptr, shape = npyFile.f.shape)

a = scipy.sparse.load_npz('forCnpy.npz')
print(a.todense())

