import numpy as np
import scipy.sparse

npyFile = np.load('forNpy.npz');
np.savez('forNpy.npz',indices=npyFile.f.indices,format=np.array('csc',dtype='|S3'),data=npyFile.f.data,indptr = npyFile.f.indptr, shape = npyFile.f.shape)

a = scipy.sparse.load_npz('forCnpy.npz')
print(a.todense())

