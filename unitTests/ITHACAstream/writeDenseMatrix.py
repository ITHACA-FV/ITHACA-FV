import numpy as np

a = np.random.rand(5,5).astype('f8')
print(a)
np.save("forNpy.npy",a)



