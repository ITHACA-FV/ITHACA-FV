from of_pybind11_system import of_pybind11_system
from scipy.sparse.linalg import spsolve

#Instantiate OF object
a = of_pybind11_system(["."])
#Get Temperature (T) Field from OF (the memory is shared with OF)
T = a.getT()
#Get Temperature (T) Field from OF (the memory is shared with OF)
S = a.getS()
# Set the first element of the source term to 1
S [0,0] = 1
# solve the problem using OF
a.solve()
# store the temperature field to folder 1
a.exportT(".","1","T")
# Set the first element of the source term to 1
S [0,0] = 4
# solve the problem using OF
a.solve()
# store the temperature field to folder 2
a.exportT(".","2","T")
# Export the system Matrix as a sparse scipy matrix (it is copy)
A = a.get_system_matrix(T,S)
# Export the rhs as a dense numpy vector (it is a copy)
b = a.get_rhs(T,S)
# Export aolve the problem using python
sol = spsolve(A, b)
# Print the solution
print(sol)
# Change the first element
sol[0] = 2
# Pass it to OF
a.setT(sol)
# Export again the solution
a.exportT(".","3","T")


