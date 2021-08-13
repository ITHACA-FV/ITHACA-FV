from of_pybind11 import *
import os
#Instantiate OF object
a = of_pybind11(["."])
#Get Temperature (T) Field from OF
c = a.getT()
#Change it with python
c[0,0] = 1
#Print The changed numpy array
print(a.getT())
#Check the OF field that is changed as well
a.printT()
