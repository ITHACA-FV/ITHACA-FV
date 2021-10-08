from files import *
#from withNet import *
from simpleTurbGeomClosed import *

nNut = [25, 25]#, 25, 25, 25, 25, 25, 25, 25]#, 10, 10, 10, 10, 10]
NU = [45, 50]#5, 10, 15, 20, 25, 35, 40, 45, 50]
NP = [45, 50]# 10, 15, 20, 25, 35, 40, 45, 50]

for nu, np, nnut in zip(NU, NP, nNut):
    Netok = Net(nu, nnut, 20000)
    Netok.read()
    Netok.train()
    Netok.save()
