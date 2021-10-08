#from withNet import *
from simpleTurbGeomClosed import *
from files import *

nNut = read_variable("NmodesNutProj", "system/ITHACAdict")
NU = read_variable("NmodesUproj", "system/ITHACAdict")
NP = read_variable("NmodesPproj", "system/ITHACAdict")

Netok = Net(NU, nNut, 20000)
Netok.read()
Netok.train()
Netok.save()

#nNut = [25, 25]
#NU = [45, 50]
#NP = [45, 50]

#for nu, np, nnut in zip(NU, NP, nNut):
#    Netok = Net(nu, nnut, 20000)
#    Netok.read()
#    Netok.train()
#    Netok.save()
