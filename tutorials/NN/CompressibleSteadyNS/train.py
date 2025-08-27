#from withNet import *
from compSteady import *
# from withNet import *
import os
# sed_variable("NmodesUproj", "system/ITHACAdict", 10)
# nNut = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
nNut = [5]#, 25]#, 25, 25, 25, 25, 25, 25, 25]#, 10, 10, 10, 10, 10]
NU = [15]#, 50]#5, 10, 15, 20, 25, 35, 40, 45, 50]
NP = [30]#, 50]# 10, 15, 20, 25, 35, 40, 45, 50]
#nNut = [10, 10, 10, 10, 10]
#NU = [30, 35, 40, 45, 50]
#NP = [30, 35, 40, 45, 50]
# nNut = [10, 10, 10, 10]
# NU = [35, 40, 45, 50]
# NP = [35, 40, 45, 50]


# for nu, np in zip(NU, NP):
#     Netok = Net(nu, nNut, 20000, "relu", 0.0)
#     Netok.read()
#     Netok.trainNet()
#     Netok.save()

for nu, np, nnut in zip(NU, NP, nNut):
    Netok = Net(nu, nnut, 20000)
    Netok.read()
    Netok.train()
    Netok.save()

# for nu, np in zip(NU,NP):
#     sed_variable("NmodesUproj", "system/ITHACAdict", nu)
#     sed_variable("NmodesPproj", "system/ITHACAdict", np)
#     sed_variable("NmodesNutproj", "system/ITHACAdict", nNut)
#     os.system("01simpleTurbGeomClosed")
