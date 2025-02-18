import numpy as np
import os
import files 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd

modesU = [10,15,20]
modesP = [10,15,20]
modesE = [10,15,20]

error_totalU=[]
error_totalP=[]
error_totalE=[]
error_totalNut=[]
errorU=[]
errorP=[]
errorE=[]
errorNut=[]

for k in modesU:
    u = "ITHACAoutput/checkOffSingle/errorU"+str(k)+"_30_mat.py"
    p = "ITHACAoutput/checkOffSingle/errorP"+str(k)+"_30_mat.py"
    e = "ITHACAoutput/checkOffSingle/errorE"+str(k)+"_30_mat.py"
    nut = "ITHACAoutput/checkOffSingle/errorNut"+str(k)+"_30_mat.py"
    m = "errorU"+str(k)+"_30"
    n = "errorP"+str(k)+"_30"
    o = "errorE"+str(k)+"_30"
    q = "errorNut"+str(k)+"_30"
    exec(open(u).read())
    exec("error_totalU.append("+m+")")
    exec(open(p).read())
    exec("error_totalP.append("+n+")")
    exec(open(e).read())
    exec("error_totalE.append("+o+")")
    exec(open(nut).read())
    exec("error_totalNut.append("+q+")")

for j in range(0,len(modesU)):
    errorU.append(np.mean(error_totalU[j]))
for k in range(0,len(modesP)):
    errorP.append(np.mean(error_totalP[k]))
for i in range(0,len(modesE)):
    errorE.append(np.mean(error_totalE[i]))
for l in range(0,len(modesU)):
    errorNut.append(np.mean(error_totalNut[l]))

print(errorU)
print(errorP)
print(errorE)
print(errorNut)

data = {'modes': modesU, 'errorU': errorU, 'errorP': errorP, 'errorE': errorE, 'errorNut': errorNut}
outputData = pd.DataFrame(data)
#outputData = outputData.drop(outputData.columns[[0]], axis=1)
print(outputData) 
outputData.to_csv("errors.dat",sep=' ',index=False)

plt.plot(modesU,errorU,'r', label='Relative error velocity')
plt.plot(modesP,errorP,'b', label='Relative error pressure')
plt.plot(modesE,errorE,'g', label='Relative error energy')
plt.plot(modesU,errorNut,'k', label='Relative error eddy viscosity')
# plt.semilogy(PRO[:,0],PRO[:,1],'k--v', label='Relative error for L2 proj.')
# plt.xlim(5,50)
plt.legend()
# plt.legend(bbox_to_anchor=(.5,  .95), loc=2, borderaxespad=0.) 
plt.xlabel("Number of modes")
plt.ylabel("$L^2$ Rel. Error.")#ndof,err=zip(*l)
plt.grid(True)
#plt.savefig("errors.pdf", bbox_inches='tight')
plt.show()

eigsU = open("ITHACAoutput/POD/CumEigenvalues_U").read()
eigsU = eigsU.splitlines()
eigsU = eigsU[2:52]
eigsP = open("ITHACAoutput/POD/CumEigenvalues_p").read()
eigsP = eigsP.splitlines()
eigsP = eigsP[2:52]
eigsE = open("ITHACAoutput/POD/CumEigenvalues_e").read()
eigsE = eigsE.splitlines()
eigsE = eigsE[2:52]
eigsNut = open("ITHACAoutput/POD/CumEigenvalues_nut").read()
eigsNut = eigsNut.splitlines()
eigsNut = eigsNut[2:52]
eigsU = np.asarray(eigsU)

xAx = np.linspace(1,50,50)

eigsData = {'xAx': xAx, 'eigsU': eigsU, 'eigsP': eigsP, 'eigsE': eigsE, 'eigsNut': eigsNut}
outputEigs = pd.DataFrame(eigsData)
print(outputEigs) 
outputEigs.to_csv("eigs.dat",sep=' ',index=False)


errorU=[]
errorP=[]
errorE=[]
errorNut=[]
parameter=[]

u = "ITHACAoutput/checkOffSingle/errorU20_30_mat.py"
p = "ITHACAoutput/checkOffSingle/errorP20_30_mat.py"
e = "ITHACAoutput/checkOffSingle/errorE20_30_mat.py"
nut = "ITHACAoutput/checkOffSingle/errorNut20_30_mat.py"
exec(open(u).read())
exec(open(p).read())
exec(open(e).read())
exec(open(nut).read())

for i in range(0,len(errorU20_30)):
    errorU.append(errorU20_30[i,0])
    errorP.append(errorP20_30[i,0])
    errorE.append(errorE20_30[i,0])
    errorNut.append(errorNut20_30[i,0])
    parameter.append(i+1)

data = {'parameter': parameter, 'errorU': errorU, 'errorP': errorP, 'errorE': errorE, 'errorNut': errorNut}
outputData = pd.DataFrame(data)
#outputData = outputData.drop(outputData.columns[[0]], axis=1)
print(outputData) 
outputData.to_csv("errors.dat",sep=' ',index=False)

#plt.figure()
#plt.plot(xAx,eigsU,'r', label='Velocity eigs')
#plt.plot(xAx,eigsP,'b', label='Pressure eigs')
#plt.plot(xAx,eigsE,'k', label='Energy eigs')
#plt.plot(xAx,eigsNut,'g', label='Nut eigs')
#plt.legend()
#plt.xlabel("Number of modes")
#plt.ylabel("Cumulative eigenvalues")#ndof,err=zip(*l)
##plt.grid(True)
##plt.savefig("eigs.pdf", bbox_inches='tight')
#plt.show()

plt.figure()
testLoss=np.load('ITHACAoutput/NN/testLoss_20_30.npy')
trainLoss=np.load('ITHACAoutput/NN/trainLoss_20_30.npy')
numLoss=testLoss.size
xAxL=np.linspace(1,numLoss,numLoss)
xAxL=xAxL*100

lossData = {'xAxL': xAxL, 'testLoss': testLoss, 'trainLoss': trainLoss}
outputLoss = pd.DataFrame(lossData)
print(outputLoss) 
outputLoss.to_csv("loss2030.dat",sep=' ',index=False)

plt.semilogy(xAxL,testLoss,'r', label='Loss on test')
plt.semilogy(xAxL,trainLoss,'b', label='Loss on train')
plt.xlabel("Epochs")
plt.ylabel("Loss")
plt.legend()
plt.grid(True)
#plt.savefig("losses.pdf", bbox_inches='tight')
plt.show()
#
#exit()
#
#eigsNutD = [0.85281804336512134768,0.06926044057470572002,0.02879636644942273893,0.01538996635181977962,0.01292705761828407719,
#0.00635956284991536554,0.00398900451129383798,0.00294413962341104037,0.00172340494094143327,0.00113019070430354701,
#0.00090328315589451979,0.00060905291174854700,0.00043208626189313374,0.00032145462155264360,0.00026468247263400244,
#0.00023887245563312948,0.00021517569576473642,0.00016408941134164339,0.00015006431970394682,0.00011378156321130438,
#0.00010347933003014966,0.00010196867607993193,0.00009452973549870097,0.00008371399188836893,0.00007528477599139524,
#0.00005848661149518058,0.00005764043079637341,0.00004815340541869198,0.00004499834307028310,0.00004012199222555402]
#
#xAxN = np.linspace(1,30,30)
#
#plt.figure()
#plt.semilogy(xAxN,eigsNutD,'g', label='Nut eigs')
##plt.legend()
#plt.xlabel("Number of modes")
#plt.ylabel("Eigenvalues for nut")#ndof,err=zip(*l)
#plt.grid(True)
#plt.savefig("eigsNut.pdf", bbox_inches='tight')
#plt.show()






