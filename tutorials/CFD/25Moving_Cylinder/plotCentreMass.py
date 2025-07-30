import matplotlib.pyplot as plt
import numpy as np


plt.rcParams.update({
    "text.usetex": True,
    "font.family": "fantasy",
     "font.size": 24, 
    "font.fantasy": 'Times New Roman',
    'figure.figsize' : (15,6),
    'figure.dpi' : 100
})
plt.style.use('seaborn-whitegrid')



yDispl1 = np.load("ITHACAoutput/DataFromFoam_1/CentreOfMassY.npy")
yDispl2 = np.load("ITHACAoutput/DataFromFoam_2/CentreOfMassY.npy")
yDispl3 = np.load("ITHACAoutput/DataFromFoam_3/CentreOfMassY.npy")
yDispl4 = np.load("ITHACAoutput/DataFromFoam_4/CentreOfMassY.npy")
yDispl5 = np.load("ITHACAoutput/DataFromFoam_5/CentreOfMassY.npy")


plt.xlabel(r'\textit{Time}~ [s]',  fontsize=20)
plt.ylabel(r'$\textit{CentreOfMass}~ [m]$',fontsize=20)

plt.plot(yDispl1, color="blue", linestyle='-', linewidth=2,   label=r"$param_{1}$")
plt.plot(yDispl2, color="red", linestyle='-', linewidth=2, label=r"$param_{2}$")
plt.plot(yDispl3, color="orange", linestyle='-', linewidth=2, label=r"$param_{3}$")
plt.plot(yDispl4, color="green", linestyle='-', linewidth=2, label=r"$param_{4}$")
plt.plot(yDispl5, color="k", linestyle='-', linewidth=2, label=r"$param_{5}$")

plt.legend(ncol=5, loc="best",fontsize=14)
plt.savefig("yDispl.pdf")
#plt.show()
