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




fomforcex1 = np.load("ITHACAoutput/DataFromFoam_1/fomforcex.npy")
fomforcex2 = np.load("ITHACAoutput/DataFromFoam_2/fomforcex.npy")
fomforcex3 = np.load("ITHACAoutput/DataFromFoam_3/fomforcex.npy")
fomforcex4 = np.load("ITHACAoutput/DataFromFoam_4/fomforcex.npy")
fomforcex5 = np.load("ITHACAoutput/DataFromFoam_5/fomforcex.npy")


plt.xlabel(r'\textit{Time} ~[s]',  fontsize=20)
plt.ylabel(r'$\textit{Drag} ~[N]$',fontsize=20)

plt.plot(fomforcex1, color="blue", linestyle='-', linewidth=2, label=r"$param_{1}$")
plt.plot(fomforcex2, color="red", linestyle='-', linewidth=2, label=r"$param_{2}$")
plt.plot(fomforcex3, color="orange", linestyle='-', linewidth=2, label=r"$param_{3}$")
plt.plot(fomforcex4, color="green", linestyle='-', linewidth=2, label=r"$param_{4}$")
plt.plot(fomforcex5, color="k", linestyle='-', linewidth=2, label=r"$param_{5}$")


plt.legend(ncol=5, loc="best",fontsize=14)
plt.savefig("Drag.pdf")
#plt.show()
