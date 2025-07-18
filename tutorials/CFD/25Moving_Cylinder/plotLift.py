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



fomforcey1 = np.load("ITHACAoutput/DataFromFoam_1/fomforcey.npy")
fomforcey2 = np.load("ITHACAoutput/DataFromFoam_2/fomforcey.npy")
fomforcey3 = np.load("ITHACAoutput/DataFromFoam_3/fomforcey.npy")
fomforcey4 = np.load("ITHACAoutput/DataFromFoam_4/fomforcey.npy")
fomforcey5 = np.load("ITHACAoutput/DataFromFoam_5/fomforcey.npy")


plt.xlabel(r'\textit{Time}~ [s]',  fontsize=20)
plt.ylabel(r'$\textit{Lift}~ [N]$',fontsize=20)

# Set up colors using a colormap
colors = plt.cm.tab10(np.linspace(0, 1, 4))  # Get 4 distinct colors from tab10 colormap

plt.plot(fomforcey1, color="blue", linestyle='-', linewidth=2, label=r"$param_{1}$")
plt.plot(fomforcey2, color="red", linestyle='-', linewidth=2, label=r"$param_{2}$")
plt.plot(fomforcey3, color="orange", linestyle='-', linewidth=2, label=r"$param_{3}$")
plt.plot(fomforcey4, color="green", linestyle='-', linewidth=2, label=r"$param_{4}$")
plt.plot(fomforcey5, color="k", linestyle='-', linewidth=2, label=r"$param_{5}$")
#plt.plot(romforcey, 'b--',  label='$ROM$')
#plt.ylim(-0.5, 0.5 ) 

plt.legend(ncol=5, loc="best", fontsize=14)
plt.savefig("Lift.pdf")
#plt.show()
