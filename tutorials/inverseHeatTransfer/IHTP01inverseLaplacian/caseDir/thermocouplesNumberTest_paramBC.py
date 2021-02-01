import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (12, 8),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

relError_L2norm = np.loadtxt("./ITHACAoutput/thermocouplesNumberTest_paramBC/relError_L2norm_mat.txt")
relError_LinfNorm = np.loadtxt("./ITHACAoutput/thermocouplesNumberTest_paramBC/relError_LinfNorm_mat.txt")
relError_L2norm_BF = np.loadtxt("./ITHACAoutput/thermocouplesNumberTest_paramBC_bestFit/relError_L2norm_mat.txt")
relError_LinfNorm_BF = np.loadtxt("./ITHACAoutput/thermocouplesNumberTest_paramBC_bestFit/relError_LinfNorm_mat.txt")
relError_L2norm_int = np.loadtxt("./ITHACAoutput/thermocouplesNumberTest_paramBC_bestInt/relError_L2norm_mat.txt")
relError_LinfNorm_int = np.loadtxt("./ITHACAoutput/thermocouplesNumberTest_paramBC_bestInt/relError_LinfNorm_mat.txt")
TCplane_Y = np.loadtxt("./ITHACAoutput/thermocouplesNumberTest_paramBC/numberTCperAxis_mat.txt")


fig, axes = plt.subplots()

plt.semilogy(TCplane_Y, relError_L2norm, "bo", markersize=15)
plt.semilogy(TCplane_Y, relError_LinfNorm, "bv", markersize=15)
plt.semilogy(TCplane_Y, relError_L2norm_BF, "ko", markersize=15,    label = r'$||\epsilon||_{L^2(\Gamma_{s_{in}})}$')
plt.semilogy(TCplane_Y, relError_LinfNorm_BF, "kv", markersize=15,  label = r'$||\epsilon||_{L^\infty(\Gamma_{s_{in}})}$')
plt.semilogy(TCplane_Y, relError_L2norm_int, "go", markersize=15)
plt.semilogy(TCplane_Y, relError_LinfNorm_int, "gv", markersize=15)
plt.xlabel(r'Number of thermocouples per axis', fontsize=25)
plt.xlim(0, TCplane_Y[-1])
plt.grid()
#plt.title(r"Parameterized BC (LU)", fontsize=25)
#plt.yscale('log')
#ax = plt.figure(5).gca()
#ax.xaxis.set_major_locator(MaxNLocator(integer=True))

leg = plt.legend(fontsize=25)
axes.add_artist(leg)
h = [plt.plot([],[], color=i, marker='o', markersize=15, ls="")[0] for i in ["b", "k", "g"]]# for j in ["-" "--"]]
plt.legend(fontsize=25, handles=h, labels=["Inverse", "BestFit", "Interpolation"], ncol = 3, bbox_to_anchor=(0., 1.1),loc=2, borderaxespad=0.)
#h = [plt.plot([],[], color=i, marker='o', markersize=15, ls="")[0] for i in ["b", "k"]]# for j in ["-" "--"]]
#plt.legend(fontsize=25, handles=h, labels=["Inverse", "Projection"], ncol = 3, bbox_to_anchor=(0., 1.1),loc=2, borderaxespad=0.)


plt.show()
