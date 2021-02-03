import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (10, 8),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

relError_L2norm = np.loadtxt("../caseDir/ITHACAoutput/thermocouplesLocationTest_paramBC/relError_L2norm_mat.txt")
relError_LinfNorm = np.loadtxt("../caseDir/ITHACAoutput/thermocouplesLocationTest_paramBC/relError_LinfNorm_mat.txt")
TCplane_Y = np.loadtxt("../caseDir/ITHACAoutput/thermocouplesLocationTest_paramBC/TCplane_Y_mat.txt")


f = plt.figure(5,figsize=(12,8))
plt.semilogy(TCplane_Y, relError_L2norm, "bo", markersize=15,label = r'$||\epsilon||_{L^2(\Gamma_{s_{in}})}$')
plt.semilogy(TCplane_Y, relError_LinfNorm, "kv", markersize=15,label = r'$||\epsilon||_{L^\infty(\Gamma_{s_{in}})}$')
plt.xlabel(r'Thermocouples distance from $\Gamma_{s_{in}}$ [m]', fontsize=25)
plt.ylim(1e-3,1e-1)
plt.grid()
plt.title(r"Parameterized BC (LU)", fontsize=25)
plt.legend(fontsize=25)
plt.yscale('log')

plt.show()
