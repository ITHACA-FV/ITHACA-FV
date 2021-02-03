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

singVal = np.loadtxt("../caseDir/ITHACAoutput/parameterizedBCtest_RBFparameter/singularValues_mat.txt")
relErr_L2norm = np.loadtxt("../caseDir/ITHACAoutput/parameterizedBCtest_RBFparameter/relError_L2norm_mat.txt")
relErr_LinfNorm = np.loadtxt("../caseDir/ITHACAoutput/parameterizedBCtest_RBFparameter/relError_LinfNorm_mat.txt")
condNumber = np.loadtxt("../caseDir/ITHACAoutput/parameterizedBCtest_RBFparameter/condNumber_mat.txt")

singVal = singVal / singVal[0,:][None,:]

xAxis = np.linspace(1, singVal.shape[0], num = singVal.shape[0])

shapePar = np.loadtxt("../caseDir/ITHACAoutput/parameterizedBCtest_RBFparameter/rbfShapeParameters_mat.txt")

#labels=[r'$\eta = 100$', r'$\eta = 10$', r'$\eta = 1$', r'$\eta = 0.3$', r'$\eta = 0.1$', r'$\eta = 0.03$', r'$\eta = 0.01$', r'$\eta = 0.0033$', r'$\eta = 0.001$']
#labels=[r'$\eta = 0.001$', r'$\eta = 0.004$', r'$\eta = 0.2$', r'$\eta = 0.07$', r'$\eta = 0.3$', r'$\eta = 1.33$', r'$\eta = 5.6$', r'$\eta = 23.7$', r'$\eta = 100$']
#
#
#
#
#f = plt.figure(2,figsize=(12,8))
#for i in range(singVal.shape[1]):
#    plt.semilogy(xAxis,singVal[:,i],'o', markersize=15, label=labels[i])
#
#
#
##plt.semilogy(xAxis, singVal, "o", markersize=15)
#plt.legend(loc='best', fontsize=15)
#plt.ylabel("Normalized singular values", fontsize=25)
#plt.xlabel("Singular values index", fontsize=25)
#plt.grid()



fig, ax1 = plt.subplots(figsize=(12,8))

color = 'tab:red'
ax1.set_xlabel(r'RBF shape parameter, $\eta$', fontsize=25)
ax1.set_ylabel('Relative error', fontsize=25, color=color)
ax1.loglog(shapePar, relErr_L2norm, 'o', markersize=10, label = r'$L^2-norm$', color=color)
ax1.loglog(shapePar, relErr_LinfNorm, 'v', markersize=10, label = r'$L^\infty-norm$', color=color)
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_ylim(1e-4, 1e6)

plt.grid()
plt.legend(loc=3, fontsize=25)
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('Condition number',fontsize=25, color=color)  # we already handled the x-label with ax1
ax2.loglog(shapePar, condNumber, 'X', markersize=10, color=color)
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_ylim(1e13, 1e20)

fig.tight_layout()  # otherwise the right y-label is slightly clipped

#plt.semilogy(xAxis, singVal, "o", markersize=15)
#plt.legend(loc='best', fontsize=15)
#plt.ylabel("Normalized singular values", fontsize=25)
#plt.grid()

#plt.loglog(shapePar, relErr_L2norm, 'bo', markersize=15, label = r'$L^2-norm$')
#plt.loglog(shapePar, relErr_LinfNorm, 'kv', markersize=15, label = r'$L^\infty-norm$')
#plt.loglog(shapePar, condNumber, 'gX', markersize=15, label = r'$L^\infty-norm$')
#plt.xlabel(r'RBF shape parameter, $\rho$', fontsize=25)
#plt.ylabel('Relative error', fontsize=25)
#plt.grid()
#plt.legend(loc='best', fontsize=25)

plt.show()
