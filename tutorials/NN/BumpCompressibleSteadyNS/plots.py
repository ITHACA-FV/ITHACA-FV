import numpy as np
import os
import files 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


errorU=[]
errorP=[]
errorE=[]

u = "./ITHACAoutput/checkOffSingle/errorU_mat.py"
p = "./ITHACAoutput/checkOffSingle/errorP_mat.py"
e = "./ITHACAoutput/checkOffSingle/errorE_mat.py"
n = "./ITHACAoutput/checkOffSingle/errorNut_mat.py"

xAx = np.linspace(1,20,20)

exec(open(u).read())

plt.semilogy(xAx,errorU,'o', label='Relative error for U')
# plt.xlim(5,50)
plt.xlabel("Solution number")
plt.ylabel("L2 Relative Error for Velocity")
# plt.legend(bbox_to_anchor=(.5,  .95), loc=2, borderaxespad=0.) 
plt.grid(True)
# f.savefig("poisson.pdf", bbox_inches='tight')
plt.show()

exec(open(p).read())

plt.semilogy(xAx,errorP,'^', label='Relative error for P')
# plt.xlim(5,50)
plt.xlabel("Solution number")
plt.ylabel("L2 Relative Error for Pressure")
# plt.legend(bbox_to_anchor=(.5,  .95), loc=2, borderaxespad=0.) 
plt.grid(True)
# f.savefig("poisson.pdf", bbox_inches='tight')
plt.show()

exec(open(e).read())

plt.semilogy(xAx,errorE,'s', label='Relative error for E')
# plt.xlim(5,50)
plt.xlabel("Solution number")
plt.ylabel("L2 Relative Error for Energy")
# plt.legend(bbox_to_anchor=(.5,  .95), loc=2, borderaxespad=0.) 
plt.grid(True)
# f.savefig("poisson.pdf", bbox_inches='tight')
plt.show()

exec(open(n).read())

plt.semilogy(xAx,errorNut,'d', label='Relative error for Nut')
# plt.xlim(5,50)
plt.xlabel("Solution number")
plt.ylabel("L2 Relative Error for Eddy viscosity")
# plt.legend(bbox_to_anchor=(.5,  .95), loc=2, borderaxespad=0.) 
plt.grid(True)
# f.savefig("poisson.pdf", bbox_inches='tight')
plt.show()

plt.figure()
plt.semilogy(xAx,errorU,'o', label='Relative error for U')
plt.semilogy(xAx,errorP,'^', label='Relative error for P')
plt.semilogy(xAx,errorE,'s', label='Relative error for E')
plt.semilogy(xAx,errorNut,'d', label='Relative error for Nut')
plt.legend()
plt.xlabel("Solution number")
plt.ylabel("L2 Relative Errors")
plt.grid(True)
plt.savefig("./ITHACAoutput/Errors.pdf", bbox_inches='tight')
plt.show()




#modes_T = [5,10,15,20,25,30,40]
#modes_DEIM = [40,40,40,40,40,40,40]

#for k,j in zip(modes_T, modes_DEIM):
#     files.sed_variable("N_modes_T","./system/ITHACAdict",k)
#     files.sed_variable("N_modes_DEIM_A","./system/ITHACAdict",j)
#     files.sed_variable("N_modes_DEIM_B","./system/ITHACAdict",j)
#     os.system("10thermalBlock_RBF")
#for k,j in zip(modes_T, modes_DEIM):
#    s = "error_"+str(k)+"_"+str(j)+"_"+str(j)+"_mat.py"
#    m = "error_"+str(k)+"_"+str(j)+"_"+str(j)
#    exec(open(s).read())
#    exec("error_total.append("+m+")")
#for j in range(0,len(modes_T)):
#    error.append(np.mean(error_total[j]))
