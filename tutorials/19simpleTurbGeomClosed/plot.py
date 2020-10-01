import numpy as np
import os
import files 
import itertools
import matplotlib.pyplot as plt

modes_U = [1,2,3,4,5,6,7,8,9,10,15,20,30,40]
modes_p = [1,2,3,4,5,6,7,8,9,10,15,20,30,40]

errorU_G = []
errorP_G = []
errorU_PG = []
errorP_PG = []


for k in modes_U:
#      files.sed_variable("NmodesUproj","./system/ITHACAdict",str(k))
#      files.sed_variable("NmodesPproj","./system/ITHACAdict",str(k))
       varnameU = "errorU_"+str(k)+"_"+str(k) 
       varnameP = "errorP_"+str(k)+"_"+str(k)
       ## Galerkin
       stringfileu = "./G/errorU_"+str(k)+"_"+str(k)+"_mat.py"
       stringfilep = "./G/errorP_"+str(k)+"_"+str(k)+"_mat.py"
       exec(open(stringfileu).read())
       exec(open(stringfilep).read())
       exec("errorU_G.append(np.mean("+varnameU+"))")
       exec("errorP_G.append(np.mean("+varnameP+"))")
       ## PG
       stringfileu = "./PG/errorU_"+str(k)+"_"+str(k)+"_mat.py"
       stringfilep = "./PG/errorP_"+str(k)+"_"+str(k)+"_mat.py"
       exec(open(stringfileu).read())
       exec(open(stringfilep).read())
       exec("errorU_PG.append(np.mean("+varnameU+"))")
       exec("errorP_PG.append(np.mean("+varnameP+"))")
       #exec(open('hello.py').read())
#      os.system("rm -r ITHACAoutput/POD")
#      os.system("20simpleGeom")


plt.semilogy(modes_U,errorU_G, label='U G')
plt.semilogy(modes_U,errorU_PG, label='U PG')

plt.semilogy(modes_U,errorP_G, label='P G')
plt.semilogy(modes_U,errorP_PG, label='P PG')
plt.legend()
plt.show()

# error_total=[]
# error=[]

# for k,j in zip(modes_T, modes_DEIM):
#      s = "error_"+str(k)+"_"+str(j)+"_"+str(j)+"_mat.py"
#      m = "error_"+str(k)+"_"+str(j)+"_"+str(j)
#      exec(open(s).read())
#      exec("error_total.append("+m+")")

# for j in range(0,len(modes_DEIM)):
#     error.append(np.mean(error_total[j]))

# print(error)

# plt.semilogy(modes_DEIM,error,':o', label='Relative error for ROM')
# # plt.semilogy(PRO[:,0],PRO[:,1],'k--v', label='Relative error for L2 proj.')
# # plt.xlim(5,50)
# plt.xlabel("$N$ of modes")
# plt.ylabel("L2 Rel. Error.")

# # plt.legend(bbox_to_anchor=(.5,  .95), loc=2, borderaxespad=0.) 
# plt.grid(True)
# # f.savefig("poisson.pdf", bbox_inches='tight')
# plt.show()
