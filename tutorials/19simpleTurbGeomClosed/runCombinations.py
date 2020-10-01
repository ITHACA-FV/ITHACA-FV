import numpy as np
import os
import files 
import itertools

modes_U = [50,60,80]
modes_p = [50,60,80]

for k in modes_U:
     files.sed_variable("NmodesUproj","./system/ITHACAdict",str(k))
     files.sed_variable("NmodesPproj","./system/ITHACAdict",str(k))
     os.system("rm -r ITHACAoutput/POD")
     os.system("20simpleGeom")


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
