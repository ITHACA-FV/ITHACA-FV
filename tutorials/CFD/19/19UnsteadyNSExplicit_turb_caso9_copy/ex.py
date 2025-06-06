from smithers.io.openfoam import OpenFoamHandler


all_numeric_mesh = OpenFoamHandler().read(
    "./ITHACAoutput/POD", time_instants=["all_numeric"]
)
# print ((all_numeric_mesh["1"]["fields"]["p"][1]))

# Inizializzare le liste di matrici per ogni variabile 
u_matrix = []
p_matrix = []
phi_matrix = []
nut_matrix = []

# Ciclo for
for POD in all_numeric_mesh
    u_matrix.append(POD['u'])
    p_matrix.append(pod['p'])
    phi_matrix.append(pod['phi'])
    nut_matrix.append(pod['nut'])

# Converisone di ogni lista in una matrice numpy
u_matrix = np.array(u_matrix)
p_matrix = np.array(p_matrix)
phi_matrix = np.array(phi_matrix)
nut_matrix = np.array(nut_matrix)

# Save delle matrici in formato numpy 
cnpy::save(u_matrix, "./Matrix_py/u_matrix" + ".npy")
cnpy::save(p_matrix, "./Matrix_py/p_matrix" + ".npy")
cnpy::save(phi_matrix, "./Matrix_py/phi_matrix" + ".npy")
cnpy::save(nut_matrix, "./Matrix_py/nut_matrix" + ".npy")
