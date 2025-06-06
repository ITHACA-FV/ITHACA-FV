from smithers.io.openfoam import OpenFoamHandler
import numpy as np

FOM_data = OpenFoamHandler().read("./ITHACAoutput/Offline/", time_instants="all_numeric")
ROM_data = OpenFoamHandler().read("./ITHACAoutput/POD/", time_instants="all_numeric")

# Lista dei time step come stringhe
time_steps_ROM = [str(i) for i in range(1, 11)]
time_steps_FOM = [str(i) for i in range(1, 2001)]

# Inizializzare le liste di matrici per ogni variabile 
u_matrix_FOM = []
p_matrix_FOM = []
# phi_matrix_FOM = []
nut_matrix_FOM = []

u_matrix_ROM = []
p_matrix_ROM = []
# phi_matrix_ROM = []
nut_matrix_ROM = []

for t in time_steps_FOM:
    fields = FOM_data[t]["fields"]

    u_data = np.array(fields["U"][1]) 
    u_flat = u_data.flatten() 
    print(f"Time {t}: U shape = {u_flat.shape}")  # DEBUG
    u_matrix_FOM.append(u_flat)    
    print(f"Appended to u_matrix_FOM (time {t}) - shape: {u_matrix_FOM[-1].shape}")


    p_data = np.array(fields["p"][1]) 
    p_flat = p_data.flatten() 
    print(f"Time {t}: p shape = {p_flat.shape}")  # DEBUG
    p_matrix_FOM.append(p_flat)
    print(f"Appended to p_matrix_FOM (time {t}) - shape: {p_matrix_FOM[-1].shape}")

    
    # phi_data = np.array(fields["phi"][1]) 
    # phi_flat = phi_data.flatten() 
    # print(f"Time {t}: phi shape = {phi_flat.shape}")  # DEBUG   
    # phi_matrix_FOM.append(phi_flat)
    # print(f"Appended to phi_matrix_FOM (time {t}) - shape: {phi_matrix_FOM[-1].shape}")


    nut_data = np.array(fields["nut"][1]) 
    nut_flat = nut_data.flatten()
    print(f"Time {t}: nut shape = {nut_flat.shape}")  # DEBUG 
    nut_matrix_FOM.append(nut_flat)
    print(f"Appended to nut_matrix_FOM (time {t}) - shape: {nut_matrix_FOM[-1].shape}")


for t in time_steps_ROM:
    fields = ROM_data[t]["fields"]

    u_data = np.array(fields["U"][1]) 
    u_flat = u_data.flatten() 
    print(f"Time {t}: U shape ROM = {u_flat.shape}")  # DEBUG
    u_matrix_ROM.append(u_flat)    
    print(f"Appended to u_matrix_ROM (time {t}) - shape: {u_matrix_ROM[-1].shape}")


    p_data = np.array(fields["p"][1]) 
    p_flat = p_data.flatten() 
    print(f"Time {t}: p shape ROM = {p_flat.shape}")  # DEBUG
    p_matrix_ROM.append(p_flat)
    print(f"Appended to p_matrix_ROM (time {t}) - shape: {p_matrix_ROM[-1].shape}")

    
    # phi_data = np.array(fields["phi"][1]) 
    # phi_flat = phi_data.flatten()
    # print(f"Time {t}: phi shape ROM = {phi_flat.shape}")  # DEBUG    
    # phi_matrix_ROM.append(phi_flat)
    # print(f"Appended to phi_matrix_ROM (time {t}) - shape: {phi_matrix_ROM[-1].shape}")


    nut_data = np.array(fields["nut"][1]) 
    nut_flat = nut_data.flatten()
    print(f"Time {t}: nut shape ROM = {nut_flat.shape}")  # DEBUG  
    nut_matrix_ROM.append(nut_flat)
    print(f"Appended to nut_matrix_ROM (time {t}) - shape: {nut_matrix_ROM[-1].shape}")


lengths = [arr.shape for arr in u_matrix_FOM]
unique_lengths = set(lengths)
print(f"Unique shapes in u_matrix_FOM: {unique_lengths}")

# Converisone di ogni lista in una matrice numpy
u_matrix_FOM = np.array(u_matrix_FOM)
p_matrix_FOM = np.array(p_matrix_FOM)
# phi_matrix_FOM = np.array(phi_matrix_FOM)
nut_matrix_FOM = np.array(nut_matrix_FOM)

u_matrix_ROM = np.array(u_matrix_ROM)
p_matrix_ROM = np.array(p_matrix_ROM)
# phi_matrix_ROM = np.array(phi_matrix_ROM)
nut_matrix_ROM = np.array(nut_matrix_ROM)

# Matrici dei coefficienti
coeff_matrix_u = u_matrix_FOM @ u_matrix_ROM
coeff_matrix_p = p_matrix_FOM @ p_matrix_ROM
# coeff_matrix_phi = phi_matrix_FOM @ phi_matrix_ROM
coeff_matrix_nut = nut_matrix_FOM @ nut_matrix_ROM

# Save delle matrici in formato numpy 
np.save("./Matrix_py/coeff_matrix_u" + ".npy", coeff_matrix_u)
np.save("./Matrix_py/coeff_matrix_p" + ".npy", coeff_matrix_p)
# np.save("./Matrix_py/coeff_matrix_phi" + ".npy", coeff_matrix_phi)
np.save("./Matrix_py/coeff_matrix_nut" + ".npy", coeff_matrix_nut)
