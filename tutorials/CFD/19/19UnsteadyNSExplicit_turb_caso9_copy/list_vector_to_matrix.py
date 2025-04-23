import numpy as np
import os
import matplotlib.pyplot as plt

########  CREAZIONE MATRICE ########

# Dimensione matrice 
n_files = 2001
n_rows = 31

# Inizializza la matrice finale
matrix = np.zeros((n_rows, n_files))

# Ciclo per leggere i file
for i in range(n_files):
    filename = f"./reduced/sol_reduced_{i}.npy"
    col_data = np.load(filename).reshape(-1)  # Assicura che sia un array 1D
    matrix[:, i] = col_data  # Inserisce come colonna

print("Shape della matrice:", matrix.shape)  

np.save("matrice_completa.npy", matrix)


########   PLOT   #######

# Carica matrice completa
matrice = np.load("matrice_completa.npy")

# Estrai la riga
riga = matrice[1, :]

tempo = np.arange(len(riga))

# Plot della riga rispetto al tempo
plt.figure(figsize=(8, 5))
plt.plot(tempo, riga, linestyle='-', color='b')         # marker='o'     label="Prima riga"
plt.xlabel("Tempo")
# plt.ylabel("")
# plt.title("")
plt.grid(True)
#plt.legend()
plt.show()
plt.savefig("plot.png")


