import numpy as np
import matplotlib.pyplot as plt

# Carica il vettore da error.npy
error = np.load('error_20_20.npy')

# Crea un vettore che va da 1 a 10 (assumendo che il numero di punti sia uguale a quello del vettore error)
x = np.linspace(1, 10, len(error))

# Crea il grafico
plt.plot(x, error, label="Error")

# Aggiungi etichette e titolo
plt.xlabel("Time Step")
plt.ylabel("Error values")
plt.title("Error between U_FOM and U_ROM")

# Griglia sullo sfondo
# plt.grid()

# Abilita i minor ticks (i "piccoli ticks")
plt.minorticks_on()

# Aggiungi la griglia maggiore e minore
plt.grid(True, which='both', linestyle='-', color='gray', linewidth=0.5)

# Puoi personalizzare la densità dei "minor ticks" per avere rettangoli più piccoli
plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(1))  # Modifica la distanza tra i minor ticks sull'asse X
plt.gca().yaxis.set_minor_locator(plt.MultipleLocator(0.5))  # Modifica la distanza tra i minor ticks sull'asse Y


# Aggiungi la legenda
plt.legend()

# Mostra il grafico
plt.show()

# Esporta il grafico in un file png
plt.savefig("errors_plot.png")