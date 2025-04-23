import numpy as np
import matplotlib.pyplot as plt

# Carica i dati dai file estratti
# Ux = np.loadtxt("logs/residuals_0")
# Uy = np.loadtxt("logs/residuals_1")
# Uz = np.loadtxt("logs/residuals_2")
p_0 = np.loadtxt("logs/p_0")
pFinalRes_0 = np.loadtxt("logs/pFinalRes_0")

# Plot dei residui
plt.figure(figsize=(8,6))
# plt.semilogy(Ux[:,0], Ux[:,1], label="Ux Residual")
# plt.semilogy(Uy[:,0], Uy[:,1], label="Uy Residual")
# plt.semilogy(Uz[:,0], Uz[:,1], label="Uz Residual")
plt.semilogy(p_0[:,0], p_0[:,1], label="Pressure Residual")
plt.semilogy(pFinalRes_0[:,0], pFinalRes_0[:,1], label="Final Pressure Residual")

plt.xlabel("Time Step")
plt.ylabel("Residuals")
plt.legend()
# plt.grid()

# Abilita i minor ticks (i "piccoli ticks")
plt.minorticks_on()

# Aggiungi la griglia maggiore e minore
plt.grid(True, which='both', linestyle='-', color='gray', linewidth=0.5)

# Puoi personalizzare la densità dei "minor ticks" per avere rettangoli più piccoli
plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(1))  # Modifica la distanza tra i minor ticks sull'asse X
plt.gca().yaxis.set_minor_locator(plt.MultipleLocator(0.5))  # Modifica la distanza tra i minor ticks sull'asse Y

plt.title("Convergence Residuals")
plt.show()
plt.savefig("residuals_plot.png")

