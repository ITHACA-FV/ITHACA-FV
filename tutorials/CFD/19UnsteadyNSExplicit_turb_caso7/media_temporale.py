import numpy as np
import matplotlib.pyplot as plt


error = np.load('error_10_10.npy')

# Calcola la media nel tempo 
mean_error = np.mean(error)

# Stampa la media
print(f"La media dell'errore nel tempo è: {mean_error}")

# Se vuoi visualizzarlo, puoi anche fare un grafico con la media (dovresti ripetere la media per ogni time step, se necessario)
plt.plot(error, label="Errore per time step")
plt.axhline(mean_error, color='red', linestyle='--', label="Media dell'errore")

plt.xlabel("Time Step")
plt.ylabel("Errore")
plt.legend()
plt.title("Errore nel tempo con la sua media")
plt.show()
