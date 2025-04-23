import numpy as np
from scipy.linalg import solve
from ReducedUnsteadyNSExplicit import ReducedUnsteadyNSExplicit


# Definire il metodo di flusso da usare: "consistent" o "inconsistent"
flux_method = "consistent"

# Velocità di input (può essere float o array, a seconda del modello)
vel = 1.0

# Istanziare il solver
solver = ReducedUnsteadyNSExplicit(flux_method)

# Esecuzione della simulazione
solver.solve_online(vel)

# Filtro per rimuovere i None da online_solution
# cleaned_solution = [arr for arr in solver.online_solution if arr is not None]

# Verifica che tutte le righe abbiano la stessa shape
# shapes = [arr.shape for arr in cleaned_solution]
shapes = [arr.shape for arr in solver.online_solution if arr is not None]

for i, arr in enumerate(solver.online_solution):
    print(f"Elemento {i}: shape = {np.shape(arr)}")

# Se tutte le righe hanno la stessa shape, si salva come array 2D
if all(s == shapes[0] for s in shapes):
    # online_solution_array = np.array(cleaned_solution)
    # np.save("online_solution_" + str(flux_method) + ".npy", online_solution_array)

    online_solution_array = np.array([arr for arr in solver.online_solution if arr is not None])
    np.save("online_solution_" + str(flux_method) + ".npy", online_solution_array)
else:
    # Altrimenti si salva come oggetto (dtype=object)
    # np.save("online_solution_" + str(flux_method) + ".npy", np.array(cleaned_solution, dtype=object))
    np.save("online_solution_" + str(flux_method) + ".npy", np.array([arr for arr in solver.online_solution if arr is not None], dtype=object))


