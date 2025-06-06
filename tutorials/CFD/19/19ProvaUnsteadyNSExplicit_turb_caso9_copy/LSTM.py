import numpy as np
import matplotlib.pyplot as plt
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import LSTM, Dense, Dropout
from sklearn.preprocessing import StandardScaler


# ==== PARAMETRI ====
lookback = 5     # fino a 20
epochs = 1000
batch_size = 16   # fino a 64

# ==== CARICAMENTO DATI POD ====
POD_u = np.load("./Matrix_py/u_matrix.npy")     # shape (T, r_u)
POD_p = np.load("./Matrix_py/p_matrix.npy")     # shape (T, r_p)
POD_phi = np.load("./Matrix_py/phi_matrix.npy") # shape (T, r_phi)
POD_nut = np.load("./Matrix_py/nut_matrix.npy") # shape (T, r_nut)

# ==== CONCATENAZIONE INPUT ====
X_all = np.hstack([POD_u, POD_p, POD_phi])  # shape (T, total_input_dim)
Y_all = POD_nut  # target

# ==== NORMALIZZAZIONE ====
x_scaler = StandardScaler()
y_scaler = StandardScaler()

X_all = x_scaler.fit_transform(X_all)
Y_all = y_scaler.fit_transform(Y_all)

# ==== CREAZIONE SEQUENZE PER LSTM ====
def create_sequences(X, Y, lookback, step=1):
    X_seq, Y_seq = [], []
    for i in range(0, len(X) - lookback, step):
        X_seq.append(X[i:i+lookback])
        Y_seq.append(Y[i+lookback])
    return np.array(X_seq), np.array(Y_seq)

X_seq, Y_seq = create_sequences(X_all, Y_all, lookback, step=1)

print("X_seq shape:", X_seq.shape)
print("Y_seq shape:", Y_seq.shape)


# ==== DEFINIZIONE MODELLO LSTM ====
model = Sequential()
model.add(LSTM(64, return_sequences=True, input_shape=(lookback, X_seq.shape[2])))  # 1° LSTM
model.add(LSTM(32, return_sequences=False))                                         # 2° LSTM
model.add(Dense(32, activation='relu'))                                             # 3° Dense
model.add(Dense(Y_seq.shape[1]))                                                    # Output layer

model.compile(optimizer='adam', loss='mse')
model.summary()

# ==== TRAINING ====
history = model.fit(X_seq, Y_seq, epochs=epochs, batch_size=batch_size, validation_split=0.2)

# ==== PREDIZIONE ====
Y_pred = model.predict(X_seq)
Y_pred_original = y_scaler.inverse_transform(Y_pred)
Y_true_original = y_scaler.inverse_transform(Y_seq)

# ==== PLOT COEFFICIENTI ====
import matplotlib.pyplot as plt
for i in range(min(3, Y_true_original.shape[1])):  # plottiamo i primi 3 modi
    plt.figure(figsize=(10, 4))
    plt.plot(Y_true_original[:, i], label=f'True α_nut_{i+1}')
    plt.plot(Y_pred_original[:, i], '--', label=f'Pred α_nut_{i+1}')
    plt.legend()
    plt.title(f'Confronto α_nut_{i+1}')
    plt.grid()
    plt.show()

# ==== (OPZIONALE) RICOSTRUZIONE CAMPI 3D ====
# Carica i modi nut e ricostruisci il campo (solo se necessario)
# try:
#     modes_nut = np.load("modes_nut.npy")  # shape (r_nut, Nx, Ny, Nz)
#     T_pred = Y_pred_original.shape[0]
#     Nx, Ny, Nz = modes_nut.shape[1:]
#     nut_reconstructed = np.zeros((T_pred, Nx, Ny, Nz))

#     for t in range(T_pred):
#         for r in range(modes_nut.shape[0]):
#             nut_reconstructed[t] += Y_pred_original[t, r] * modes_nut[r]
    
#     print("Ricostruzione nut completata. Shape:", nut_reconstructed.shape)
# except FileNotFoundError:
#     print("File modes_nut.npy non trovato: salto ricostruzione campi.")

