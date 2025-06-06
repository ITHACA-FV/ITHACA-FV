import numpy as np
import matplotlib.pyplot as plt
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.models import Sequential, load_model
from tensorflow.keras.layers import LSTM, Dense, Dropout
from tensorflow.keras.callbacks import EarlyStopping
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
import os


###################    CREAZIONE MODELLO RETE NEURALE LSTM     ###################


# === CREA CARTELLE RISULTATI ===
os.makedirs("Eigen/plots", exist_ok=True)

# ==== PARAMETRI ====
lookback = 15     # fino a 20
epochs = 1300
batch_size = 64   # fino a 64

# ==== CARICAMENTO DATI POD ====
coeffs_u = np.load("./Coeffs/u_coeffs.npy")     # shape (T, r_u)
coeffs_p = np.load("./Coeffs/p_coeffs.npy")     # shape (T, r_p)
coeffs_nut = np.load("./Coeffs/nut_coeffs.npy") # shape (T, r_nut)

eigen_u = np.loadtxt("./ITHACAoutput/POD/CumEigenvalues_U", skiprows=2)
eigen_p = np.loadtxt("./ITHACAoutput/POD/CumEigenvalues_p", skiprows=2)
eigen_nut = np.loadtxt("./ITHACAoutput/POD/CumEigenvalues_nut", skiprows=2)

# ==== CONSIDERARE SOLO E PRIME 1800 COLONNE (SNAPSHOTS) PER IL TRAINING ====
u_coeffs_1800 = coeffs_u[:, :1800]
p_coeffs_1800 = coeffs_p[:, :1800]
nut_coeffs_1800 = coeffs_nut[:, :1800]

first_eigen_u = eigen_u[:1800]
first_eigen_p = eigen_p[:1800]
first_eigen_nut = eigen_nut[:1800]

train_coeff_u = u_coeffs_1800 * first_eigen_u
train_coeff_p = p_coeffs_1800 * first_eigen_p
train_coeff_nut = nut_coeffs_1800 * first_eigen_nut

print("train_coeff_u shape", train_coeff_u.shape)
print("train_coeff_p shape", train_coeff_p.shape)
print("train_coeff_nut shape", train_coeff_nut.shape)

# ==== CONSIDERARE LE ULTIME 201 COLONNE (SNAPSHOTS) PER LA VALIDATION ====
u_coeffs_last201 = coeffs_u[:, -201:]
p_coeffs_last201 = coeffs_p[:, -201:]
nut_coeffs_last201 = coeffs_nut[:, -201:]

final_eigen_u = eigen_u[-201:]
final_eigen_p = eigen_p[-201:]
final_eigen_nut = eigen_nut[-201:]

val_coeff_u = u_coeffs_last201 * final_eigen_u
val_coeff_p = p_coeffs_last201 * final_eigen_p
val_coeff_nut = nut_coeffs_last201 * final_eigen_nut

print("val_coeff_u shape", val_coeff_u.shape)
print("val_coeff_p shape", val_coeff_p.shape)
print("val_coeff_nut shape", val_coeff_nut.shape)

# ==== CONCATENAZIONE INPUT ====
X_all = np.hstack([train_coeff_u.T, train_coeff_p.T])  
Y_all = train_coeff_nut.T  

# X_all = np.hstack([u_coeffs_1800.T, p_coeffs_1800.T])  
# Y_all = nut_coeffs_1800.T  

# ==== NORMALIZZAZIONE PER TRAINING E VALIDAZIONE ====
x_scaler = StandardScaler()
y_scaler = StandardScaler()

X_all = x_scaler.fit_transform(X_all)                         # Training
Y_all = y_scaler.fit_transform(Y_all)                         # Training

X_val_raw = np.hstack([val_coeff_u.T, val_coeff_p.T])   # Validazione
Y_val_raw = val_coeff_nut.T                                # Validazione

# X_val_raw = np.hstack([u_coeffs_last201.T, p_coeffs_last201.T])   # Validazione
# Y_val_raw = nut_coeffs_last201.T                                # Validazione

X_val = x_scaler.transform(X_val_raw)                         # Validazione
Y_val = y_scaler.transform(Y_val_raw)                         # Validazione

# ==== CREAZIONE SEQUENZE PER LSTM ====
def create_sequences(X, Y, lookback, step=1):
    X_seq, Y_seq = [], []
    for i in range(0, len(X) - lookback, step):
        X_seq.append(X[i:i+lookback])
        Y_seq.append(Y[i+lookback])
    return np.array(X_seq), np.array(Y_seq)

X_train_seq, Y_train_seq = create_sequences(X_all, Y_all, lookback, step=1)

X_val_seq, Y_val_seq = create_sequences(X_val, Y_val, lookback=lookback, step=1)

print("X_seq shape:", X_train_seq.shape)
print("Y_seq shape:", Y_train_seq.shape)

# ==== DEFINIZIONE MODELLO LSTM ====
# model = Sequential()
# model.add(LSTM(64, return_sequences=False, input_shape=(lookback, X_train_seq.shape[2])))
# model.add(Dense(32, activation='relu'))  # Hidden Dense
# model.add(Dense(Y_train_seq.shape[1]))   # Output

model = Sequential()
model.add(LSTM(64, return_sequences=True, input_shape=(lookback, X_train_seq.shape[2])))  # 1° LSTM
model.add(LSTM(32, return_sequences=False))                                               # 2° LSTM
model.add(Dense(32, activation='relu'))                                                   # 3° Dense
model.add(Dense(Y_train_seq.shape[1]))                                                    # Output layer

opt = Adam(learning_rate=1.6e-5)

model.compile(optimizer=opt, loss='mse')
model.summary()

early_stop = EarlyStopping(
    monitor='val_loss',        # monitora la loss di validazione
    patience=5,                # numero di epoche senza miglioramento prima di fermare
    restore_best_weights=True # ripristina i pesi migliori al termine
)


###################    TRAINING     ###################


history = model.fit(X_train_seq, Y_train_seq, epochs=epochs, batch_size=batch_size, validation_data=(X_val_seq, Y_val_seq))

# ==== PREDIZIONE DOPO IL TRAINING ====
Y_pred = model.predict(X_train_seq)
Y_pred_original = y_scaler.inverse_transform(Y_pred)
Y_true_original = y_scaler.inverse_transform(Y_train_seq)

# ==== SALVATAGGIO MODELLO ====
model.save("./Eigen/trained_model.keras")

# ==== PLOT ERRORE Relativo TRA VALORE PREDETTO E QUELLO ORIGINALE DOPO IL TRAINING ====
for i in range(Y_true_original.shape[1]): 
    rel_error_train = np.abs(Y_true_original[:, i] - Y_pred_original[:, i])/np.abs(Y_pred_original[:, i])

    plt.figure(figsize=(10, 4))
    plt.plot(rel_error_train, label=f'|Errore training|')
    plt.title(f'Errore Relativo (Training)')
    plt.xlabel("Time step")
    plt.ylabel("Errore relativo")
    # plt.yscale("log") 
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(f"Eigen/plots/rel_error_train.png")
    plt.show()
    plt.close()

# ==== PLOT ERRORE Assoluto TRA VALORE PREDETTO E QUELLO ORIGINALE DOPO IL TRAINING ====
for i in range(min(3, Y_true_original.shape[1])):  
    abs_error_train = np.abs(Y_true_original[:, i] - Y_pred_original[:, i])

    plt.figure(figsize=(10, 4))
    plt.plot(abs_error_train, label=f'|Errore training|')
    plt.title(f'Errore Assoluto (Training)')
    plt.xlabel("Time step")
    plt.ylabel("Errore assoluto")
    # plt.yscale("log")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(f"Eigen/plots/abs_error_train.png")
    plt.show()
    plt.close()



###################    VALIDAZIONE     ###################

# X_val_seq, Y_val_seq = create_sequences(X_val, Y_val, lookback=lookback, step=1)

Y_val_pred = model.predict(X_val_seq)
Y_val_pred_orig = y_scaler.inverse_transform(Y_val_pred)
Y_val_true_orig = y_scaler.inverse_transform(Y_val_seq)

# === SALVA PREDICTION E VERITÀ ===
np.save("Eigen/Y_val_pred.npy", Y_val_pred_orig)
np.save("Eigen/Y_val_true.npy", Y_val_true_orig)

# === SALVA HISTORY (loss e val_loss) ===
np.save("Eigen/training_history.npy", history.history)

plt.figure(figsize=(8, 4))
plt.plot(history.history['loss'], label='Train Loss')
plt.plot(history.history['val_loss'], label='Val Loss')
plt.xlabel("Epoch")
plt.ylabel("Loss (MSE)")
plt.title("Training History")
plt.yscale("log") 
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("Eigen/plots/training_history.png")
plt.close()

# === METRICHE PER OGNI MODO POD ===
errors = {
    'mse': [],
    'mae': [],
}

for i in range(Y_val_true_orig.shape[1]):
    true_vals = Y_val_true_orig[:, i]
    pred_vals = Y_val_pred_orig[:, i]

    errors['mse'].append(mean_squared_error(true_vals, pred_vals))
    errors['mae'].append(mean_absolute_error(true_vals, pred_vals))

errors_array = np.column_stack((errors['mse'], errors['mae']))
np.savetxt("Eigen/validation_errors.csv", errors_array, delimiter=",", header="mse,mae", comments='')

# === ERRORE PER ISTANTE TEMPORALE (SNAPSHOT) ===
snapshot_mse = np.mean((Y_val_true_orig - Y_val_pred_orig)**2, axis=1)
np.save("Eigen/snapshot_mse.npy", snapshot_mse)

plt.figure(figsize=(10, 4))
plt.plot(snapshot_mse)
plt.title("Snapshot-wise MSE (Validation Set)")
plt.xlabel("Time step")
plt.ylabel("MSE")
plt.grid()
plt.tight_layout()
plt.savefig("Eigen/plots/snapshot_mse.png")
plt.close()

# === PLOT ERRORE Relativo DOPO LA VALIDAZIONE ===
for i in range(Y_val_true_orig.shape[1]):
    rel_error_val = np.abs(Y_val_true_orig[:, i] - Y_val_pred_orig[:, i])/np.abs(Y_val_pred_orig[:, i])

    plt.figure(figsize=(10, 4))
    plt.plot(rel_error_val, label='|True - Pred|')
    plt.title(f'Errore Relativo (Validation)')
    plt.xlabel("Time step")
    plt.ylabel("Errore Relativo") 
    plt.yscale("log")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(f"Eigen/plots/rel_error_val.png")
    plt.close()

# === PLOT ERRORE Assoluto DOPO LA VALIDAZIONE ===
for i in range(min(3, Y_val_true_orig.shape[1])):
    abs_error_val = np.abs(Y_val_true_orig[:, i] - Y_val_pred_orig[:, i])

    plt.figure(figsize=(10, 4))
    plt.plot(abs_error_val, label='|True - Pred|')
    plt.title(f'Errore Assoluto (Validation)')
    plt.xlabel("Time step")
    plt.ylabel("Errore Assoluto") 
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(f"Eigen/plots/abs_error_val.png")
    plt.close()