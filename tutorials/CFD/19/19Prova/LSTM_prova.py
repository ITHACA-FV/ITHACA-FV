import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.models import Sequential, load_model
from tensorflow.keras.layers import LSTM, Dense, Dropout
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
import os


###################    CREAZIONE MODELLO RETE NEURALE LSTM     ###################


# === CREA CARTELLE RISULTATI ===
os.makedirs("Prova/plots", exist_ok=True)

# ==== PARAMETRI ====
lookback = 5        # fino a 20
epochs = 1000       # epoche
batch_size = 32     # fino a 64
dt = 1.0            # Timestep
lambda_phys = 0.01  # Peso della loss fisica

# ==== CARICAMENTO DATI POD ====
coeffs_u = np.load("./Coeffs/u_coeffs.npy")     # shape (T, r_u)
coeffs_p = np.load("./Coeffs/p_coeffs.npy")     # shape (T, r_p)
coeffs_nut = np.load("./Coeffs/nut_coeffs.npy") # shape (T, r_nut)

# ==== CONSIDERARE SOLO E PRIME 1800 COLONNE (SNAPSHOTS) PER IL TRAINING ====
u_coeffs_1800 = coeffs_u[:, :1800]
p_coeffs_1800 = coeffs_p[:, :1800]
nut_coeffs_1800 = coeffs_nut[:, :1800]

# ==== CONSIDERARE LE ULTIME 201 COLONNE (SNAPSHOTS) PER LA VALIDATION ====
u_coeffs_last201 = coeffs_u[:, -201:]
p_coeffs_last201 = coeffs_p[:, -201:]
nut_coeffs_last201 = coeffs_nut[:, -201:]

# ==== CONCATENAZIONE INPUT ====
X_all = np.hstack([u_coeffs_1800.T, p_coeffs_1800.T])             # Training
Y_all = nut_coeffs_1800.T                                         # Training

X_val_raw = np.hstack([u_coeffs_last201.T, p_coeffs_last201.T])   # Validation
Y_val_raw = nut_coeffs_last201.T                                  # Validation

# ==== NORMALIZZAZIONE PER TRAINING E VALIDAZIONE ====
x_scaler = StandardScaler()
y_scaler = StandardScaler()

X_all = x_scaler.fit_transform(X_all)                         # Training
Y_all = y_scaler.fit_transform(Y_all)                         # Training

X_val = x_scaler.transform(X_val_raw)                         # Validazione
Y_val = y_scaler.transform(Y_val_raw)                         # Validazione

# ==== CREAZIONE SEQUENZE PER LSTM ====
def create_sequences(X, Y, lookback, step=1):
    X_seq, Y_seq = [], []
    for i in range(0, len(X) - lookback, step):
        X_seq.append(X[i:i+lookback])
        Y_seq.append(Y[i+lookback])
    return np.array(X_seq), np.array(Y_seq)

X_train_seq, Y_train_seq = create_sequences(X_all, Y_all, lookback, step=1)       # Training

X_val_seq, Y_val_seq = create_sequences(X_val, Y_val, lookback=lookback, step=1)  # Validation

print("X_seq shape:", X_train_seq.shape)
print("Y_seq shape:", Y_train_seq.shape)

# ==== COEFF. PRESSIONE ====
X_train_p_seq = p_coeffs_1800.T[lookback:]
X_val_p_seq = p_coeffs_last201.T[lookback:]

# ==== MATRICI PROIEZIONI GALERKIN ====
B = np.load("./ITHACAoutput/Matrices/B_0_10_0.npy")
BP = np.load("./ITHACAoutput/Matrices/BP_10.npy")
C = np.load("./ITHACAoutput/Matrices/C_0_10_0_t.npy")
Cf = np.load("./ITHACAoutput/Matrices/Cf_0_10_0_t.npy")
Ci = np.load("./ITHACAoutput/Matrices/Ci_0_10_0_t.npy")
DF = np.load("./ITHACAoutput/Matrices/DF_10_0.npy")
I = np.load("./ITHACAoutput/Matrices/I_10.npy")
K = np.load("./ITHACAoutput/Matrices/K_0_10_0_10.npy")
KF = np.load("./ITHACAoutput/Matrices/KF_10_0.npy")
P = np.load("./ITHACAoutput/Matrices/P_0_10_0_10.npy")
W = np.load("./ITHACAoutput/Matrices/W_10.npy")
RC = np.load("./ITHACAoutput/Matrices/RC/RC0_10_0.npy")
RD = np.load("./ITHACAoutput/Matrices/RD/RD0_10_0.npy")
SC = np.load("./ITHACAoutput/Matrices/SC/SC0_10_0.npy")
SD = np.load("./ITHACAoutput/Matrices/SD/SD0_10_0.npy")

# ==== CONVERSIONE IN TF ====
B_tf = tf.constant(B, dtype=tf.float32)
BP_tf = tf.constant(BP, dtype=tf.float32)
C_tf = tf.constant(C, dtype=tf.float32)
Cf_tf = tf.constant(Cf, dtype=tf.float32)
Ci_tf = tf.constant(Ci, dtype=tf.float32)
DF_tf = tf.constant(DF, dtype=tf.float32)
I_tf = tf.constant(I, dtype=tf.float32)
K_tf = tf.constant(K, dtype=tf.float32)
KF_tf = tf.constant(KF, dtype=tf.float32)
P_tf = tf.constant(P, dtype=tf.float32)
W_tf = tf.constant(W, dtype=tf.float32)
RC_tf = tf.constant(RC, dtype=tf.float32)
RD_tf = tf.constant(RD, dtype=tf.float32)
SC_tf = tf.constant(SC, dtype=tf.float32)
SD_tf = tf.constant(SD, dtype=tf.float32)

# ==== GALERKIN RHS ====
def galerkin_rhs(a, p):
    B_term   = tf.einsum('ij,bj->bi', B_tf, a)     
    BP_term  = tf.einsum('ij,bj->bi', BP_tf, p)
    DF_term  = tf.einsum('ij,bj->bi', DF_tf, a)
    C_term   = tf.einsum('ijk,bj,bk->bi', C_tf, a, a)
    Cf_term  = tf.einsum('ijk,bj,bk->bi', Cf_tf, a, a)
    Ci_term  = tf.einsum('ijk,bj,bk->bi', Ci_tf, a, a)
    K_term   = tf.einsum('ij,bj->bi', K_tf, a)
    P_term   = tf.einsum('ij,bj->bi', P_tf, p)

    # batch_size = tf.shape(a)[0]
    # RC_term = tf.tile(tf.reshape(RC_tf, (1, -1)), [batch_size, 1])
    # RD_term = tf.tile(tf.reshape(RD_tf, (1, -1)), [batch_size, 1])
    # SC_term = tf.tile(tf.reshape(SC_tf, (1, -1)), [batch_size, 1])
    # SD_term = tf.tile(tf.reshape(SD_tf, (1, -1)), [batch_size, 1])

    return B_term + BP_term + DF_term + C_term + Cf_term + Ci_term + K_term + P_term  # + RC_term + RD_term + SC_term + SD_term

# ==== CUSTOM LOSS ====
def make_custom_loss(p_seq_tf):
    def loss(y_true, y_pred):
        # p_seq_batch va selezionato in base al batch
        batch_size = tf.shape(y_pred)[0]
        p_seq_batch = p_seq_tf[:batch_size]

        # Calcola la loss fisica
        galerkin_pred = galerkin_rhs(y_pred, p_seq_batch)
        galerkin_true = galerkin_rhs(y_true, p_seq_batch)
        mse_data = tf.reduce_mean(tf.square(y_true - y_pred))
        mse_phys = tf.reduce_mean(tf.square(galerkin_pred - galerkin_true))
        return mse_data + lambda_phys * mse_phys
    return loss

# def custom_loss(y_true, y_pred, p_seq):

#     tf.print("y_true shape:", tf.shape(y_true))
#     tf.print("y_pred shape:", tf.shape(y_pred))
#     tf.print("p_seq shape:", tf.shape(p_seq))


#     mse_data = tf.reduce_mean(tf.square(y_true - y_pred))
#     a_t = y_pred[:-1]
#     a_tp1 = y_pred[1:]
#     a_dot_pred = (a_tp1 - a_t) / dt

#     batch_size = tf.shape(y_pred)[0]
#     p_seq_batch = p_seq_tf[:batch_size]       
#     p_seq_t = p_seq_batch[:-1]           
    
#     f_rhs = galerkin_rhs(a_t, p_seq[:-1])
#     mse_phys = tf.reduce_mean(tf.square(a_dot_pred - f_rhs))
#     return mse_data + lambda_phys * mse_phys

# ==== WRAPPER ====
# def loss_wrapper(p_seq_tf):
#     def loss_fn(y_true, y_pred):
#         return custom_loss(y_true, y_pred, p_seq_tf)
#     return loss_fn

p_seq_tf = tf.constant(X_train_p_seq, dtype=tf.float32)

# ==== DEFINIZIONE MODELLO LSTM ====
model = Sequential()
model.add(LSTM(32, return_sequences=True, input_shape=(lookback, X_train_seq.shape[2])))
model.add(LSTM(16, return_sequences=False))
model.add(Dense(Y_train_seq.shape[1]))

opt = Adam(learning_rate=0.001)
# model = Sequential()
# model.add(LSTM(64, return_sequences=True, input_shape=(lookback, X_train_seq.shape[2])))  # 1° LSTM
# model.add(LSTM(32, return_sequences=False))                                               # 2° LSTM
# model.add(Dense(32, activation='relu'))                                                   # 3° Dense
# model.add(Dense(Y_train_seq.shape[1]))                                                    # Output layer

model.compile(optimizer=opt, loss=make_custom_loss(p_seq_tf))
# model.compile(optimizer='adam', loss=loss_wrapper(p_seq_tf))
model.summary()


###################    TRAINING     ###################

history = model.fit(X_train_seq, Y_train_seq, epochs=epochs, batch_size=batch_size, validation_data=(X_val_seq, Y_val_seq))

# ==== PREDIZIONE DOPO IL TRAINING ====
Y_pred = model.predict(X_train_seq)
Y_pred_original = y_scaler.inverse_transform(Y_pred)
Y_true_original = y_scaler.inverse_transform(Y_train_seq)

# ==== SALVATAGGIO MODELLO ====
model.save("Prova/trained_model.keras")

# ==== PLOT ERRORE Relativo TRA VALORE PREDETTO E QUELLO ORIGINALE DOPO IL TRAINING ====
for i in range(Y_true_original.shape[1]): 
    abs_error_train = np.abs(Y_true_original[:, i] - Y_pred_original[:, i])/np.abs(Y_pred_original[:, i])

    plt.figure(figsize=(10, 4))
    plt.plot(abs_error_train, label=f'|Errore training|')
    plt.title(f'Errore Relativo (Training)')
    plt.xlabel("Time step")
    plt.ylabel("Errore relativo")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(f"Prova/plots/rel_error_train.png")
    plt.show()
    plt.close()

# ==== PLOT ERRORE Assoluto TRA VALORE PREDETTO E QUELLO ORIGINALE DOPO IL TRAINING ====
for i in range(min(3, Y_true_original.shape[1])):  
    rel_error_train = np.abs(Y_true_original[:, i] - Y_pred_original[:, i])

    plt.figure(figsize=(10, 4))
    plt.plot(rel_error_train, label=f'|Errore training|')
    plt.title(f'Errore Assoluto (Training)')
    plt.xlabel("Time step")
    plt.ylabel("Errore assoluto")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(f"Prova/plots/abs_error_train.png")
    plt.show()
    plt.close()



###################    VALIDAZIONE     ###################

Y_val_pred = model.predict(X_val_seq)
a_t_val = Y_val_pred[:-1]
a_tp1_val = Y_val_pred[1:]
a_dot_val = (a_tp1_val - a_t_val) / dt
p_val_t = tf.constant(X_val_p_seq[:-1], dtype=tf.float32)
f_rhs_val = galerkin_rhs(tf.constant(a_t_val), p_val_t)
physics_loss_val = tf.reduce_mean(tf.square(a_dot_val - f_rhs_val))
print("Validation Physics Loss:", physics_loss_val.numpy())


Y_val_pred_orig = y_scaler.inverse_transform(Y_val_pred)
Y_val_true_orig = y_scaler.inverse_transform(Y_val_seq)

# === SALVATAGGI ===
np.save("Prova/Y_val_pred.npy", Y_val_pred_orig)
np.save("Prova/Y_val_true.npy", Y_val_true_orig)

# === SALVA HISTORY (loss e val_loss) ===
np.save("Prova/training_history.npy", history.history)

plt.figure(figsize=(8, 4))
plt.plot(history.history['loss'], label='Train Loss')
plt.plot(history.history['val_loss'], label='Val Loss')
plt.xlabel("Epoch")
plt.ylabel("Loss (MSE)")
plt.title("Training History")
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("Prova/plots/training_history.png")
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
np.savetxt("Prova/validation_errors.csv", errors_array, delimiter=",", header="mse,mae", comments='')

# === ERRORE PER ISTANTE TEMPORALE (SNAPSHOT) ===
snapshot_mse = np.mean((Y_val_true_orig - Y_val_pred_orig)**2, axis=1)
np.save("Prova/snapshot_mse.npy", snapshot_mse)

plt.figure(figsize=(10, 4))
plt.plot(snapshot_mse)
plt.title("Snapshot-wise MSE (Validation Set)")
plt.xlabel("Time step")
plt.ylabel("MSE")
plt.grid()
plt.tight_layout()
plt.savefig("Prova/plots/snapshot_mse.png")
plt.close()

# === PLOT ERRORE Relativo DOPO LA VALIDAZIONE ===
for i in range(Y_val_true_orig.shape[1]):
    rel_error_val = np.abs(Y_val_true_orig[:, i] - Y_val_pred_orig[:, i])/np.abs(Y_val_pred_orig[:, i])

    plt.figure(figsize=(10, 4))
    plt.plot(rel_error_val, label='|True - Pred|')
    plt.title(f'Errore Relativo (Validation)')
    plt.xlabel("Time step")
    plt.ylabel("Errore Relativo")
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(f"Prova/plots/rel_error_val.png")
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
    plt.savefig(f"Prova/plots/abs_error_val.png")
    plt.close()