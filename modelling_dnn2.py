# -*- coding: utf-8 -*-
"""
Takes as input peptides_array.npy

achieves 80% with HemoPi2
"""
from IPython import get_ipython
get_ipython().magic('reset -sf') 
import numpy as np
from keras.models import Sequential
from keras.layers import Dense, Activation, Dropout, Conv1D, Embedding, MaxPooling1D, LSTM, Merge, GRU, Bidirectional
from sklearn.model_selection import train_test_split
from keras.callbacks import EarlyStopping, ModelCheckpoint, TensorBoard
from keras import optimizers
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.model_selection import StratifiedKFold
from aa_indeces import aai_to_get

from tensorflow.python.client import device_lib

print(device_lib.list_local_devices())



inp = np.load("inputs/peptides_array_7_array.npy", allow_pickle=False)

X = np.array([i for i in inp])

Y = np.array([i[-1:] for i in inp])

kfold = StratifiedKFold(n_splits=5, shuffle=True, random_state=441653 +1)
cvsscores = []

####Scaling#####

scaler = MinMaxScaler(feature_range=(0, 1))
X_ann = scaler.fit_transform(X)

scaler = StandardScaler(with_mean=True, with_std=True)
scaler.fit(X_ann)
X_ann = scaler.transform(X_ann)

#print(X)
#######
X_rnn_0 = np.load("inputs/peptides_array_7_alpha2num_veltri.npy")
X_rnn_1 = np.load("inputs/peptides_array_7_pos_first_aa_1.npy")
X_rnn_2 = np.load("inputs/peptides_array_7_pos_first_veltri.npy")
X_rnn_3 = np.load("inputs/peptides_array_7_pos_last_aa_1.npy")
X_rnn_4 = np.load("inputs/peptides_array_7_pos_last_veltri.npy")
X_rnn_5 = np.load("inputs/peptides_array_7_property_pos_first_1.npy")
X_rnn_6 = np.load("inputs/peptides_array_7_property_pos_last_1.npy")

embed_vocab_0 = int(max(list(X_rnn_0.flatten()))) + 1
embed_vocab_1 = int(max(list(X_rnn_1.flatten()))) + 1
embed_vocab_2 = int(max(list(X_rnn_2.flatten()))) + 1
embed_vocab_3 = int(max(list(X_rnn_3.flatten()))) + 1
embed_vocab_4 = int(max(list(X_rnn_4.flatten()))) + 1

dims = X.shape[1]

X_rnn_aai_x = np.load("inputs/peptides_array_7_property_pos_convx.npy")
X_rnn_aai_shape = X_rnn_aai_x.shape

print(dims, 'dims')
print("Building model...")

nb_classes = Y.shape[1]
print(nb_classes, 'classes')

rate = 0.1

iteration = 0
for train, test in kfold.split(X_ann, Y):
    
    dims_ann = X_ann[train].shape[1]
    
    models = []
    for embed_vocab, X_rnn, inp_len, kernel_size, num_conv in [
            (embed_vocab_0, X_rnn_0, 200, 2, 4), 
            (embed_vocab_1, X_rnn_1, 15, 2, 2), 
            (embed_vocab_2, X_rnn_2, 15, 2, 2), 
            (embed_vocab_3, X_rnn_3, 15, 2, 2), 
            (embed_vocab_4, X_rnn_4, 15, 2, 2), ]:
        print(inp_len)
        model_rnn = Sequential()
        #model.add(Dense(1024, input_shape=(dims,), init='uniform', activation='relu'))
        #model.add(Dropout(0.94, noise_shape=None, seed=42))
        model_rnn.add(Embedding(embed_vocab, 128, input_length=inp_len))
        for i_num_conv in range(num_conv):
            model_rnn.add(Conv1D(64, kernel_size, strides=1, init="normal"))
            model_rnn.add(MaxPooling1D(pool_size=2, strides=2))
            model_rnn.add(Dropout(0.5))
            
        model_rnn.add(GRU(32, unroll = True, stateful = False, dropout = 0.5))
        #model_rnn.add(Dense(1, init='normal', activation="sigmoid"))
        models.append(model_rnn)
    
    model_aai = Sequential()
    model_aai.add(Conv1D(64, 1, strides=1, init="normal", activation="tanh", input_shape=(X_rnn_aai_shape[1], X_rnn_aai_shape[2])))  # X_rnn_aai_shape[1:]
    model_aai.add(MaxPooling1D(pool_size=2, strides=2))
    model_aai.add(Dropout(0.2))
    model_aai.add(Conv1D(64, 1, strides=1, init="normal", activation="tanh"))  # X_rnn_aai_shape[1:]
    model_aai.add(MaxPooling1D(pool_size=2, strides=2))
    model_aai.add(Dropout(0.2))
    model_aai.add(Conv1D(64, 1, strides=1, init="normal", activation="tanh"))  # X_rnn_aai_shape[1:]
    model_aai.add(MaxPooling1D(pool_size=2, strides=2))
    model_aai.add(Dropout(0.2))
    model_aai.add(Conv1D(64, 1, strides=1, init="normal", activation="tanh"))  # X_rnn_aai_shape[1:]
    model_aai.add(MaxPooling1D(pool_size=2, strides=2))
    model_aai.add(Dropout(0.2))
    model_aai.add(Conv1D(64, 1, strides=1, init="normal", activation="tanh"))  # X_rnn_aai_shape[1:]
    model_aai.add(MaxPooling1D(pool_size=2, strides=2))
    model_aai.add(Dropout(0.2))
    model_aai.add(Bidirectional(GRU(32, unroll = True, stateful = False, dropout = 0.88)))
    models.append(model_aai)

    model_ann = Sequential()
    model_ann.add(Dense(1024, input_shape=(dims,), init='uniform', activation='relu'))
    model_ann.add(Dropout(0.85, input_shape=(dims_ann,), noise_shape=None, seed=42))
    model_ann.add(Dense(64,  init='normal', activation='relu'))
    model_ann.add(Dropout(0.85, noise_shape=None, seed=42))
    model_ann.add(Dense(32,  init='normal', activation='relu'))
    model_ann.add(Dropout(0.8, seed=42, noise_shape=None))
    models.append(model_ann)
        
    model = Sequential()
    model.add(Merge(models, mode="concat"))
    model.add(Dense(16, init="normal", activation="sigmoid"))
    model.add(Dropout(0.4, seed=42, noise_shape=None))
    model.add(Dense(1, init='normal', activation="sigmoid"))
    
    sgd = optimizers.SGD(lr=0.01, decay=0.005, momentum=0.5, nesterov=True)
    rmsprop = optimizers.RMSprop(lr=0.02, rho=0.9, epsilon=1e-08, decay=0.01)
    #adagrad = optimizers.Adagrad(lr=0.005)
    adagrad = optimizers.Adagrad(lr=0.005)
    adadelta = optimizers.Adadelta(lr=0.01, rho=0.95, epsilon=1e-08, decay=0.005)
    adamax = optimizers.Adamax(lr=0.02, beta_1=0.9, beta_2=0.999, epsilon=1e-08, decay=0.001)
    adam = optimizers.Adam(lr=0.02, beta_1=0.9, beta_2=0.999, epsilon=1e-08, decay=0.001)
    
    model.compile(optimizer=adagrad, loss='binary_crossentropy', metrics=["binary_accuracy"])
    model.summary()
    
    print("Model fitting...")
    
    fBestModel = 'best_modelxdx5_' + str(iteration) + '.h5' 
    early_stop = EarlyStopping(monitor='val_loss', patience=150, verbose=1) 
    best_model = ModelCheckpoint(fBestModel, verbose=0, save_best_only=True)
    
    tensorboard_callback = TensorBoard(log_dir='./tb_logs_250430', histogram_freq=0,  
          write_graph=True, write_images=True)
    
    model.fit(x = [X_rnn_0[train], X_rnn_1[train], X_rnn_2[train], X_rnn_3[train], X_rnn_4[train], X_rnn_aai_x[train], X_ann[train]], 
              y = Y[train], 
              validation_data = (
                      [X_rnn_0[test], X_rnn_1[test], X_rnn_2[test], X_rnn_3[test], X_rnn_4[test],  X_rnn_aai_x[test], X_ann[test]], 
                      Y[test],
                      ), 
              epochs=500, batch_size=405, verbose=True,
              #shuffle=True, 
              callbacks=[
                      ##best_model, 
                      early_stop,
                      #tensorboard_callback,
                      ])
    print("Model predicting...")
    
    iteration = iteration + 1
