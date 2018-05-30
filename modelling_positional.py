# -*- coding: utf-8 -*-
"""
Takes as input peptides_array.npy

achieves 80% with HemoPi2
"""
from IPython import get_ipython
get_ipython().magic('reset -sf') 
import numpy as np
from keras.models import Sequential
from keras.layers import Dense, Activation, Dropout, Conv1D, Embedding, MaxPooling1D, LSTM, Merge, Bidirectional, GRU, SimpleRNN, GaussianNoise, Flatten
from sklearn.model_selection import train_test_split
from keras.callbacks import EarlyStopping, ModelCheckpoint, TensorBoard
from keras import optimizers
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.model_selection import StratifiedKFold
from aa_indeces import aai_to_get
from scalers_3d import standard_scaler_3d, minmax_scaler_3d

from tensorflow.python.client import device_lib

print(device_lib.list_local_devices())



inp = np.load("inputs/peptides_array_7_array.npy", allow_pickle=False)

X = np.array([i for i in inp])



Y = np.array([i[-1:] for i in inp])

print(Y)

kfold = StratifiedKFold(n_splits=5, shuffle=True, random_state=441653)
cvsscores = []

####Scaling#####

scaler = MinMaxScaler(feature_range=(-1, 1))
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


X_rnn_aais = [np.load("inputs/peptides_array_7_" + identifier+"_"+ mode+"_" + which+".npy") for identifier, mode in aai_to_get for which in ["first","last"]]
X_rnn_aai_x = np.load("inputs/peptides_array_7_property_pos_convx.npy")

X_rnn_aai_shape = X_rnn_aai_x.shape
print(X_rnn_aai_shape)
#X = np.concatenate((X_rnn_0, X_rnn_1, X_rnn_2, X_rnn_3, X_rnn_4, X_rnn_5, X_rnn_6, X_ann), axis=1)

if False:
    for i in range(X_rnn_aai_x.shape[2]):
        scaler = MinMaxScaler(feature_range=(-1, 1))
        X_rnn_aai_x[:,:,i] = scaler.fit_transform(X_rnn_aai_x[:,:,i])
        
        scaler = StandardScaler(with_mean=True, with_std=True)
        scaler.fit(X_rnn_aai_x[:,:,i])
        X_rnn_aai_x[:,:,i] = scaler.transform(X_rnn_aai_x[:,:,i])
 
#SCALING NOW DONE BY descriptors7.py       
#print("MinMax Scaling...")
#X_rnn_aai_x = minmax_scaler_3d(X_rnn_aai_x, feature_range=(-1,1))
#print("Standard Scaling...")
#X_rnn_aai_x = standard_scaler_3d(X_rnn_aai_x, with_std=True, with_mean=True)

dims = X.shape[1]

print(dims, 'dims')
print("Building model...")

nb_classes = Y.shape[1]
print(nb_classes, 'classes')

rate = 0.1

iteration = 0
for train, test in kfold.split(X_rnn_aai_x, Y):
    
    print(X_rnn_aai_x[train])
    print(X_rnn_aai_x[train].ndim)
    
    model_aai = Sequential()
    model_aai.add(Conv1D(128, 3, strides=1,  init="normal", activation="tanh", input_shape=(X_rnn_aai_shape[1], X_rnn_aai_shape[2])))  # X_rnn_aai_shape[1:]
    model_aai.add(MaxPooling1D(pool_size=2, strides=2))
    #model_aai.add(Dropout(0.25))
   
    model_aai.add(Conv1D(128, 3, strides=1, init="normal", activation="tanh"))  # X_rnn_aai_shape[1:]
    model_aai.add(MaxPooling1D(pool_size=2, strides=2))
    #model_aai.add(Dropout(0.25))
   
    model_aai.add(Conv1D(128, 3, strides=1, init="normal", activation="tanh"))  # X_rnn_aai_shape[1:]
    model_aai.add(MaxPooling1D(pool_size=2, strides=2))
    #model_aai.add(Dropout(0.2))
    
    model_aai.add(Conv1D(128, 3, strides=1, init="normal", activation="tanh"))  # X_rnn_aai_shape[1:]
    model_aai.add(MaxPooling1D(pool_size=2, strides=2))
    #model_aai.add(Dropout(0.25))

    model_aai.add(Conv1D(128, 3, strides=1, init="normal", activation="tanh"))  # X_rnn_aai_shape[1:]
    model_aai.add(MaxPooling1D(pool_size=2, strides=2))
    #model_aai.add(Dropout(0.25))

    #model_aai.add(GaussianNoise(0.05))
    #model_aai.add(Bidirectional(GRU(32, unroll = True, stateful = False, dropout = 0.88)))
    model_aai.add(Flatten())
    model_aai.add(Dropout(0.933))
    model_aai.add(Dense(1, init='normal', activation="sigmoid"))
    
    adam = optimizers.Adam(lr=0.02, beta_1=0.9, beta_2=0.999, epsilon=1e-08, decay=0.001)
    
    model_aai.compile(optimizer=optimizers.Adagrad(lr=0.01, decay=0.005), loss='binary_crossentropy', metrics=["binary_accuracy"])
    model_aai.summary()
    
    print("Model fitting...")
    
    fBestModel = 'best_modelxdx5_' + str(iteration) + '.h5' 
    early_stop = EarlyStopping(monitor='val_loss', patience=100, verbose=1) 
    best_model = ModelCheckpoint(fBestModel, verbose=0, save_best_only=True)
    
    tensorboard_callback = TensorBoard(log_dir='./tb_logs_010529', histogram_freq=0,  
          write_graph=True, write_images=True)
    
    
    
    model_aai.fit(x = X_rnn_aai_x[train], 
              y = Y[train], 
              validation_data = (
                      X_rnn_aai_x[test], 
                      Y[test],
                      ), 
              epochs=500, batch_size=32, verbose=True,
              #shuffle=True, 
              callbacks=[
                      ##best_model, 
                      early_stop,
                      tensorboard_callback,
                      ])
    print("Model predicting...")
    score = model_aai.evaluate(X_rnn_aai_x[test], Y[test], batch_size=256, verbose=1)
    cvsscores.append(score[1] * 100)
    print("Iteration No. " + str(iteration))
    print("Evaluation Score: "+ str(score))
    iteration = iteration + 1

print(cvsscores)