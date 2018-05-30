# -*- coding: utf-8 -*-
"""
Takes as input peptides_array.npy

"""
from IPython import get_ipython
get_ipython().magic('reset -sf') 
import numpy as np
from keras.models import Sequential
from keras.layers import Dense, Activation, Dropout, Conv1D, Embedding
from sklearn.model_selection import train_test_split
from keras.callbacks import EarlyStopping, ModelCheckpoint, TensorBoard
from keras import optimizers
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.model_selection import StratifiedKFold

from tensorflow.python.client import device_lib

print(device_lib.list_local_devices())

inp = np.load("peptides_array_cdhit.npy", allow_pickle=False)

X = np.array([i[201:-1] for i in inp])
Y = np.array([i[-1:] for i in inp])
print(Y)
kfold = StratifiedKFold(n_splits=5, shuffle=True, random_state=441653)
cvsscores = []

####Scaling#####
print(X)
scaler = MinMaxScaler(feature_range=(-1, 1))
X_1 = scaler.fit_transform(X)

scaler = StandardScaler(with_mean=True, with_std=True)
scaler.fit(X_1)
X = scaler.transform(X_1)

print(X)
#######


dims = X.shape[1]

print(dims, 'dims')
print("Building model...")

nb_classes = Y.shape[1]
print(nb_classes, 'classes')

rate = 0.1

iteration = 0
for train, test in kfold.split(X, Y):

    model = Sequential()
    #model.add(Dense(1024, input_shape=(dims,), init='uniform', activation='relu'))
    model.add(Dropout(0.85, input_shape=(dims,), noise_shape=None, seed=42))
    model.add(Dense(64,  init='normal', activation='relu'))
    model.add(Dense(32,  init='normal', activation='relu'))
    model.add(Dropout(0.5, noise_shape=None, seed=42))
    model.add(Dense(1,  init='normal', activation='sigmoid'))
    
    
    sgd = optimizers.SGD(lr=0.01, decay=0.005, momentum=0.5, nesterov=True)
    rmsprop = optimizers.RMSprop(lr=0.02, rho=0.9, epsilon=1e-08, decay=0.01)
    #adagrad = optimizers.Adagrad(lr=0.005)
    adagrad = optimizers.Adagrad(lr=0.02)
    adadelta = optimizers.Adadelta(lr=0.01, rho=0.95, epsilon=1e-08, decay=0.005)
    adamax = optimizers.Adamax(lr=0.02, beta_1=0.9, beta_2=0.999, epsilon=1e-08, decay=0.001)
    adam = optimizers.Adam(lr=0.02, beta_1=0.9, beta_2=0.999, epsilon=1e-08, decay=0.001)
    
    model.compile(optimizer=adagrad, loss='binary_crossentropy', metrics=["binary_accuracy"])
    model.summary()
    
    print("Model fitting...")
    
    fBestModel = 'best_modelxdx5_' + str(iteration) + '.h5' 
    early_stop = EarlyStopping(monitor='val_loss', patience=150, verbose=1) 
    best_model = ModelCheckpoint(fBestModel, verbose=0, save_best_only=True)
    
    tensorboard_callback = TensorBoard(log_dir='./tb_logs_030430', histogram_freq=0,  
          write_graph=True, write_images=True)
    
    model.fit(X[train], Y[train], validation_data = (X[test], Y[test]), 
              epochs=10000, batch_size=dims, verbose=True,
              #shuffle=True, 
              callbacks=[
                      #best_model, 
                      early_stop,
                      tensorboard_callback,
                      ])
    print("Model predicting...")
    score = model.evaluate(X[test], Y[test], batch_size=32, verbose=1)
    cvsscores.append(score[1] * 100)
    print("Iteration No. " + str(iteration))
    print("Evaluation Score: "+ str(score))
    iteration = iteration + 1
