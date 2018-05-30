# -*- coding: utf-8 -*-
"""
Takes as input peptides_array.npy

achieves 80% with HemoPi2
"""
from IPython import get_ipython
get_ipython().magic('reset -sf') 
import numpy as np
from keras.models import Sequential
from keras.layers import Dense, Activation, Dropout, Conv1D, Embedding, MaxPooling1D, LSTM, Merge, GRU, Bidirectional, Flatten, GaussianNoise
from sklearn.model_selection import train_test_split
from keras.callbacks import EarlyStopping, ModelCheckpoint, TensorBoard
from keras import optimizers
from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler
from sklearn.decomposition import PCA
from sklearn.model_selection import StratifiedKFold
from aa_indeces import aai_to_get

from tensorflow.python.client import device_lib

print(device_lib.list_local_devices())



inp = np.load("inputs/peptides_array_7_array.npy", allow_pickle=False)

X_ann = np.array([i[:-1] for i in inp])

Y = np.array([i[-1:] for i in inp])

kfold = StratifiedKFold(n_splits=5, shuffle=True, random_state=441653 +1)
cvsscores = []

####Scaling#####

#scaler = RobustScaler(with_centering=True, with_scaling=True)
scaler = MinMaxScaler(feature_range=(0,1))
scaler.fit(X_ann)
X_ann = scaler.transform(X_ann)
print(X_ann.shape)
pca = PCA(n_components=0.98, svd_solver="full")
X_ann = pca.fit_transform(X_ann)
print(X_ann.shape)
#######

dims = X_ann.shape[1]
print(dims, 'dims')
print("Building model...")

nb_classes = Y.shape[1]
print(nb_classes, 'classes')

rate = 0.1

iteration = 0

dropout=0.89

cvsscores = []

    
dims_ann = X_ann.shape[1]

models = []

model_ann = Sequential()
model_ann.add(Dense(512, input_shape=(dims,), init='uniform', activation="relu"))
model_ann.add(Dropout(dropout, noise_shape=None, seed=42))
model_ann.add(Dense(512,  init='uniform', activation='relu'))
model_ann.add(Dropout(dropout, noise_shape=None, seed=42))
model_ann.add(Dense(512,  init='uniform', activation='relu'))
model_ann.add(Dropout(dropout, seed=42, noise_shape=None))

model_ann.add(Dense(1, activation="sigmoid"))


sgd = optimizers.SGD(lr=0.01, decay=0.005, momentum=0.5, nesterov=True)
rmsprop = optimizers.RMSprop(lr=0.01, rho=0.9, epsilon=1e-08, decay=0.005)
#adagrad = optimizers.Adagrad(lr=0.005)
adagrad = optimizers.Adagrad(lr=0.01, decay=0.005)
adadelta = optimizers.Adadelta(lr=0.01, rho=0.95, epsilon=1e-08, decay=0.005)
adamax = optimizers.Adamax(lr=0.01, beta_1=0.9, beta_2=0.999, epsilon=1e-08, decay=0.005)
adam = optimizers.Adam(lr=0.01, beta_1=0.9, beta_2=0.999, epsilon=1e-08, decay=0.005)

model_ann.compile(optimizer=adagrad, loss='binary_crossentropy', metrics=["binary_accuracy"])
model_ann.summary()

print("Model fitting...")

fBestModel = 'best_model_prod.h5' 

early_stop = EarlyStopping(monitor='val_loss', patience=100, verbose=1) 
best_model = ModelCheckpoint(fBestModel, verbose=0, save_best_only=True)

tensorboard_callback = TensorBoard(log_dir='./tb_logs_020502', histogram_freq=0,  
      write_graph=True, write_images=True)

model_ann.fit(x = X_ann, 
          y = Y, 
          validation_data = (
                  X_ann, 
                  Y,
                  ), 
          epochs=1000, batch_size=202, verbose=True,
          #shuffle=True, 
          callbacks=[
                  best_model, 
                  early_stop,
                  tensorboard_callback,
                  ])
print("Model predicting...")
cvsscores.append(model_ann.evaluate(x=X_ann[test], y=Y[test]))
iteration = iteration + 1

for i in cvsscores:
    print(i)