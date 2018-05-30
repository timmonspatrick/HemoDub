# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 16:33:36 2018

@author: Patrick
"""

from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler
import numpy as np
from tqdm import tqdm

def minmax_scaler_3d(X, feature_range=(-1,1)):
    for i in tqdm(range(X.shape[2])):
        scaler = MinMaxScaler(feature_range=feature_range)
        
        X_flat = np.array([ [x] for x in X[:,:,i].flatten() if x != 0.9191919191919156 ] )
        #print(X_flat)
        #print("_________")
        scaler.fit(X_flat)
        for a in range(X.shape[0]):
            for b in range(X.shape[1]):
                if X[a,b,i] != 0.9191919191919156:
                    X[a,b,i] = scaler.transform(X[a,b,i])
                else:
                    X[a,b,i] = 0.9191919191919156
    return X

#X = np.random.rand(2, 4, 6)
#print(X)
#print("_"*40)
#print(minmax_scaler_3d(X))


def standard_scaler_3d(X, with_mean=True, with_std=True):
    for i in tqdm(range(X.shape[2])):
        scaler = StandardScaler(with_mean=with_mean, with_std=with_std)
        
        X_flat = np.array([ [x] for x in X[:,:,i].flatten() if x != 0.9191919191919156 ] )
        #print(X_flat)
        #print("_________")
        scaler.fit(X_flat)
        for a in range(X.shape[0]):
            for b in range(X.shape[1]):
                if X[a,b,i] != 0.9191919191919156:
                    X[a,b,i] = scaler.transform(X[a,b,i])
                else:
                    X[a,b,i] = scaler.mean_[0]
    return X

def robust_scaler_3d(X, with_centering=True, with_scaling=True, quantile_range=(25.0, 75.0), copy=True):
    for i in tqdm(range(X.shape[2])):
        scaler = RobustScaler(with_centering=with_centering, with_scaling=with_scaling, quantile_range=quantile_range, copy=copy)
        
        X_flat = np.array([ [x] for x in X[:,:,i].flatten() if x != 0.9191919191919156 ] )
        #print(X_flat)
        #print("_________")
        scaler.fit(X_flat)
        for a in range(X.shape[0]):
            for b in range(X.shape[1]):
                if X[a,b,i] != 0.9191919191919156:
                    X[a,b,i] = scaler.transform(X[a,b,i])
                else:
                    X[a,b,i] = scaler.center_[0]
    return X