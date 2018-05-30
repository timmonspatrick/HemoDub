# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 12:20:31 2018

@author: Patrick
"""
import numpy as np

def filter_array_by_column(X, cut_off=0.95):
    n_cols = X.shape[1]
    bad_cols = set()
    for n in range(n_cols):
        Y = X[:,n]
        unique, counts = np.unique(Y, return_counts=True)
        
        counts_sum = sum(counts)
        counts = [i / counts_sum for i in counts]
        
        if len([i for i in counts if i >= cut_off]) > 0:
            bad_cols.add(n)

    good_cols = [i for i in range(n_cols) if i not in bad_cols]
    
    X_new = X[:,good_cols]
    return X_new