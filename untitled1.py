# -*- coding: utf-8 -*-
"""
Created on Tue May  1 09:00:10 2018

@author: Patrick
"""

import numpy as np
from sklearn.preprocessing import StandardScaler, MinMaxScaler

a = np.array([[2], [4] ,[6]])

scaler = StandardScaler(with_mean=False, with_std=True)
a = scaler.fit_transform(a)

print(a)
