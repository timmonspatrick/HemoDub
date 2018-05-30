# -*- coding: utf-8 -*-
"""
Created on Thu May 03 21:07:32 2018

@author: Patrick
"""
import aaindex
from aaindex import grep, search
aaindex.init(path='.')
for r in search(""):
    print(r)