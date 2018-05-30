# -*- coding: utf-8 -*-
"""
Created on Wed Dec 20 18:07:22 2017

@author: Patrick
"""
from __future__ import division
import aaindex
from aa_indeces import aai_to_get
from math import cos, sin, radians, degrees
aaindex.init(path='.')

def aaf(sequence, identifier, mode):
    x = aaindex.get(identifier)
    
    total = 0
    total_list = []
    for aa in sequence:
        total = total + x.get(aa) 
        total_list.append(x.get(aa))
        
    if mode == "mean":
        return total / len(sequence)
    elif mode == "total":
        return total
    elif mode == "max":
        return max(total_list)


def aaf_angular(sequence, identifier, mode):
    x = aaindex.get(identifier)
    
    total_cos = 0
    total_sin = 0
    total_list = []
    for aa_n in range(len(sequence)):
        
        aa = sequence[aa_n]
        
        angle = float(90/(len(sequence) - 1)) * aa_n
        
        total_cos = total_cos + ( x.get(aa) * cos(radians(angle)) )
        total_sin = total_sin + ( x.get(aa) * sin(radians(angle)) )
        total_list.append(x.get(aa) * cos(radians(angle)))
        

        
    if mode == "mean":
        try:
            return total_cos / (total_cos**2 + total_sin**2)**0.5
        except:
            return 0
    elif mode == "total":
        return total_cos
    elif mode == "max":
        return max(total_list)


aaf_angular("FLGALLGPLMNLLQ", "KLEP840101", "mean")

