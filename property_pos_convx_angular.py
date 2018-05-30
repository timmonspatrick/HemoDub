from __future__ import division, print_function
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 18:20:10 2018

@author: Patrick
"""
import numpy as np
from helper_functions import aa_letters, di_letters, unique_everseen
from aa_indeces import aai_to_get
import aaindex
from residue_distribution import alphabets
from math import cos, sin, degrees, radians

aaindex.init(path='.')

first_no = 10 # Number of amino acids to take from the beginning of the sequence
last_no = 10 # Number of amino acids to take from the end of the sequence

desired_identifiers = [i[0] for i in unique_everseen(aai_to_get)]

def property_pos_convx_angular(seq, f_length, number = 200):

    first_no = number
    len_seq = len(seq)
    seq = (seq + "X"*first_no)[:first_no]
    seq_list = [seq[a:a+f_length] for a in range(0, number)]
    
    l = []
    #print("WARNING: ONLY 7 AAIs selected. Remove before final.")
    for aa_n in range(len(seq_list)):
        l2 = []
        
        aa = seq[aa_n]
        angle = float(90/(len_seq - 1)) * aa_n
        
        if "X" not in aa:        
            for identifier in desired_identifiers:
                
                identifiers = []
                for letter in aa:
                    x = aaindex.get(identifier)
                    y = x.get(letter)
                    
                    z = y * cos(radians(angle))
                    identifiers.append(z)
                avg_val = sum(identifiers) / f_length
                
                l2.append(avg_val)
                
        else:
            for identifier in desired_identifiers:
                l2.append(0.9191919191919156)

        l.append(np.array(l2))

    l = np.array([l])
    
    return l

        

print("___________")
print(property_pos_convx_angular("ACEACE",1))