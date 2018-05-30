# -*- coding: utf-8 -*-
#
 
from __future__ import division

"""
Created on Sun Mar 04 15:30:02 2018
Hemolytic potency descriptor
@author: Patrick
"""

from helper_functions import aa_letters
import numpy as np

hemopotent_difs = {'A': -20, 'C': 17, 'E': -16, 'D': -7, 
                   'G': 17, 'F': 19, 'I': 24, 'H': -24, 
                   'K': -15, 'M': -7, 'L': 118, 'N': -6, 
                   'Q': -18, 'P': 7, 'S': 6, 'R': -78, 
                   'T': 15, 'W': 8, 'V': -20, 'Y': -18}

hemolytic = {'A': 155, 'C': 94, 'E': 19, 'D': 23, 
                  'G': 150, 'F': 136, 'I': 138, 'H': 27, 
                  'K': 378, 'M': 14, 'L': 381, 'N': 36, 
                  'Q': 20, 'P': 62, 'S': 94, 'R': 130, 
                  'T': 57, 'W': 100, 'V': 112, 'Y': 30}

nonhemolytic = {'A': 175, 'C': 77, 'E': 35, 'D': 30, 
             'G': 133, 'F': 117, 'I': 114, 'H': 51, 
             'K': 393, 'M': 21, 'L': 263, 'N': 42, 
             'Q': 38, 'P': 55, 'S': 88, 'R': 208, 
             'T': 42, 'W': 92, 'V': 132, 'Y': 48}

hemopotent_dict = {a: (hemolytic[a] - nonhemolytic[a] ) / ((hemolytic[a] + nonhemolytic[a])/2) for a in aa_letters}

print(hemopotent_dict)
def hemopotent(seq):
    hemolytic_abs = sum([hemopotent_dict[a] for a in seq])
    hemolytic_avg = hemolytic_abs / len(seq)
    h_r = np.array([[hemolytic_avg]])
    return h_r

