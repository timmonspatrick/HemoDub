# -*- coding: utf-8 -*-
"""
Created on Thu Mar 01 19:39:53 2018

@author: Patrick
"""
import numpy as np
from functools import lru_cache

aa_letters_1 = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
aa_letters_2 = ["%s%s" % (a, b) for a in aa_letters_1 for b in aa_letters_1]
aa_letters_3 = ["%s%s%s" % (a, b, c) for a in aa_letters_1 for b in aa_letters_1 for c in aa_letters_1]
conjoint_letters_1 = ["A", "I", "Y", "H", "R", "D", "C"]
conjoint_letters_2 = ["%s%s" % (a, b) for a in conjoint_letters_1 for b in conjoint_letters_1]
conjoint_letters_3 = ["%s%s%s" % (a, b, c) for a in conjoint_letters_1 for b in conjoint_letters_1 for c in conjoint_letters_1]
TD_3_letters_1 = ["V", "C", "E"]
TD_3_letters_2 = ["%s%s" % (a, b) for a in TD_3_letters_1 for b in TD_3_letters_1]
TD_3_letters_3 = ["%s%s%s" % (a, b, c) for a in TD_3_letters_1 for b in TD_3_letters_1 for c in TD_3_letters_1]
TD_5_letters_1 = ["V", "C", "E", "R", "G"]
TD_5_letters_2 = ["%s%s" % (a, b) for a in TD_5_letters_1 for b in TD_5_letters_1]
TD_5_letters_3 = ["%s%s%s" % (a, b, c) for a in TD_5_letters_1 for b in TD_5_letters_1 for c in TD_3_letters_1]
TD_10_letters_1 = ['A', 'C', 'E', 'G', 'H', 'P', 'S', 'R', 'T', 'W', 'V']
TD_10_letters_2 = ["%s%s" % (a, b) for a in TD_10_letters_1 for b in TD_10_letters_1]
TD_10_letters_3 = ["%s%s%s" % (a, b, c) for a in TD_10_letters_1 for b in TD_10_letters_1 for c in TD_10_letters_1]
veltri_letters_1 = ['H', 'Y', 'G', 'A', 'L', 'Q', 'D', 'C']
veltri_letters_2 = ["%s%s" % (a, b) for a in veltri_letters_1 for b in veltri_letters_1]
veltri_letters_3 = ["%s%s%s" % (a, b, c) for a in veltri_letters_1 for b in veltri_letters_1 for c in veltri_letters_1]


alphabets = {"aa" : {1 : aa_letters_1, 2 : aa_letters_2, 3 : aa_letters_3},
             "conjoint" : {1 : conjoint_letters_1, 2: conjoint_letters_2, 3 : conjoint_letters_3},
             "TD_3" : {1 : TD_3_letters_1, 2 : TD_3_letters_2, 3 : TD_3_letters_3},
             "TD_5" : {1 : TD_5_letters_1, 2 : TD_5_letters_2, 3 : TD_5_letters_3},
             "TD_10" : {1 : TD_10_letters_1, 2 : TD_10_letters_2, 3 : TD_10_letters_3},
             "veltri" : {1 : veltri_letters_1, 2 : veltri_letters_2, 3 : veltri_letters_3}
             }

@lru_cache(maxsize=30)
def counter(sequence, seq_type, f_length):
    l = len(sequence)
    d = {i : 0 for i in alphabets[seq_type][f_length]}
    
    for a in range(l-f_length+1):
        s = sequence[a:a+f_length]
        try:
            d[s] += 1.0
        except KeyError:
            d[s] = 1.0
    
    return d

def residue_distribution(sequence, seq_type, f_length, typ, constraint=None, constraint_val = None):
    '''an improved residue distribution function
    typ = 'distribution', 'absolute', 'boolean'
    '''
    
    d = counter(sequence, seq_type, f_length)
    if constraint:
        residue_counts = list(sorted([ (i, d[i]) for i in alphabets[seq_type][f_length] if constraint[i] >= constraint_val]))
    else:
        residue_counts = list(sorted([(i, d[i]) for i in alphabets[seq_type][f_length]]))
    
    if typ == "distribution":
        residue_counts = [(i[0], i[1] / ( len(sequence) - f_length + 1) ) for i in residue_counts]

    elif typ == "boolean":
        residue_counts = [(i[0], 1 if i[1] > 0 else 0) for i in residue_counts]
    elif typ == "absolute":
        residue_counts = residue_counts
    
    r_c = [i[1] for i in residue_counts]
    dis = np.array([r_c,])
    return dis