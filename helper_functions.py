# -*- coding: utf-8 -*-
"""
A set of helper functions imported by descriptors.py
"""
import numpy as np
from functools import lru_cache
aa_letters = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
di_letters = ["%s%s" % (a, b) for a in aa_letters for b in aa_letters]
tri_letters = ["%s%s%s" % (a, b, c) for a in aa_letters for b in aa_letters for c in aa_letters]
conjoint_letters = ["A", "I", "Y", "H", "R", "D", "C"]
veltri_letters =["H", "Y", "G", "A", "L", "Q", "D", "C"]
veltri_letters_di = ["%s%s" % (a, b) for a in veltri_letters for b in veltri_letters]
letters = {1 : aa_letters, 2 : di_letters, 3 : tri_letters, 4 : conjoint_letters}
from reduced_alphabets import conjoint_dict, TD_3, TD_5, TD_10, veltri

@lru_cache(maxsize=10)
def alpha2num_dict_f(alphabet=aa_letters):
    return {(["X"] + alphabet)[i] : i for i in range(0,len(alphabet)+1)}


def alpha2num(seq, alphabet=aa_letters, padding=200, just="right"):
    if just == "right":
        padded_seq = seq.rjust(padding, "X")
    elif just == "left":
        padded_seq = seq.ljust(padding, "X")
    return np.array([[alpha2num_dict_f(alphabet)[s] for s in padded_seq]])

def feature_occurence(peptides):
    '''
    A function that returns the total number of peptides in which a certain di-peptide or tri-peptide occurs.
    Returns two dicts eg. {"TW" : 48, "AQ" : 90} , {"GIG" : 20, "GAG" : 30,}
    '''
    dp = {i: 0 for i in letters[2]}
    tp = {i: 0 for i in letters[3]}

    for peptide in peptides:
        temp_set_2 = set() # Temp_set of di-amino-acids for each peptide, removes duplicates from each peptide
        temp_set_3 = set() # Temp_set of tri-amino-acids for each peptide, removes duplicates from each peptide
        seq = peptide["seq"]
        l = len(seq)
        
        for a in range(l-1):
            s = seq[a:a+2]
            temp_set_2.add(s)

        for s in temp_set_2:
            dp[s] = dp[s] + 1
       
        for a in range(l-2):
            s = seq[a:a+3]
            temp_set_3.add(s)
 
        for s in temp_set_3:
            tp[s] = tp[s] + 1
    
    return dp, tp

@lru_cache(maxsize=30)
def conjoint_seq(seq, alphabet = conjoint_dict):
    '''A function to return the conjoint sequence of a peptide sequence.
    This reduces the 20 amino acids to an alphabet of 7 amino acids'''
    
    #Conjoint src = https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0828-1
    
    return "".join([alphabet[letter] for letter in seq])

def unique_everseen(iterable, key=None):
    """
    Yield unique elements, preserving order.
        >>> list(unique_everseen('AAAABBBCCDAABBB'))
        ['A', 'B', 'C', 'D']
        >>> list(unique_everseen('ABBCcAD', str.lower))
        ['A', 'B', 'C', 'D']
    Sequences with a mix of hashable and unhashable items can be used.
    The function will be slower (i.e., `O(n^2)`) for unhashable items.
    """
    seenset = set()
    seenset_add = seenset.add
    seenlist = []
    seenlist_add = seenlist.append
    if key is None:
        for element in iterable:
            try:
                if element not in seenset:
                    seenset_add(element)
                    yield element
            except TypeError:
                if element not in seenlist:
                    seenlist_add(element)
                    yield element
    else:
        for element in iterable:
            k = key(element)
            try:
                if k not in seenset:
                    seenset_add(k)
                    yield element
            except TypeError:
                if k not in seenlist:
                    seenlist_add(k)
                    yield element
