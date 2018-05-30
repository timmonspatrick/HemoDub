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

aaindex.init(path='.')

first_no = 10 # Number of amino acids to take from the beginning of the sequence
last_no = 10 # Number of amino acids to take from the end of the sequence

def res_pos_first(seq, seq_type, f_length, number=first_no):
    seq = (seq + "X"*number)[:number] # Pad short sequences with X letter, and cut all sequences to 30 aa's
    seq_list = [seq[a:a+f_length] for a in range(0, number)]
    r = np.array([[1 if s == a else 0 for s in seq_list for a in alphabets[seq_type][f_length]]])
    
    return r

def res_pos_last(seq, seq_type, f_length, number = last_no):
    seq = ("X"*number + seq)[-number:] # Pad short sequences with X letter, and cut all sequences to last_no aa's
    seq_list = [seq[a:a+f_length] for a in range(0, number)]
    r = np.array([[1 if s == a else 0 for s in seq_list for a in alphabets[seq_type][f_length]]])
    
    return r

desired_identifiers = [i[0] for i in unique_everseen(aai_to_get)]

def property_pos_first(seq, f_length, number = first_no):
    seq = (seq + "X"*first_no)[:first_no]
    seq_list = [seq[a:a+f_length] for a in range(0, number)]
    
    l = []
    for aa in seq_list:
        if "X" not in aa:        
            for identifier in desired_identifiers:
                identifiers = []
                for letter in aa:
                    x = aaindex.get(identifier)
                    y = x.get(letter)
                    identifiers.append(y)
                avg_val = sum(identifiers) / f_length
                
                l.append(avg_val)
                
        else:
            for identifier in desired_identifiers:
                l.append(0)
    l = np.array([l])  
    return l

def property_pos_convx(seq, f_length, number = 200):
    first_no = number
    seq = (seq + "X"*first_no)[:first_no]
    seq_list = [seq[a:a+f_length] for a in range(0, number)]
    
    l = []
    #print("WARNING: ONLY 7 AAIs selected. Remove before final.")
    for aa in seq_list:
        l2 = []
        if "X" not in aa:        
            for identifier in desired_identifiers:
                
                identifiers = []
                for letter in aa:
                    x = aaindex.get(identifier)
                    y = x.get(letter)
                    identifiers.append(y)
                avg_val = sum(identifiers) / f_length
                
                l2.append(avg_val)
                
        else:
            for identifier in desired_identifiers:
                l2.append(0.9191919191919156)

        l.append(np.array(l2))

    l = np.array([l])
    #print(l)
    return l

print("___________")
print(property_pos_convx("ACEACE",1))

def property_pos_first_conv(seq, identifier, f_length, number = first_no):
    seq = (seq + "X"*first_no)[:first_no]
    #seq = seq # No Padding, pad in Keras instead
    seq_list = [seq[a:a+f_length] for a in range(0, number)]
    
    l = []
    for aa in seq_list:
        if "X" not in aa:        
            identifiers = []
            for letter in aa:
                x = aaindex.get(identifier)
                y = x.get(letter)
                identifiers.append(y)
            avg_val = sum(identifiers) / f_length
            
            l.append(avg_val)
                
        else:            
            l.append(0)
    l = np.array([l])  
    return l

def property_pos_last(seq, f_length, number = last_no):
    seq = ("X"*last_no + seq)[-last_no:] # Pad short sequences with X letter, and cut all sequences to last_no aa's
    seq_list = [seq[a:a+f_length] for a in range(0, number)]
    
    l = []
    for aa in seq_list:
        if "X" not in aa:        
            for identifier in desired_identifiers:
                identifiers = []
                for letter in aa:
                    x = aaindex.get(identifier)
                    y = x.get(letter)
                    identifiers.append(y)
                avg_val = sum(identifiers) / f_length
                
                l.append(avg_val)
                
        else:
            for identifier in desired_identifiers:
                l.append(0)
    l = np.array([l])  
    return l


def property_pos_last_conv(seq, identifier, f_length, number = last_no):
    seq = ("X"*last_no + seq)[-last_no:] # Pad short sequences with X letter, and cut all sequences to last_no aa's
    seq_list = [seq[a:a+f_length] for a in range(0, number)]
    
    l = []
    for aa in seq_list:
        if "X" not in aa:        
            identifiers = []
            for letter in aa:
                x = aaindex.get(identifier)
                y = x.get(letter)
                identifiers.append(y)
            avg_val = sum(identifiers) / f_length
            
            l.append(avg_val)
                
        else:
            l.append(0)
    l = np.array([l])  
    return l

def property_di_first(seq, f_length = 2, freq = 3, number = first_no):
    seq = (seq + "X"*first_no)[:first_no]
    
    seq_list = []
    for a in range(len(seq)):
        try:
            x = "".join([seq[x] for x in range(a, a  + f_length * freq, freq)])
            seq_list.append(x)
        except IndexError:
            pass
    
    l = []
    for aa in seq_list:
        if "X" not in aa:        
            for identifier in desired_identifiers:
                identifiers = []
                for letter in aa:
                    x = aaindex.get(identifier)
                    y = x.get(letter)
                    identifiers.append(y)
                avg_val = sum(identifiers) / len(identifiers)
                
                l.append(avg_val)
                
        else:
            for identifier in desired_identifiers:
                l.append(0)
    l = np.array([l])  
    return l

def property_di_last(seq, f_length = 2, freq = 3, number = last_no):
    seq = ("X"*last_no + seq)[-last_no:] # Pad short sequences with X letter, and cut all sequences to last_no aa's
    
    seq_list = []
    for a in range(len(seq)):
        try:
            x = "".join([seq[x] for x in range(a, a  + f_length * freq, freq)])
            seq_list.append(x)
        except IndexError:
            pass
    
    l = []
    for aa in seq_list:
        if "X" not in aa:        
            for identifier in desired_identifiers:
                identifiers = []
                for letter in aa:
                    x = aaindex.get(identifier)
                    y = x.get(letter)
                    identifiers.append(y)
                avg_val = sum(identifiers) / len(identifiers)
                
                l.append(avg_val)
                
        else:
            for identifier in desired_identifiers:
                l.append(0)
    l = np.array([l])  
    return l
