# -*- coding: utf-8 -*-

from __future__ import print_function
from IPython import get_ipython
get_ipython().magic('reset -sf') 
import json
from modlamp.descriptors import PeptideDescriptor, GlobalDescriptor
import numpy as np
from di2 import di2, di3
from pydpi.pypro import PyPro
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from aaf import aaf
from PyBioMed.PyProtein import CTD, ConjointTriad
import random
from helper_functions import conjoint_seq, feature_occurence, alpha2num, veltri_letters, aa_letters, conjoint_letters
from residue_distribution import residue_distribution
from positional import res_pos_first, property_pos_first, res_pos_last, property_pos_last, property_di_first, property_di_last, property_pos_first_conv, property_pos_last_conv, property_pos_convx
from aa_indeces import aai_to_get
from pep_descriptor import pep_descriptor
from tqdm import tqdm
from hemopotent import hemopotent
from reduced_alphabets import TD_3, TD_5, TD_10, veltri, conjoint_dict
from scalers_3d import standard_scaler_3d, minmax_scaler_3d, robust_scaler_3d
from D_descriptors import D_descriptor
from math import cos, radians
import aaindex
random.seed(31) ## Set a random seed for reproducible results

path = "peptides_hemopi3.json"   #This file is produced by running read_files.py
with open(path, "r") as f:
    text = f.read()
    peptides = eval(text)["Peptides"]
    random.shuffle(peptides)

dp, tp = feature_occurence(peptides) ##A function that returns the number of peptides that each di or tri-peptide occurs in
                                         ##Useful to prevent overfitting later on.

def val(angle):
    return max(0, cos(radians(angle)) + 0.5)

def wheel(seq, identifier):
    all_totals = []
    for starting_degree in range(0, 360, 10):
        total = 0
        angle = starting_degree
        for aa in seq:
            x = aaindex.get(identifier)
            value = x.get(aa)
            value = value * val(angle)
            total = total + value
            angle = angle + 100
        all_totals.append(total)
    
    return max(all_totals) / len(seq)


def describe_sequence(peptide):
    #print(peptide["seq"])
     ###Calculate all the Global descriptors##################################
   
    
    if peptide["activity"] == "YES":
        pepact = 1
    else:
        pepact = 0
    pepact = np.array([[pepact]])
    
    
    ############################################################################################################
    peptide["array"] = np.concatenate((
           np.array([[wheel(peptide["seq"], identifier) for identifier, mode in aai_to_get]]),   #aminoacidindeces,

           pepact,
           ), axis=1)

    
        
    return peptide

for peptide in tqdm(peptides):
    peptide = describe_sequence(peptide)


x = np.concatenate([peptide["array"] for peptide in peptides], axis=0)
np.save("inputs/peptides_array_aahydrophobicity", x, allow_pickle=False)
print("________________\n"*5)
print(x)    
    
print("complete")