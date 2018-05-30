# -*- coding: utf-8 -*-
from __future__ import print_function

from IPython import get_ipython

get_ipython().magic('reset -sf') 

from describe_sequences import describe_sequence
import pickle as pickle
import json
from modules.modlamp.descriptors import PeptideDescriptor, GlobalDescriptor
import numpy as np
from di2 import di2, di3
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from aaf import aaf, aaf_angular
from modules.PyBioMed.PyProtein import CTD, ConjointTriad
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
from property_pos_convx_angular import property_pos_convx_angular
from modules.pydpi.pypro import PyPro
random.seed(31) ## Set a random seed for reproducible results

path = "peptides_hemopi3.json"   #This file is produced by running read_files.py
with open(path, "r") as f:
    text = f.read()
    peptides = eval(text)["Peptides"]
    random.shuffle(peptides)

dp, tp = feature_occurence(peptides) ##A function that returns the number of peptides that each di or tri-peptide occurs in
                                         ##Useful to prevent overfitting later on.

pickle.dump( dp, open("dp.p", "wb"))
pickle.dump( tp, open("tp.p", "wb"))

for peptide in tqdm(peptides):
    peptide = describe_sequence(peptide)

for key in ["array", "alpha2num_veltri", "pos_first_aa_1", "pos_first_veltri", "pos_last_aa_1", "pos_last_veltri", "property_pos_first_1", "property_pos_last_1"]:
    x = np.concatenate([peptide[key] for peptide in peptides], axis=0)
    np.save("inputs/peptides_array_7_" + key, x, allow_pickle=False)
    
for identifier, mode in aai_to_get:
    x1 = np.concatenate([peptide[identifier+"_"+mode+"_first"] for peptide in peptides], axis=0)
    np.save("inputs/peptides_array_7_" + identifier+"_"+mode+"_first", x1, allow_pickle=False)
    x2 = np.concatenate([peptide[identifier+"_"+mode+"_last"] for peptide in peptides], axis=0)
    np.save("inputs/peptides_array_7_" + identifier+"_"+mode+"_last", x1, allow_pickle=False)

x = np.concatenate([peptide["property_pos_convx"] for peptide in peptides], axis=0)
print("MinMax Scaling...")
#x = minmax_scaler_3d(x, feature_range=(-1,1))
print("Standard Scaling...")
x = robust_scaler_3d(x, with_scaling=True, with_centering=True)
np.save("inputs/peptides_array_7_property_pos_convx", x, allow_pickle=False)
print("________________\n"*5)


x = np.concatenate([peptide["property_pos_convx_angular"] for peptide in peptides], axis=0)
print("MinMax Scaling...")
#x = minmax_scaler_3d(x, feature_range=(-1,1))
print("Standard Scaling...")
x = robust_scaler_3d(x, with_scaling=True, with_centering=True)
np.save("inputs/peptides_array_7_property_pos_convx_angular", x, allow_pickle=False)

print(x)    
    
print("complete")