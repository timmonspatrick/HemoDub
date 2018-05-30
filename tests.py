# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 23:05:37 2018

@author: Patrick
"""
from __future__ import print_function

from helper_functions import conjoint_seq, feature_occurence, residue_distribution, residue_boolean, residue_abs
from positional import single_res_pos_first, dipeptide_res_pos_first, property_pos_first, single_res_pos_last, dipeptide_res_pos_last, property_pos_last
from aa_indeces import aai_to_get

path = "peptides.json"   #This file is produced by running read_files.py
with open(path, "r") as f:
    text = f.read()
    peptides = eval(text)["Peptides"]
    random.shuffle(peptides)

dp, tp = feature_occurence(peptides) ##A function that returns the number of peptides that each di or tri-peptide occurs in
print(tp)
freq_3d = residue_distribution("FFKFFK", 3, dp, tp)
print(freq_3d)

freq_1dbool = residue_boolean("CCCA", 1, dp)
print(freq_1dbool)

freq_2dbool = residue_boolean("ACCCA", 2, dp)
print(freq_2dbool)