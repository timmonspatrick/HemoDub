# -*- coding: utf-8 -*-
"""
Created on Wed May 30 20:56:29 2018

@author: Patrick
"""

from __future__ import print_function 
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
import pickle
random.seed(31) ## Set a random seed for reproducible results

dp = pickle.load( open( "dp.p", "rb"))
tp = pickle.load( open( "tp.p", "rb"))

def describe_sequence(peptide):
    #print(peptide["seq"])
     ###Calculate all the Global descriptors##################################
    globdesc = GlobalDescriptor(peptide["seq"])
    globdesc.calculate_charge(ph=7.4, amide=False, append=False)
    #print(globdesc.descriptor)
    
    globdesc.calculate_all(amide = peptide["cTer"] == "Amidation")
    ## Calculate all the peptide descriptors############################################################
    pepdesc = pep_descriptor(peptide["seq"])
    #########################################################################
    
    ##Composition, Transition and Distribution descriptors based on the different properties of AADs####
    ctdc = CTD.CalculateC(peptide["seq"])
    ctdc_keys = list(sorted(list([key for key in ctdc])))
    ctdc_vals = np.array([[ctdc[key] for key in ctdc_keys]]) 
    ####################################################################################################
    
    ##Descriptors based on the Conjoint alphabet########################################################
    
    conjointtriad = ConjointTriad.CalculateConjointTriad(peptide["seq"])
    conjointtriad_keys = list(sorted(list([key for key in conjointtriad])))
    conjointtriad_vals = np.array([[conjointtriad[key] for key in conjointtriad_keys]]) 
    #Count of conjoint triples####    
    
    
    ### Pseudo amino acid composition descriptors##########################################################
    protein = PyPro()
    protein.ReadProteinSequence(peptide["seq"])
    paac=protein.GetPAAC(lamda=1, weight=0.05)
    paac2 = [[
            paac[a] for a in list(sorted([k for k in paac], key = lambda x : int(x.replace("PAAC",""))))
            ]]        
    paac = np.array(paac2)
    #######################################################################################################
    
    cTer = np.array([[1 if peptide["cTer"] == "Amidation" else 0]])
    
    
    ### List of the fraction of amino acids which tend to be in helix, turn or sheet#######################
    analysed_seq = ProteinAnalysis(peptide["seq"])
    secondary_structure_fraction = np.array([analysed_seq.secondary_structure_fraction()])
    #######################################################################################################
    
    try:
        pepid = np.array([[int(peptide["id"].replace("HEMOLYTIK","99909").replace("DRAMP","99919").replace("DBAASP","99929"))]])
    except KeyError:
        pepid = 0
    
    len_peptide = np.array([[len(peptide["seq"])]])
    
    if peptide["activity"] == "YES":
        pepact = 1
    else:
        pepact = 0
    pepact = np.array([[pepact]])
    
    ########################################################################################################
    
    ## Tetra-peptide analysis
    ## Binary descriptor describing whether these tetra-peptides are present
    ## Previous paper has found these to be contributory to being hemolytic
    tetra_peptides = ["KLLL",    # src = https://github.com/riteshcanfly/Hemopi/blob/master/tetrapos.txt
                      "GCSC",
                      "AAAK",
                      "KLLS",
                      "LGKL",
                      "VLKA",
                      "LLGK",
                      "LVGA",
                      "LSDF",
                      "SDFK",
                      "SWLR",
                      "WLRD",]
    
    tp_bin = [] 
    for t_p in tetra_peptides:
        if t_p in peptide["seq"]:
            tp_bin.append(1)
        else:
            tp_bin.append(0)
    tp_bin = np.array([tp_bin])


    ############################################################################################################
    peptide["array"] = np.concatenate((
           pepdesc.descriptor, globdesc.descriptor, len_peptide, 
           cTer,
           secondary_structure_fraction,
           np.array([[aaf(peptide["seq"], identifier, mode) for identifier, mode in aai_to_get]]),   #aminoacidindeces,
           np.array([[aaf_angular(peptide["seq"], identifier, mode) for identifier, mode in aai_to_get]]),   #aminoacidindeces,
           ctdc_vals,
           tp_bin, 
           
           np.array([[D_descriptor(peptide["seq"])[0]]]), 
           np.array([[D_descriptor(peptide["seq"])[1]]]), 
           np.array([[D_descriptor(peptide["seq"])[2]]]), 
           np.array([[D_descriptor(peptide["seq"])[3]]]), 
           np.array([[D_descriptor(peptide["seq"])[4]]]), 
           np.array([[D_descriptor(peptide["seq"])[5]]]), 
           np.array([[D_descriptor(peptide["seq"])[6]]]), 

           residue_distribution(peptide["seq"], "aa", 1, "distribution"),  #freq_1d, 
           residue_distribution(peptide["seq"], "aa", 2, "distribution", constraint=dp, constraint_val = 50),  #freq_2d,
           #residue_distribution(peptide["seq"], "aa", 3, "distribution", constraint=tp, constraint_val = 20),
           residue_distribution(peptide["seq"], "aa", 1, "boolean"),  #freq_1d, 
           #residue_distribution(peptide["seq"], "aa", 2, "boolean", constraint=dp, constraint_val = 50),  #freq_2d,
           #residue_distribution(peptide["seq"], "aa", 3, "boolean", constraint=tp, constraint_val = 20),
           residue_distribution(peptide["seq"], "aa", 1, "absolute"),  #freq_1d, 
           #residue_distribution(peptide["seq"], "aa", 2, "absolute", constraint=dp, constraint_val = 50),  #freq_2d,
           #residue_distribution(peptide["seq"], "aa", 3, "absolute", constraint=tp, constraint_val = 20),
           
           residue_distribution(peptide["conjoint_seq"], "conjoint", 1, "distribution"),  #freq_1d, 
           residue_distribution(peptide["conjoint_seq"], "conjoint", 2, "distribution", ),  #freq_2d,
           #residue_distribution(peptide["conjoint_seq"], "conjoint", 3, "distribution", ),
           residue_distribution(peptide["conjoint_seq"], "conjoint", 1, "boolean"),  #freq_1d, 
           residue_distribution(peptide["conjoint_seq"], "conjoint", 2, "boolean", ),  #freq_2d,
           #residue_distribution(peptide["conjoint_seq"], "conjoint", 3, "boolean",),
           residue_distribution(peptide["conjoint_seq"], "conjoint", 1, "absolute"),  #freq_1d, 
           residue_distribution(peptide["conjoint_seq"], "conjoint", 2, "absolute", ),  #freq_2d,
           #residue_distribution(peptide["conjoint_seq"], "conjoint", 3, "absolute", ),
           
           #residue_distribution(peptide["TD_3_seq"], "TD_3", 1, "distribution"),  #freq_1d, 
           #residue_distribution(peptide["TD_3_seq"], "TD_3", 2, "distribution", ),  #freq_2d,
           #residue_distribution(peptide["TD_3_seq"], "TD_3", 3, "distribution", ),
           #residue_distribution(peptide["TD_3_seq"], "TD_3", 1, "boolean"),  #freq_1d, 
           #residue_distribution(peptide["TD_3_seq"], "TD_3", 2, "boolean", ),  #freq_2d,
           #residue_distribution(peptide["TD_3_seq"], "TD_3", 3, "boolean",),
           #residue_distribution(peptide["TD_3_seq"], "TD_3", 1, "absolute"),  #freq_1d, 
           #residue_distribution(peptide["TD_3_seq"], "TD_3", 2, "absolute", ),  #freq_2d,
           #residue_distribution(peptide["TD_3_seq"], "TD_3", 3, "absolute", ),
           
           residue_distribution(peptide["TD_5_seq"], "TD_5", 1, "distribution"),  #freq_1d, 
           residue_distribution(peptide["TD_5_seq"], "TD_5", 2, "distribution", ),  #freq_2d,
           #residue_distribution(peptide["TD_5_seq"], "TD_5", 3, "distribution", ),
           #residue_distribution(peptide["TD_5_seq"], "TD_5", 1, "boolean"),  #freq_1d, 
           #residue_distribution(peptide["TD_5_seq"], "TD_5", 2, "boolean", ),  #freq_2d,
           #residue_distribution(peptide["TD_5_seq"], "TD_5", 3, "boolean",),
           #residue_distribution(peptide["TD_5_seq"], "TD_5", 1, "absolute"),  #freq_1d, 
           #residue_distribution(peptide["TD_5_seq"], "TD_5", 2, "absolute", ),  #freq_2d,
           #residue_distribution(peptide["TD_5_seq"], "TD_5", 3, "absolute", ),
           
           #residue_distribution(peptide["TD_10_seq"], "TD_10", 1, "distribution"),  #freq_1d, 
           #residue_distribution(peptide["TD_10_seq"], "TD_10", 2, "distribution", ),  #freq_2d,
           #residue_distribution(peptide["TD_10_seq"], "TD_10", 3, "distribution", ),
           #residue_distribution(peptide["TD_10_seq"], "TD_10", 1, "boolean"),  #freq_1d, 
           #residue_distribution(peptide["TD_10_seq"], "TD_10", 2, "boolean", ),  #freq_2d,
           #residue_distribution(peptide["TD_10_seq"], "TD_10", 3, "boolean",),
           #residue_distribution(peptide["TD_10_seq"], "TD_10", 1, "absolute"),  #freq_1d, 
           #residue_distribution(peptide["TD_10_seq"], "TD_10", 2, "absolute", ),  #freq_2d,
           #residue_distribution(peptide["TD_10_seq"], "TD_10", 3, "absolute", ),
           
           residue_distribution(peptide["veltri_seq"], "veltri", 1, "distribution"),  #freq_1d, 
           residue_distribution(peptide["veltri_seq"], "veltri", 2, "distribution", ),  #freq_2d,
           #residue_distribution(peptide["veltri_seq"], "veltri", 3, "distribution", ),
           #residue_distribution(peptide["veltri_seq"], "veltri", 1, "boolean"),  #freq_1d, 
           #residue_distribution(peptide["veltri_seq"], "veltri", 2, "boolean", ),  #freq_2d,
           #residue_distribution(peptide["veltri_seq"], "veltri", 3, "boolean",),
           #residue_distribution(peptide["veltri_seq"], "veltri", 1, "absolute"),  #freq_1d, 
           #residue_distribution(peptide["veltri_seq"], "veltri", 2, "absolute", ),  #freq_2d,
           #residue_distribution(peptide["veltri_seq"], "veltri", 3, "absolute", ),
           
           di2(peptide["seq"], "aa"),  #Every third amino acid, useful for alpha-helices,
           di2(peptide["conjoint_seq"], "conjoint"),
           #di2(peptide["TD_10_seq"], "TD_10"),
           
           di3(peptide["conjoint_seq"]),  #Conjoint Alphabet #Every third amino acid, useful for alpha-helices,

           paac,
           hemopotent(peptide["seq"]),
           ##Positional information####################################################################################
           
           pepact,
           ), axis=1)
    
    #print(peptide["array"].shape)
    #print(residue_distribution(peptide["veltri_seq"], "veltri", 3, "distribution", ).shape)
    
    
    positional_size = 15
    peptide["alpha2num_veltri"] = alpha2num(peptide["veltri_seq"], alphabet=veltri_letters, padding=200, just="right")

    peptide["pos_first_aa_1"] = alpha2num(peptide["seq"][:positional_size], alphabet=aa_letters, padding = positional_size, just="left")

    peptide["pos_first_veltri"] = alpha2num(peptide["veltri_seq"][:positional_size], alphabet=veltri_letters, padding=15, just="left")
    
    peptide["pos_last_aa_1"] = alpha2num(peptide["seq"][-positional_size:], alphabet=aa_letters, padding=15, just="right")
    peptide["pos_last_veltri"] = alpha2num(peptide["veltri_seq"][-positional_size:], alphabet=veltri_letters, padding=15, just="right")
    
    peptide["property_pos_first_1"] = property_pos_first(peptide["seq"], 1, number=10)
    
    peptide["property_pos_last_1"] = property_pos_last(peptide["seq"], 1, number=10)
    
    peptide["property_pos_convx"] = property_pos_convx(peptide["seq"], 1, number=200)
    peptide["property_pos_convx_angular"] = property_pos_convx_angular(peptide["seq"], 1, number=200)
    
    for identifier, mode in aai_to_get:
        peptide[identifier+"_"+mode+"_first"] = property_pos_first_conv(peptide["seq"], identifier, 1, number=10)
        peptide[identifier+"_"+mode+"_last"] = property_pos_last_conv(peptide["seq"], identifier, 1, number=10) 
        
    return peptide