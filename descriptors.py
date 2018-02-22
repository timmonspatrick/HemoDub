# -*- coding: utf-8 -*-

from __future__ import print_function
import json
from modlamp.descriptors import PeptideDescriptor, GlobalDescriptor
import numpy as np
from di2 import di2, di3
from pydpi.pypro import PyPro
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from aaf import aaf
from PyBioMed.PyProtein import CTD, ConjointTriad
import random
from helper_functions import conjoint_seq, feature_occurence, residue_distribution, residue_boolean, residue_abs

random.seed(31) ## Set a random seed for reproducible results

path = "peptides.json"   #This file is produced by running read_files.py
with open(path, "r") as f:
    text = f.read()
    peptides = eval(text)["Peptides"]
    random.shuffle(peptides)



def describe_sequences(peptides):    
            
    for peptide in peptides:
        peptide["conjoint_seq"] = conjoint_seq(peptide["seq"])
    
    dp, tp = feature_occurence(peptides) ##A function that returns the number of peptides that each di or tri-peptide occurs in
                                         ##Useful to prevent overfitting later on.
    
    for peptide in peptides:
        
        ###Calculate all the Global descriptors##################################
        globdesc = GlobalDescriptor(peptide["seq"])
        globdesc.calculate_all(amide = peptide["cTer"] == "Amidation")
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
        
        conjoint_dis = residue_distribution(peptide["conjoint_seq"], 4, None, None)
        #Frequency of the conjoint amino acids in the peptide sequence
        ####################################################################################################
        
        ## Calculate all the peptide descriptors############################################################
        pepdesc = PeptideDescriptor(peptide["seq"], "eisenberg")
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        pepdesc.calculate_moment(modality="max", append=True)
        pepdesc.calculate_moment(modality="mean", append=True)
        #pepdesc.calculate_profile(append=True, prof_type = "uH")
        
        pepdesc.load_scale("Ez")
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        
        pepdesc.load_scale("aasi")
        pepdesc.calculate_global(append=True)
        pepdesc.calculate_moment(modality="max", append=True)
        pepdesc.calculate_moment(modality="mean", append=True)
        
        pepdesc.load_scale("abhprk")
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        
        pepdesc.load_scale("charge_acid")
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        pepdesc.calculate_moment(modality="max", append=True)
        pepdesc.calculate_moment(modality="mean", append=True)
        
        pepdesc.load_scale("cougar")
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        
        pepdesc.load_scale("gravy")
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        pepdesc.calculate_moment(modality="max", append=True)
        pepdesc.calculate_moment(modality="mean", append=True)
        
        pepdesc.load_scale("hopp-woods")
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        pepdesc.calculate_moment(modality="max", append=True)
        pepdesc.calculate_moment(modality="mean", append=True)
        
        pepdesc.load_scale("kytedoolittle")
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        pepdesc.calculate_moment(modality="max", append=True)
        pepdesc.calculate_moment(modality="mean", append=True)
        
        pepdesc.load_scale("ppcali")
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        
        pepdesc.load_scale("msw")
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        
        pepdesc.load_scale("charge_phys")
        pepdesc.calculate_moment(modality="max", append=True)
        pepdesc.calculate_moment(modality="mean", append=True)
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        
        pepdesc.load_scale("flexibility")
        pepdesc.calculate_moment(modality="max", append=True)
        pepdesc.calculate_moment(modality="mean", append=True)
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        
        pepdesc.load_scale("bulkiness")
        pepdesc.calculate_moment(modality="max", append=True)
        pepdesc.calculate_moment(modality="mean", append=True)
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        
        pepdesc.load_scale("TM_tend")
        pepdesc.calculate_moment(modality="max", append=True)
        pepdesc.calculate_moment(modality="mean", append=True)
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        
        pepdesc.load_scale("mss")
        pepdesc.calculate_moment(modality="max", append=True)
        pepdesc.calculate_moment(modality="mean", append=True)
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        
        pepdesc.load_scale("t_scale")
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        
        pepdesc.load_scale("peparc")
        pepdesc.calculate_arc(modality="max", append=True)
        pepdesc.calculate_arc(modality="mean", append=True)        

        pepdesc.load_scale("msw")
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        
        pepdesc.load_scale("polarity")
        pepdesc.calculate_moment(modality="max", append=True)
        pepdesc.calculate_moment(modality="mean", append=True)
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        
        pepdesc.load_scale("pepcats")
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        
        pepdesc.load_scale("isaeci")
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
    
        pepdesc.load_scale("refractivity")
        pepdesc.calculate_moment(modality="max", append=True)
        pepdesc.calculate_moment(modality="mean", append=True)
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        
        pepdesc.load_scale("z3")
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        
        pepdesc.load_scale("z5")
        pepdesc.calculate_global(modality="mean", append=True)
        pepdesc.calculate_global(modality="max", append=True)
        
        
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
        
        freq_1d = residue_distribution(peptide["seq"], 1, dp, tp)
        freq_2d = residue_distribution(peptide["seq"], 2, dp, tp)
        freq_3d = residue_distribution(peptide["seq"], 3, dp, tp)
        freq_1dbool = residue_boolean(peptide["seq"], 1, dp)
        freq_2dbool = residue_boolean(peptide["seq"], 2, dp)
        freq_1dabs = residue_abs(peptide["seq"], 1, dp)
        freq_2dabs = residue_abs(peptide["seq"], 2, dp)
        
        len_peptide = np.array([[len(peptide["seq"])]])
        
        if peptide["activity"] == "YES":
            pepact = 1
        else:
            pepact = 0
        pepact = np.array([[pepact]])
        
        peptide_di2 = di2(peptide["seq"])
        peptide_di3 = di3(peptide["conjoint_seq"])
        
        ####################### AAindex #########################
        to_get = [("CHAM810101", "mean"), #Steric Hinderance
                ("CHAM810101", "total"), #Steric Hinderance
                  ("KYTJ820101", "mean"), #Hydropathy
                  ("KLEP840101", "total"), #Charge
                  ("KLEP840101", "mean"), #Charge
                  ("MITS020101", "mean"), #Amphiphilicity
                  ("FAUJ830101", "mean"), #Hydrophobic parameter pi
                  ("GOLD730102", "total"), #Residue volume
                  ("MEEJ800101", "mean"), #Retention coefficient in HPLC
                  ("OOBM850105", "mean"), #Optimized side chain interaction parameter
                  ("OOBM850105", "total"), #Optimized side chain interaction parameter
                  ("VELV850101", "total"), #Electron-ion interaction parameter
                  ("VELV850101", "mean"), #Electron-ion interaction parameter
                  ("PUNT030102", "mean"), #Knowledge-based membrane-propensity scale from 3D_Helix
                  ("BHAR880101", "mean"), #Average flexibility indeces
                  ("KRIW790102", "mean"), #Fraction of site occupied by water
                  ("PLIV810101", "mean"), #Partition coefficient
                  ("ZIMJ680102", "mean"), #Bulkiness
                  ("ZIMJ680102", "total"), #Bulkiness
                  ("ZHOH040101", "mean"), #Stability scale
                  ("CHAM820102", "total"), #Free energy solubility in water
                                          #From HemoPi: src = https://github.com/riteshcanfly/Hemopi/blob/master/pcCalculator.java
                  ("HOPT810101", "mean"), #Hydrophilicity 
                  ("EISD840101", "mean"), #Hydrophobicity
                  ("FAUJ880109", "total"), #Net Hydrogen
                  ("EISD860101", "mean"), #Solvation
                ]
        
        aminoacidindeces = np.array([[aaf(peptide["seq"], identifier, mode) for identifier, mode in to_get]])
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
        
        peptide["array"] = np.concatenate((pepid, pepdesc.descriptor, globdesc.descriptor, len_peptide, 
               cTer,
               secondary_structure_fraction,
               aminoacidindeces,
               ctdc_vals,
               conjointtriad_vals, 
               tp_bin, 
               conjoint_dis,
               freq_1d, 
               freq_2d,
               freq_3d,
               freq_1dbool,
               freq_2dbool,
               freq_1dabs,
               freq_2dabs,
               peptide_di2,  #Every third amino acid, useful for alpha-helices,
               peptide_di3,  #Conjoint Alphabet #Every third amino acid, useful for alpha-helices,
               paac,
               pepact,), axis=1)
        #print(peptide["TotalDescriptor"])
    
    
    x = np.concatenate([peptide["array"] for peptide in peptides], axis=0)
    
    np.save("peptides_array", x, allow_pickle=False)
    
describe_sequences(peptides)
print("complete")