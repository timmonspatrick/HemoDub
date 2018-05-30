# -*- coding: utf-8 -*-
"""
Takes as input the HTML pages in the subdirectories HEMOLYTIK and DBAASP.
Gives as output a file named peptides.json
"""

import os
from comment_parser import hemolytik_parser, dbaasp_parser
import json
import random
from Bio.SeqUtils import molecular_weight
from helper_functions import conjoint_seq
from residue_distribution import alphabets
from reduced_alphabets import conjoint_dict, TD_3, TD_5, TD_10, veltri
from tqdm import tqdm

hemopi = "hemopi3"

filenames = [(hemopi + "_pos.fa", "YES", True), (hemopi+"_neg.fa", "NO", True,), (hemopi+"_pos_val.fa", "YES", False), (hemopi+"_neg_val.fa", "NO", False)]

##########HEMOLYTIK########################
fasta_peptides = []



for filename, activity, train in tqdm(filenames):

    with open(filename,"r", encoding="utf-8") as f:
        text = f.read()
        
        for line in text.split("\n"):
            print(line)
            if ">" not in line and len(line) > 2:
                line = line.replace(" ","")
                fasta_peptides.append({"seq" : line, "nTer" : 0, "cTer" : 0, 
                   "seq_len" : len(line), "activity" : activity, 
                   "id" : "0",
                   "conjoint_seq" : conjoint_seq(line, alphabet=conjoint_dict),
                   "TD_3_seq": conjoint_seq(line, alphabet=TD_3),
                   "TD_5_seq": conjoint_seq(line, alphabet=TD_5),
                   "TD_10_seq": conjoint_seq(line, alphabet=TD_10),
                   "veltri_seq": conjoint_seq(line, alphabet=veltri)
                   } )

random.shuffle(fasta_peptides)

j_s = {"Peptides" : fasta_peptides}
with open("peptides_%s.json" % (hemopi,),"w") as f:
	f.write(json.dumps(j_s))
    