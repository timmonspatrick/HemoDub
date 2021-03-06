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
from reduced_alphabets import conjoint_dict, TD_3, TD_5, TD_10

##########HEMOLYTIK########################
hemolytik_peptides = {}

for filename in os.listdir("HEMOLYTIK"):
    if "HEMOLYTIK" in filename:
        with open("HEMOLYTIK\{}".format(filename),"r", encoding="utf-8") as f:
            text = f.read()
            natural = "Non-Natural</b></td><td><FONT  COLOR='black'>None"
            if "Primary information" in text and natural in text and "Ψ" not in text:
                comments = text.split("ACTIVITY</b></td><td><FONT  COLOR='black'>")[1].split("</b></td></tr>")[0]
                sequence = text.split("SEQUENCE</b></td><td><FONT  COLOR='black'>")[1].split("</b></td></tr>")[0]
                sequence = sequence.replace("/","")
                
                if "X" not in sequence and "O" not in sequence and "Z" not in sequence:

                    sequence_len = int(text.split("LENGTH</b></td><td><FONT  COLOR='black'>")[1].split("</b></td>")[0])

                    cTer = text.split("C-ter Modification</b></td><td><FONT  COLOR='black'>")[1].split("</b></td>")[0]
                    nTer = text.split("N-ter Modification</b></td><td><FONT  COLOR='black'>")[1].split("</b></td>")[0]
                    
                    linear = "COLOR='#5D325A'>Linear/Cyclic</b></td><td><FONT  COLOR='black'>Linear</b>" in text
                    l_stereo = "Stereochemistry</b></td><td><FONT  COLOR='black'>L</b>" in text
                    antimicrobial = "NATURE</b></td><td><FONT  COLOR='black'>Antimicrobial</b>" in text
                    
                    if l_stereo and antimicrobial and "*" not in sequence and "spacer" not in sequence and "-" not in sequence and sequence == sequence.upper() and len(sequence) > 6:
                        mw = molecular_weight(sequence, seq_type="protein")
                        hemolytic_category = hemolytik_parser(comments, mw)

                        #hemolytik_peptides.append({"seq" : sequence, "nTer" : nTer, "cTer" : cTer, "seq_len" : sequence_len, "activity" : hemolytic_category, "id" : filename.replace(".html","")})
                        key = (sequence, nTer, cTer, sequence_len)
                        
                        if key not in hemolytik_peptides:
                            hemolytik_peptides[key] = {"activity" : hemolytic_category, "id" : filename.replace(".html","")}
                        elif key in hemolytik_peptides:
                            if hemolytik_peptides[key]["activity"] == "YES":
                                pass
                            elif hemolytik_peptides[key]["activity"] == "NO":
                                if hemolytic_category == "YES":
                                    hemolytik_peptides[key] = {"activity" : hemolytic_category, "id" : filename.replace(".html","")}
                        

hemolytik_peptides = [{"seq" : k[0], "nTer" : k[1], "cTer" : k[2], 
                       "seq_len" : k[3], "activity" : hemolytik_peptides[k]["activity"], 
                       "id" : hemolytik_peptides[k]["id"],}
                        for k in hemolytik_peptides]

random.shuffle(hemolytik_peptides)

j_s = {"Peptides" : hemolytik_peptides}


#################DBAASP##############################
dbaasp_peptides = []

for f in os.listdir('DBAASP'):
    if "DBAASP" in f:
        with open("DBAASP\{}".format(f), "r", encoding="utf-8") as file:
            
            try:
                text = file.read()
                data = json.loads(text)
                assert "seq" in data["peptideCard"]
                                     
                sequence = data["peptideCard"]["seq"]
                
                if "X" not in sequence and "O" not in sequence and "Z" not in sequence and sequence == sequence.upper() and len(sequence) > 6:
                    # If the sequence length can't be parsed to the integer, it's obviously an empty record
                    
                    if "monomers" in data["peptideCard"]:
                        assert len(data["peptideCard"]["monomers"]) < 2
                    
                                        
                    sequence_len = data["peptideCard"]["seqLength"]            
                
                    if "hemoliticCytotoxicActivities" in data["peptideCard"]:                                            
                        
                        hemolytic_activities = data["peptideCard"]["hemoliticCytotoxicActivities"]
                        
                        hemolytic_activities = [c for c in hemolytic_activities if "erythrocyte" in c["targetCell"].lower()]
                        
                        
                        
                        if len(hemolytic_activities) != 0:      

                            mw = molecular_weight(sequence, seq_type="protein")
                            hemolytic_category = dbaasp_parser(hemolytic_activities, mw)
                            
                            if data["peptideCard"]["cTerminus"] == "AMD":
                                cTer = "Amidation" 
                            else:
                                cTer = "None"
                                
                            dbaasp_peptides.append({"seq" : sequence, 
                                                   "nTer" : "None", 
                                                   "cTer" : cTer ,
                                                   "seq_len" : sequence_len, 
                                                   "activity" : hemolytic_category, 
                                                   "id" : f.replace(".html","")})
            except AssertionError:
                pass
                

################FASTA FILES FROM HEMOPI##############################################
filenames = [("pos.fa", "YES", True), ("neg.fa", "NO", True,), ("pos_val.fa", "YES", False), ("neg_val.fa", "NO", False)]
fasta_peptides = []



for filename, activity, train in filenames:

    with open(filename,"r", encoding="utf-8") as f:
        text = f.read()
        
        for line in text.split("\n"):
            print(line)
            if ">" not in line and len(line) > 2:
                line = line.replace(" ","")
                fasta_peptides.append({"seq" : line, "nTer" : 0, "cTer" : 0, 
                   "seq_len" : len(line), "activity" : activity, 
                   "id" : "0",
                   } )
            
##################MERGE BOTH SETS###################
all_peptides = fasta_peptides.copy()

if True:
    fasta_seqs = [e["seq"] for e in fasta_peptides]
    
    for peptide in dbaasp_peptides:
        if peptide["seq"] not in fasta_seqs:
            all_peptides.append(peptide)

##################MERGE BOTH SETS###################

if True:
    all_seqs = [e["seq"] for e in all_peptides]
    
    for peptide in hemolytik_peptides:
        if peptide["seq"] not in all_seqs:
            all_peptides.append(peptide)
    
random.shuffle(all_peptides)

for peptide in all_peptides:
    peptide["conjoint_seq"] = conjoint_seq(peptide["seq"], alphabet=conjoint_dict)
    peptide["TD_3_seq"] = conjoint_seq(peptide["seq"], alphabet=TD_3)
    peptide["TD_5_seq"] = conjoint_seq(peptide["seq"], alphabet=TD_5)
    peptide["TD_10_seq"] = conjoint_seq(peptide["seq"], alphabet=TD_10)

j_s = {"Peptides" : all_peptides}
with open("peptides.json","w") as f:
	f.write(json.dumps(j_s))
    