# -*- coding: utf-8 -*-
"""
Created on Thu Mar 01 19:14:52 2018

@author: Patrick

A set of translation tables for reduced alphabets
"""

TD_3 = {'A': 'V', 'C': 'C', 'E': 'E', 'D': 'E', 'G': 'E', 'F': 'V', 'I': 'V', 'H': 'E', 'K': 'E', 'M': 'V', 'L': 'V', 'N': 'E', 'Q': 'E', 'P': 'E', 'S': 'E', 'R': 'E', 'T': 'E', 'W': 'V', 'V': 'V', 'Y': 'V'}

TD_5 = {'A': 'V', 'C': 'C', 'E': 'E', 'D': 'E', 'G': 'G', 'F': 'V', 'I': 'V', 'H': 'G', 'K': 'R', 'M': 'V', 'L': 'V', 'N': 'G', 'Q': 'G', 'P': 'G', 'S': 'G', 'R': 'R', 'T': 'G', 'W': 'V', 'V': 'V', 'Y': 'V'}

TD_10 = {"V" : "V", "I" : "V", "L" : "V", "M" : "V", "F" : "V", "H" : "H", "Q" :"H", "N" : "H", "C" :"C", "E": "E","D":"E","R":"R","K":"R","A":"A","G":"G","W":"W","Y":"W","P":"P","S":"S","T":"T"}

conjoint_dict = {"A" : "A", "G" : "A", "V" : "A",               #A dictionary to translate amino acids
                 "I" : "I", "L" : "I", "F" : "I", "P" : "I",        #to conjoint amino acids
                 "Y" : "Y", "M" : "Y", "T" : "Y", "S" : "Y",
                 "H" : "H", "N" : "H", "Q" : "H", "W" : "H",
                 "R" : "R", "K" : "R",
                 "D" : "D", "E" : "D",
                 "C" : "C",}
				 
veltri = {"H" : "H", "N" : "H",
		  "Y" : "Y", "P" : "Y", "R" : "Y", "K" : "Y",
		  "G" : "G", "W" : "G",
		  "A" : "A", "V" : "A", "F" : "A",
		  "L" : "L", "I" : "L", "M" : "L",
		  "Q" : "Q", "S" : "Q", "T" : "Q",
		  "D" : "D", "E" : "D",
         "C" : "C",}

all_alphabets = {"TD_3" : TD_3,
                 "TD_5" : TD_5,
                 "TD_10" : TD_10,
                 "conjoint_dict" : conjoint_dict,
                 "veltri" : veltri}