# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 11:46:24 2018

@author: Patrick
"""
import aaindex
aaindex.init(path='.')

aai_to_get = [("CHAM810101", "mean"), #Steric Hinderance
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
                  
                  #("CHOP780212", "mean"),  #Frequency of first residue in turn
                  #("CHOP780213", "mean"),  #Frequency of 2nd residue in turn
                  #("CHOP780214", "mean"),  #Frequency of 3rd residue in turn
                  #("CHOP780215", "mean"),  #Frequency of 4th residue in turn
                  #("CHOP780212", "total"),  #Frequency of first residue in turn
                  #("CHOP780213", "total"),  #Frequency of 2nd residue in turn
                  #("CHOP780214", "total"),  #Frequency of 3rd residue in turn
                  #("CHOP780215", "total"),  #Frequency of 4th residue in turn
                  
                  
                                            #From HemoPi: src = https://github.com/riteshcanfly/Hemopi/blob/master/pcCalculator.java
                  ("HOPT810101", "mean"), #Hydrophilicity 
                  ("EISD840101", "mean"), #Hydrophobicity
                  ("FAUJ880109", "total"), #Net Hydrogen
                  ("EISD860101", "mean"), #Solvation
                  
                                          #From HemoPred
                  ("PONP800104", "mean"), #Hydrophobicity
                  ("NADH010103", "mean"), #Hydrophobicity
                  ("NADH010104", "mean"), #Hydrophobicity
                  ("MONM990101", "mean"), #Turn propernsity scale
                  ("AURR980106", "mean"), #Normalizied positional residue frequency at turn helix
                  
                  ("GUYH850104", "mean"),
                  ("GUYH850104", "max"), #Hydrophobicity, file:///C:/Users/Patrick/OneDrive%20-%20University%20College%20Dublin/Papers/680096.Juretic-Bioinformatics-of-Antimicrobial-Peptides.pdf
                  ("JANJ790101", "mean"),
                  ("JANJ790101", "max"),
                  
                ]