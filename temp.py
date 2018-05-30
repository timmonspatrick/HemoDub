# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 18:43:56 2018

@author: Patrick
"""

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

new = [i[0] for i in to_get]
new = list(set(new))
print(new)