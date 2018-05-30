# -*- coding: utf-8 -*-
"""
Created on Thu May 03 15:10:45 2018

@author: Patrick

Source idea: file:///C:/Users/Patrick/OneDrive%20-%20University%20College%20Dublin/Papers/680096.Juretic-Bioinformatics-of-Antimicrobial-Peptides.pdf
"""
from __future__ import division
from aaf import aaf
from math import cos, sin, degrees, radians
import aaindex

aaindex.init(path='.')

def D_descriptor(seq):
    aaf_Jvs = []
    aaf_Gvs = []
    for aa_n in range(len(seq)):
        aa = seq[aa_n]
        angle = float(90/(len(seq) - 1)) * aa_n
        aaf_J = aaf(aa, "JANJ790101", "mean")
        
        aaf_Jv = (aaf_J * cos(radians(angle)) , aaf_J * sin(radians(angle)))
 
        
        aaf_Jvs.append(aaf_Jv)
        
        aaf_G = aaf(aa, "GUYH850104", "mean")
        aaf_Gv = (aaf_G * cos(radians(angle)) , aaf_J * sin(radians(angle)))
        
        aaf_Gvs.append(aaf_Gv)
    
    aaf_GvT = ( sum([i[0] for i in aaf_Gvs]), sum([i[1] for i in aaf_Gvs]) )
    aaf_GvTl = (aaf_GvT[0]**2 + aaf_GvT[1]**2)**0.5
    aaf_JvT = ( sum([i[0] for i in aaf_Jvs]), sum([i[1] for i in aaf_Jvs]) )
    aaf_JvTl = (aaf_JvT[0]**2 + aaf_JvT[1]**2)**0.5

    
    D_desc = (aaf_GvT[0] * aaf_JvT[0] + aaf_GvT[1] * aaf_JvT[1])/(aaf_GvTl*aaf_JvTl)
 
    return D_desc, aaf_JvT[0], aaf_JvT[1], aaf_JvTl, aaf_GvT[0], aaf_GvT[1], aaf_GvTl
D_descriptor("FGPAC")

