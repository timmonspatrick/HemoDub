# -*- coding: utf-8 -*-
"""
Created on Thu May 03 15:10:45 2018

@author: Patrick

Source idea: file:///C:/Users/Patrick/OneDrive%20-%20University%20College%20Dublin/Papers/680096.Juretic-Bioinformatics-of-Antimicrobial-Peptides.pdf
"""
from D_descriptors import D_descriptor
import pytest

def test_D_descriptor():
    assert D_descriptor("FGPAC") == (0.7760238852276552, 4.937809062252897, 7.283689452238277, 8.799664219270753, -0.6330958041584763, 7.283689452238277, 7.311151915662141)