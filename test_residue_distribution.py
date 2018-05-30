# -*- coding: utf-8 -*-
"""
Created on Wed May 30 18:33:14 2018

@author: Patrick
"""

from residue_distribution import *
import pytest
import numpy as np

def test_counter():
    assert counter("ACEACEGP", "aa", 1) == {"A" : 2.0, "C" : 2.0, "E" : 2.0, "G" : 1.0, "P" : 1.0,
                                            "D" : 0.0, "F" : 0.0, "H" : 0.0, "I" : 0.0, "K" : 0.0,
                                            "V" : 0.0, "Y" : 0.0, "L" : 0.0, "M" : 0.0, "N" : 0.0,
                                            "Q" : 0.0, "R" : 0.0, "S" : 0.0, "T" : 0.0, "W" : 0.0,
                                            }
    
def test_residue_distribution():
    assert list(residue_distribution("ACEACEGP", "aa", 1, "distribution")[0]) == [ 0.25,0.25,0, 0.25,0,0.125,0,0,0,0,0,0,0.125,0,0,0,0,0,0,0]
    assert list(residue_distribution("ACEACEGP", "aa", 1, "boolean")[0]) == [1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]
    
    assert list(residue_distribution("VCVCVC", "TD_3", 1, "distribution")[0]) == [0.5, 0, 0.5]
    
    assert 0.99 < sum(list(residue_distribution("ACACEE", "aa", 2, "distribution")[0])) < 1.01
