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
                                            }
    
def test_residue_distribution():
    assert residue_distribution("ACEACEGP", "aa", 1, "distribution") == np.array([[ 0.25,0.25,0, 0.25,0,0.125,0,0,0,0,0,0,0.125,0,0,0,0,0,0,0]])
    assert residue_distribution("ACEACEGP", "aa", 1, "boolean") == np.array([[1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]])
    
    assert residue_distribution("VCVCVC", "TD_3", 1, "distribution") == np.array([[0,0.5,0.5]])