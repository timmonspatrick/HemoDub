# -*- coding: utf-8 -*-
"""
Created on Wed May 30 18:50:15 2018

@author: Patrick
"""
from di2 import *
import pytest

def test_di2():  ##Need to specify alphabet
    assert sum(di2("ACDEFH", alphabet="aa")[0]) == 1
    
def test_di3():  #Use only conjoint letters
    assert 0.99 < sum(di3("AIYHRDCIYIYHRDC")[0]) < 1.01