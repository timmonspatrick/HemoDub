# -*- coding: utf-8 -*-
"""
Created on Wed May 30 22:28:35 2018

@author: Patrick
"""
from helper_functions import *
import pytest

def test_alpha2num():
    assert list(alpha2num("ACE", padding=4, just="right")[0]) == [0, 1, 2, 4]
    
