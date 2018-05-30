# -*- coding: utf-8 -*-
"""
Created on Wed May 30 18:26:35 2018

@author: Patrick
"""

from aaf import aaf
import pytest

def test_aaf():
    assert aaf("GHACEF", "KLEP840101", "mean") == -0.16666666666666666
    assert aaf("GHACEF", "KLEP840101", "total") == -1.0
    assert aaf("GHACEF", "KLEP840101", "max") == 0.0
    
