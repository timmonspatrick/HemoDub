# -*- coding: utf-8 -*-
"""
Created on Sat Feb 24 13:35:49 2018

@author: Patrick
"""
from modules.modlamp.descriptors import PeptideDescriptor, GlobalDescriptor


def pep_descriptor(seq):

    pepdesc = PeptideDescriptor(seq, "eisenberg")
    pepdesc.calculate_global(modality="mean", append=True)
    pepdesc.calculate_global(modality="max", append=True)
    pepdesc.calculate_moment(modality="max", append=True)
    pepdesc.calculate_moment(modality="mean", append=True)
    #pepdesc.calculate_profile(append=True, prof_type = "uH")
    
    pepdesc.load_scale("Ez")
    pepdesc.calculate_global(modality="mean", append=True)
    pepdesc.calculate_global(modality="max", append=True)
    
    pepdesc.load_scale("aasi")
    pepdesc.calculate_global(append=True)
    pepdesc.calculate_moment(modality="max", append=True)
    pepdesc.calculate_moment(modality="mean", append=True)
    
    pepdesc.load_scale("abhprk")
    pepdesc.calculate_global(modality="mean", append=True)
    pepdesc.calculate_global(modality="max", append=True)
    
    pepdesc.load_scale("charge_acid")
    pepdesc.calculate_global(modality="mean", append=True)
    pepdesc.calculate_global(modality="max", append=True)
    pepdesc.calculate_moment(modality="max", append=True)
    pepdesc.calculate_moment(modality="mean", append=True)
    
    pepdesc.load_scale("cougar")
    pepdesc.calculate_global(modality="mean", append=True)
    pepdesc.calculate_global(modality="max", append=True)
    
    pepdesc.load_scale("gravy")
    pepdesc.calculate_global(modality="mean", append=True)
    pepdesc.calculate_global(modality="max", append=True)
    pepdesc.calculate_moment(modality="max", append=True)
    pepdesc.calculate_moment(modality="mean", append=True)
    
    pepdesc.load_scale("hopp-woods")
    pepdesc.calculate_global(modality="mean", append=True)
    pepdesc.calculate_global(modality="max", append=True)
    pepdesc.calculate_moment(modality="max", append=True)
    pepdesc.calculate_moment(modality="mean", append=True)
    
    pepdesc.load_scale("kytedoolittle")
    pepdesc.calculate_global(modality="mean", append=True)
    pepdesc.calculate_global(modality="max", append=True)
    pepdesc.calculate_moment(modality="max", append=True)
    pepdesc.calculate_moment(modality="mean", append=True)
    
    pepdesc.load_scale("ppcali")
    pepdesc.calculate_global(modality="mean", append=True)
    pepdesc.calculate_global(modality="max", append=True)
    
    pepdesc.load_scale("msw")
    pepdesc.calculate_global(modality="mean", append=True)
    pepdesc.calculate_global(modality="max", append=True)
    
    pepdesc.load_scale("charge_phys")
    pepdesc.calculate_moment(modality="max", append=True)
    pepdesc.calculate_moment(modality="mean", append=True)
    pepdesc.calculate_global(modality="mean", append=True)
    pepdesc.calculate_global(modality="max", append=True)
    
    pepdesc.load_scale("flexibility")
    pepdesc.calculate_moment(modality="max", append=True)
    pepdesc.calculate_moment(modality="mean", append=True)
    pepdesc.calculate_global(modality="mean", append=True)
    pepdesc.calculate_global(modality="max", append=True)
    
    pepdesc.load_scale("bulkiness")
    pepdesc.calculate_moment(modality="max", append=True)
    pepdesc.calculate_moment(modality="mean", append=True)
    pepdesc.calculate_global(modality="mean", append=True)
    pepdesc.calculate_global(modality="max", append=True)
    
    pepdesc.load_scale("TM_tend")
    pepdesc.calculate_moment(modality="max", append=True)
    pepdesc.calculate_moment(modality="mean", append=True)
    pepdesc.calculate_global(modality="mean", append=True)
    pepdesc.calculate_global(modality="max", append=True)
    
    pepdesc.load_scale("mss")
    pepdesc.calculate_moment(modality="max", append=True)
    pepdesc.calculate_moment(modality="mean", append=True)
    pepdesc.calculate_global(modality="mean", append=True)
    pepdesc.calculate_global(modality="max", append=True)
    
    pepdesc.load_scale("t_scale")
    pepdesc.calculate_global(modality="mean", append=True)
    pepdesc.calculate_global(modality="max", append=True)
    
    pepdesc.load_scale("peparc")
    pepdesc.calculate_arc(modality="max", append=True)
    pepdesc.calculate_arc(modality="mean", append=True)        

    pepdesc.load_scale("msw")
    pepdesc.calculate_global(modality="mean", append=True)
    pepdesc.calculate_global(modality="max", append=True)
    
    pepdesc.load_scale("polarity")
    pepdesc.calculate_moment(modality="max", append=True)
    pepdesc.calculate_moment(modality="mean", append=True)
    pepdesc.calculate_global(modality="mean", append=True)
    pepdesc.calculate_global(modality="max", append=True)
    
    pepdesc.load_scale("pepcats")
    pepdesc.calculate_global(modality="mean", append=True)
    pepdesc.calculate_global(modality="max", append=True)
    
    pepdesc.load_scale("isaeci")
    pepdesc.calculate_global(modality="mean", append=True)
    pepdesc.calculate_global(modality="max", append=True)

    pepdesc.load_scale("refractivity")
    pepdesc.calculate_moment(modality="max", append=True)
    pepdesc.calculate_moment(modality="mean", append=True)
    pepdesc.calculate_global(modality="mean", append=True)
    pepdesc.calculate_global(modality="max", append=True)
    
    pepdesc.load_scale("z3")
    pepdesc.calculate_global(modality="mean", append=True)
    pepdesc.calculate_global(modality="max", append=True)
    
    pepdesc.load_scale("z5")
    pepdesc.calculate_global(modality="mean", append=True)
    pepdesc.calculate_global(modality="max", append=True)
    
    return pepdesc