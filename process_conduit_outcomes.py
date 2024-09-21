#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 21 12:32:45 2024

@author: crrowell
"""

import numpy as np

# Processing rules
 # 1 Cut to a specified max pressure
 # 2 replace invalid effusive with effusive
 # 3 Nearest neighbour interp to replace errored results?
 # 4 Conncomp to get remaining chunks
 #   -> get outlines
 #   -> Any gaps to fill?
 #   -> smooth outlines
 #   -> exclude regions


# Master function to run the workflow on each array of codes
def process_outcome_codes(outcome_codes: np.ndarray, plot_codes: np.ndarray, P, max_P=2e7) -> (np.ndarray, dict):
    
    # 1 Cut to a specified max pressure
    outcome_codes, plot_codes, P = clip_max_p(outcome_codes, plot_codes, P)
    
    # 2 replace invalid effusive with effusive
    invalid_effusive = outcome_codes==-1
    outcome_codes[invalid_effusive] = 1
    plot_codes[invalid_effusive] = 1
    
    return outcome_codes, plot_codes, P


# 1
def clip_max_p(outcome_codes: np.ndarray, plot_codes: np.ndarray, P, max_P=2e7):
    
    P_keep = P <= max_P
    outcome_codes = outcome_codes[P_keep,:]
    plot_codes = plot_codes[P_keep,:]
    P = P[P_keep]
    
    
    return outcome_codes, plot_codes, P

# 2
def reset_invalid_effusive