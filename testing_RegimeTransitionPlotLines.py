#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 21 12:23:48 2024

@author: crrowell
"""

import mat_tools as mat
import process_conduit_outcomes as po
from os.path import join
import matplotlib.pyplot as plt

dataDir   = '/Users/crrowell/Kahuna/data/myojin/mainSweep2/refinedSweep/'
codesFile = 'outcomeCodeSummary.mat'

# % Test file
# % foi = '2024-04-26_myojin_Q8_Z0900_Zw500_357n_dP_21_n0_excess_17.mat';
foi = 'myojin_Q8_Z0900_Zw500';

max_P = 2e7 # Max overpressure
## ------------------------------------------------------------------------##


matfile = join(dataDir,codesFile)

all_outcome_codes = mat.matfile_struct_to_dict(matfile, 'allOutcomeCodes')
simple_plot_codes = mat.matfile_struct_to_dict(matfile, 'simplePlotIndex')

n0 = all_outcome_codes['n0_excess']
P = all_outcome_codes['dP']


outcome_codes = all_outcome_codes[foi]
plot_codes = simple_plot_codes[foi]

outcome_codes, plot_codes, P = po.process_outcome_codes(outcome_codes, plot_codes, P, max_P=max_P)

print(outcome_codes.shape)
print(P.shape)

# plt.imshow(plot_codes)

plt.pcolor(n0,P/1e6,plot_codes)
plt.xlabel('$\\Delta P$ (MPA)')
plt.ylabel('$n_{ex}$')

plt.show()
