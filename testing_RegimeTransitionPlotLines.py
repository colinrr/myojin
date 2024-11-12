#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 21 12:23:48 2024

@author: crrowell
"""


from os.path import join

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
# from scipy.ndimage import gaussian_filter
# from skimage import measure
# from skimage.transform import resize
# from scipy.interpolate import RegularGridInterpolator

import myojin_python.mat_tools as mat
import myojin_python.process_conduit_outcomes as po
from myojin_python.config import DATA_DIR

dataDir   = DATA_DIR / 'refinedSweep/'
codesFile = 'outcomeCodeSummary.mat'

# % Test file
# % foi = '2024-04-26_myojin_Q8_Z0900_Zw500_357n_dP_21_n0_excess_17.mat'
# foi = 'myojin_Q8_Z0900_Zw500'
# foi = 'myojin_Q8_Z01100_Zw300'
# foi = 'myojin_Q9_Z01100_Zw300'
foi = 'myojin_Q9_Z05500_Zw500'


max_P = 2e7 # Max overpressure
sigma = 2.0 # Gaussian filter size (pix)


## ------------------------------------------------------------------------##

# --------  PULL DATA --------
matfile = join(dataDir,codesFile)

all_outcome_codes = mat.matfile_struct_to_dict(matfile, 'allOutcomeCodes')
simple_plot_codes = mat.matfile_struct_to_dict(matfile, 'simplePlotIndex')

n0 = all_outcome_codes['n0_excess']
P = all_outcome_codes['dP']


outcome_codes = all_outcome_codes[foi]
plot_codes = simple_plot_codes[foi]

# ------- PARSE SWEEP NAMES --------
# sweep_names = [key for key in outcome_codes.keys() if 'myojin' in key]

# for foi in sweep_names:
#     chunks = foi.split('_')
    


# --------  PROCESSING WORKFLOW --------
outcome_codes2, plot_codes2, P2 = po.process_outcome_codes(outcome_codes, plot_codes, P, max_P=max_P)

# Testing method to rescale image and get contours
# -> Method works, just needs repeating for each unique plot code, plus smoothing
# -> One method could be to:
    # a) get binary mask for single code
    # b) upscale - run same upscale on coordinate meshgrids
    # c) gaussian filter
    # d) get contour
    # e) assign n0,P coordinates using upscaled grids
    
contours, contours_xy, unique_codes = po.get_label_contours(plot_codes2, n0, P2, sigma=sigma)




# -------- PLOTTING/REPORTING --------
# print(outcome_codes.shape)
# print(P.shape)

# plt.imshow(plot_codes)
mesh_cmap = mpl.cm.Dark2.resampled(len(unique_codes))
alpha = 0.4

fig,ax = plt.subplots(1,2)

ax[0].pcolor(n0,P/1e6,plot_codes, cmap=mesh_cmap,alpha=alpha)
ax[0].set_ylabel('$\\Delta P$ (MPA)')
ax[0].set_xlabel('$n_{ex}$')

# CS = ax[1].pcolormesh(n0,P2/1e6,plot_codes2, shading='auto', edgecolor = '#eeefff')
# ax[1].contour(n0,P2/1e6,gaussian_filter(plot_codes2, .7), [0.5, 1.5], colors='k' )

CS = ax[1].pcolormesh(n0,P2/1e6,plot_codes2, shading='auto', linewidth=0.1, cmap=mesh_cmap, alpha=alpha)
# ax[1].contour(gaussian_filter(plot_codes2, .7), [0.5, 1.5], colors='k' )

line_colors = mpl.colormaps['tab10'](np.linspace(0,1,len(unique_codes)))
# for code,contour,xy in zip(unique_codes,contours,contours_xy):
for col,(code, con_group) in zip(mesh_cmap.colors,contours.items()):
    if con_group is not None:
        for con in con_group:
            # plt.plot(n0[contour[:, 1]], P[contour[:, 0]]/1e6, colors = 'b')
            ax[1].plot(con[:, 0], con[:,1]/1e6, color = col) # NOTE dividing by two because I haven't rescaled the contours

# ax[1].set_ylabel('$\\Delta P$ (MPA)')
ax[1].set_xlabel('$n_{ex}$')
plt.colorbar(CS)

# fig = plt.figure()
# plt.pcolor(n0,P2/1e6,lin_grid)
# plt.colorbar()

plt.show()
