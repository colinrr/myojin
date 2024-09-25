#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 21 12:32:45 2024

@author: crrowell
"""

import numpy as np
# from scipy.interpolate import NearestNDInterpolator
from scipy.interpolate import RegularGridInterpolator, griddata
from scipy.ndimage import gaussian_filter
from skimage.measure import find_contours
# from skimage.measure import label, regionprops
from skimage.transform import resize

# Processing rules
 # 1 replace invalid effusive with effusive
 # 2 Nearest neighbour interp to replace errored results?
 # 3 Cut to a specified max pressure
 # 4 Conncomp to get remaining chunks
 #   -> get outlines
 #   -> Any gaps to fill?
 #   -> smooth outlines
 #   -> exclude regions


# Master function to run the workflow on each array of codes
def process_outcome_codes(outcome_codes: np.ndarray, plot_codes: np.ndarray, P, max_P=2e7) \
    -> (np.ndarray, np.ndarray, np.ndarray):
    
    
    # 1 replace invalid effusive with effusive
    invalid_effusive = outcome_codes==-1
    outcome_codes[invalid_effusive] = 1
    plot_codes[invalid_effusive] = 1
    
    # 2 Nearest neighbour interp to replace errored results?
    plot_codes = clear_error_codes(outcome_codes, plot_codes)
    
    # 3 Cut to a specified max pressure
    outcome_codes, plot_codes, P = clip_max_p(outcome_codes, plot_codes, P)
    
    # Conncomp retrieval
    # region_label = label(plot_codes)
    # region_props
    
    return outcome_codes, plot_codes, P, #region_label

# 2
def clear_error_codes(outcome_codes: np.ndarray, plot_codes: np.ndarray):
    
    get_interp = (outcome_codes <= -3)
    mask = np.where(~get_interp)
    
    # perform nearest neighhbour interpolation on missing values
    x_grid,y_grid = np.mgrid[0:outcome_codes.shape[0],0:outcome_codes.shape[1]]
    plot_codes = griddata(np.transpose(mask), plot_codes[mask],(x_grid,y_grid), method='nearest')
    # lin_grid = griddata(np.transpose(mask), plot_codes[mask],(x_grid,y_grid), method='linear')


    return plot_codes

# 3
def clip_max_p(outcome_codes: np.ndarray, plot_codes: np.ndarray, P, max_P=2e7):
    
    P_keep = P <= max_P
    outcome_codes = outcome_codes[P_keep,:]
    plot_codes = plot_codes[P_keep,:]
    P = P[P_keep]
    
    
    return outcome_codes, plot_codes, P


# 4 Get contours of connnected components
def get_label_contours(img: np.ndarray, x: np.ndarray, y: np.ndarray, sigma: float = 3.):
    
    def normalizer(A,a):
        A = (A - a.min())/np.ptp(a)
        return A, a.min(), np.ptp(a)

    def denormalizer(A, a_min, a_rng):
        A = A * a_rng + a_min
        return A
    
#### 1 interpolate the image to a finer resolution to sample between actual regions
    # > Get old and new grids
    X,Y = np.meshgrid(x,y) 
    
    # Upscale code image
    img_us = resize(img, np.array(img.shape)*2+1,mode='edge', order=0)
    ni,nj = img_us.shape
    
    # Normalize x and y values for skimage functions, then resize
    X_scaled, x_min, x_rng = normalizer(X, x)
    Y_scaled, y_min, y_rng = normalizer(Y, y)
    
    X_us = resize(X_scaled, np.array(img.shape)*2+1,mode='edge', order=1)
    Y_us = resize(Y_scaled, np.array(img.shape)*2+1,mode='edge', order=1)

#### Interpolate x,y meshgrids for later retrieval of contoured coordinates
    X_interp = RegularGridInterpolator((np.linspace(0, ni-1, ni), np.linspace(0, nj-1, nj)), X_us)
    Y_interp = RegularGridInterpolator((np.linspace(0, ni-1, ni), np.linspace(0, nj-1, nj)), Y_us)
    
#### Cycle through unique codes to get bounding contours for each, applying a 
    # gaussian filter to smooth the contours
    unique_codes = np.unique(img)
    contours = []
    contours_xy = []
    for code in unique_codes:
        mask = img_us==code
        contour   = find_contours(gaussian_filter(mask.astype(float), sigma), 0.5)
        contours += contour
        xy = []
        for subcon in contour:
            xy.append(np.array([ denormalizer(X_interp(subcon) , x_min, x_rng) , 
                                denormalizer(Y_interp(subcon) , y_min, y_rng) ]).transpose())
        contours_xy += xy
    

    return contours, contours_xy, unique_codes
