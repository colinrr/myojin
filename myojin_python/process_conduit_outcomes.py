#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 21 12:32:45 2024

@author: crrowell
"""

from typing import Union, List

import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator, griddata
from scipy.ndimage import gaussian_filter
from skimage.measure import find_contours
from skimage.transform import resize
# from scipy.interpolate import NearestNDInterpolator
# from skimage.measure import label, regionprops

import myojin_python.config as config  


# Lookup functions for conduit outcome codes as a reference
def get_outcome_code_table(display=True, simplified=True):
    if simplified:
        table = pd.read_csv(config.SIMPLIFIED_OUTCOME_CSV, index_col='Row')
    else:
        table = pd.read_csv(config.OUTCOME_CSV, index_col='Row')
    
    
    if display:
        print(table)
    return table

# Super function to process a suite of files and produce contours
def sweep_set_codes(sweep_list: List[str], outcome_codes: dict, plot_codes: dict):
    
    # Parse name(s)

    
    # Process codes
    # Extract 2 contour lines for each
    # Make contour lines plot
    pass


# Processing rules for each parameter sweep file
 # 1 replace invalid effusive with effusive
 # 2 Nearest neighbour interp to replace errored results?
 # 3 Cut to a specified max pressure
 # 4 Conncomp to get remaining chunks
 #   -> get outlines
 #   -> Any gaps to fill?
 #   -> smooth outlines
 #   -> exclude regions of overlap?

# Function to run the workflow on each array of codes
def process_outcome_codes(outcome_codes: np.ndarray, plot_codes: np.ndarray, P, max_P=2e7) \
    -> tuple[np.ndarray,np.ndarray,np.ndarray]:
    """
    Take a 2D array of outcome codes from conduit a parameter sweep, and process
    them to create a 'clean' image
      - scrub/fill invalid and errored codes where appropriate (e.g. for known issue
        with the conduit search routine)
      - crop the regime space to a more imited pressure range

    Parameters
    ----------
    outcome_codes : np.ndarray
        Array of parameter sweep outcome codes.
    plot_codes : np.ndarray
        Array of simplified outcome codes used for convenience plotting.
    P : np.ndarray
        Vector of pressure values. (Sweeping over excess pressure is assumed)
    max_P : TYPE, optional
        Upper limit of P for cropping the regime space. The default is 2e7.

    Returns
    -------
    outcome_codes : np.ndarray
        Process outcome codes.
    plot_codes : np.ndarray
        Processes simplified codes.
    P : TYPE
        Cropped P vector.

    """
    
    # 1 replace invalid effusive with effusive - based on verification showing 
    # most/all invalid effusive cases resulted from a lack of search precision,
    # not lack of viable solution - 2024-09-21
    invalid_effusive = outcome_codes==-1
    outcome_codes[invalid_effusive] = 1
    plot_codes[invalid_effusive] = 1
    
    # 2 replace validFragPressBalanced with validExplosive, since we are mainly
    # interested in the explosive-effusive transition.
    valid_press_bal = outcome_codes==2
    outcome_codes[valid_press_bal]==3
    
    # 2 Nearest neighbour interp to replace errored results?
    plot_codes = clear_error_codes(outcome_codes, plot_codes)
    
    # 3 Cut to a specified max pressure
    outcome_codes, plot_codes, P = clip_max_p(outcome_codes, plot_codes, P)
    
    # Conncomp retrieval
    # region_label = label(plot_codes)
    # region_props
    
    return outcome_codes, plot_codes, P, #region_label

# 2
def clear_error_codes(outcome_codes: np.ndarray, plot_codes: np.ndarray) -> np.ndarray:
    """
    clear_error_codes(outcome_codes: np.ndarray, plot_codes: np.ndarray)
    Nearest neighbour interpolation to replace errored results.

    Parameters
    ----------
    outcome_codes : np.ndarray
        Conduit model OutcomeCodes.
    plot_codes : np.ndarray
        Simplified conduit OutcomeCode indices for plotting.

    Returns
    -------
    plot_codes : np.ndarray
        Plot codes cleared of errors for simple plotting.

    """
    get_interp = (outcome_codes <= -3)
    mask = np.where(~get_interp)
    
    # perform nearest neighhbour interpolation on missing values
    x_grid,y_grid = np.mgrid[0:outcome_codes.shape[0],0:outcome_codes.shape[1]]
    plot_codes = griddata(np.transpose(mask), plot_codes[mask],(x_grid,y_grid), method='nearest')
    # lin_grid = griddata(np.transpose(mask), plot_codes[mask],(x_grid,y_grid), method='linear')


    return plot_codes

# 3
def clip_max_p(outcome_codes: np.ndarray, plot_codes: np.ndarray, P: np.ndarray, max_P: float=2e7) \
    -> tuple[np.ndarray,np.ndarray,np.ndarray]:
    """
    Cut arrays to a specified max pressure
    Parameters
    ----------
    outcome_codes : np.ndarray
        Conduit model OutcomeCodes.
    plot_codes : np.ndarray
        Simplified conduit OutcomeCode indices for plotting.
    P: np.ndarray
        Vector of overpressure values.
    max_P: float
        Maximum pressure above which crop values from arrays.

    Returns
    -------
    outcome_codes : np.ndarray
        Conduit model OutcomeCodes cropped to max pressure.
    plot_codes : np.ndarray
        Plot codes cropped to max pressure.
    P : np.ndarray
        Cropped overpressure vector.

    """
    
    P_keep = P <= max_P
    outcome_codes = outcome_codes[P_keep,:]
    plot_codes = plot_codes[P_keep,:]
    P = P[P_keep]
    
    
    return outcome_codes, plot_codes, P


# 4 Get contours of connected components
def get_label_contours(img: np.ndarray, x: np.ndarray, y: np.ndarray, sigma: float = 2.) \
    -> tuple[np.ndarray,np.ndarray,np.ndarray]:
    """
    Get smoothed contours for the major indices in each plot code "image".

    Parameters
    ----------
    img : np.ndarray
        Image array of plot code indices.
    x : np.ndarray
        Image x vector.
    y : np.ndarray
        Image y vector.
    sigma : float, optional
        Standard deviation (in pixels) of gaussian filter for smoothing 
        boundaries. The default is 2.0.

    Returns
    -------
    contours_xy: np.ndarray
        [x y] column vectors, one for each unique value in img, giving
        coordinates of smoothed boundaries.
    unique_codes: np.ndarray
        Vector of unique values in img.

    """
    
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

