# myojin_knoll

Scripts and function to process output of volcanic conduit model parameter sweeps, for a project on Myojin Knoll volcano.

## Project Organization

```
├── LICENSE            <- Open-source license
├── Makefile           <- Makefile with convenience commands like `make create_environment`
├── README.md          <- The top-level README for developers using this project.
├── data               <- Secondary data from processing model sweep outputs.            
│
├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
│                         the creator's initials, and a short `-` delimited description, e.g.
│                         `1.0-jqp-initial-data-exploration`.
│
├── pyproject.toml     <- Project configuration file with package metadata for 
│                         myojin_python and configuration for tools like black
│
├── references         <- Data dictionaries, manuals, and all other explanatory materials.
│   ├── OutcomeCodeTable.csv            <- Raw conduit model outcome code descriptions
│   └── simplifiedOutComeCodeTable.csv  <- Simplified/plot outcome code descriptions
│
├── # reports          <- Generated analysis as HTML, PDF, LaTeX, etc.
│   └── figures        <- Generated graphics and figures to be used in reporting
│
├── environment.yml    <- The requirements file for reproducing the analysis environment
│
├── config.m            <- Configuring matlab directories, path, other useful variables
│
├── conduit_sweep_scripts <- MATLAB scripts and configuration used to run the conduit model
│                            parameter sweeps. NOTE: the conduit model itself is not included.
│
├── sandbox_scripts <- Matlab/python scripts for general calculations and data exploration
│
├── plot_tools      <- Plotting help functions and colormaps
│
└── myojin_python   <- Source code for use in this project.
    │
    ├── __init__.py                 <- Makes myojin_python a Python module
    │
    ├── config.py                   <- Store useful variables and configuration
    │
    ├── # dataset.py                  <- Scripts to download or generate data
    │
    ├── mat_tools.py                <- Utility functions to transfer matlab data to usable python format
    │
    ├── process_conduit_outcomes.py <- Functions to extract conduit model outcome codes for plotting
    │
    └── # plots.py                    <- Code to create visualizations
```

--------

### Python environment
To build the python environment using conda, run the following command in the main project directory (See `Makefile`):

`make create_environment`

This will create a conda environment named 'hydroplume_emulator'. You can update the environement using
`make requirements` if you need to add new packages.

If not using conda, refer to the `environment.yml` file to build your requirements.

### Directories
To work with the raw conduit model output, change the DATA_DIR variable in `myojin_python/config.py` to your own data directory.

Secondary data sets generated from processing are in the `/data` folder.

# Technical Notes

## Scripts to get started
`sandbox_scripts/MyojinSamplesData.m`
 -> Run this to generate the myojin sample data table (`data/MyojinSampleDataTable.mat`)
 -> Will also plot solubility for the Myojin sample melt inclusions, and shows the Liu et al. (2005) solubility equation curve.

`conduit_sweep_scripts/getAllResultsTable.m`
 -> get a suite of summary data for the full set of parameter sweeps, and outputs the resulting data table to `data/ConduitSweepsDataTable`. 
 -> DATA VARIABLES:
  ---- Scalar input parameters ---- 
    Q           : Mass discharge rate (kg/s)
    Z0          : Conduit length (m)
    Zw          : Water depth (m)
    dP          : Chamber excess pressure (Pa)
    n0_excess   : Excess exsolved gas fraction (wt.%)
    N0          : Initial Bubble Number Density (must be non-zero for n0_excess>0)
    conduit_radius : (m)

  ---- Scalar output values ---- 
    OutcomeCode, simpleCode  : see notes on Model Outcome Codes below
    Zf          : Fragmentation depth (m)
    C0          : Initial water solubility of melt (wt.%)


  ---- Calculated output fields ---- 
    Z_total       : Conduit length plus water depth (m)
    Csol          : NOT CURRENTLY USED. Can revisit if solublity info
                   is needed to compare with melt inclusion data
    Cm_0          : INITIAL dissolved H20 (wt.%)
    Psat          : Max. supersaturation pressure = initial saturation pressure - p_final (Pa)
    t_ascent      : Total modeled ascent time (s)
    dP_dt_bar     : AVERAGE decompression rate over modeled rise time (Pa/s)
    BND_f         : FINAL bubble number density (m^-3)
    N_nucl        : Number of discrete nucleation events
    Z_nucl        : Onset depth of discrete nucleation events (m)
    dZ_nucl       : Difference between frag. depth and nucl. depth (if frag. occurred) (m)
    phi_0         : INITIAL porosity (fraction)
    phi_f         : FINAL porosity (at surfac/end) (fraction)

`conduit_sweep_scripts/makeSweepScatterPlots.m`
 -> spits out a series of scatter plots to summarize the total set of parameter sweep data

## Conduit Model Parameter Sweeps

## Conduit Model Outcome Codes

The conduit model contains a class called OutcomeCode, which interprets the physical result of the simulation. E.g. as explosive and overpressured at the vent, explosive and pressure balanced at the vent, effusive, or an as erroneous results, and so on. The codes come in 2 verions:
  1) Raw model outcome codes. These are fairly fine grained and allow for detailed interpretation of results, but are a bit difficult to visualize.
  2) Simplified or Plot codes, in which various raw codes have been lumped together to create a simplified description for easy interpretation and plotting.

For the complete table of codes and their description, see references/OutcomeCodeTable.csv and references/simplifiedOutcomeCodeTable.csv.
-> You can also access these in matlab by running ConduitOutcome.getTable
-> Or in python using the module function process_conduit_outcomes.get_outcome_code_table(display=True, simplified=<True/False>). (simplified = False to show raw model codes, True to show simplified plot codes)

### On simplifying outcome codes for plotting purposes
To create the simplified version of outcome codes for data plots, the following rules are applied to the raw model outcome codes:
    1)...
    2) ...