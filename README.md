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
│
├── # reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
│   └── figures        <- Generated graphics and figures to be used in reporting
│
├── environment.yml    <- The requirements file for reproducing the analysis environment
│
├── setup.cfg          <- Configuration file for flake8
│
├── conduit_sweep_scripts <- MATLAB scripts and configuration used to run the conduit model
│                            parameter sweeps. NOTE: the conduit model itself is not included.
│
├── sandbox_scripts <- Matlab/python scripts for general calculations and data exploration
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