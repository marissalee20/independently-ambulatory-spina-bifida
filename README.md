# About the Project
This repository contains de-identified simulated walking data for 22 barefoot, independently ambulatory children with spina bifida and 17 children with typical development. It also contains code to generate simulations for each walking trial and code to analyze kinematic and kinetic differences between high-functioning individuals with spina bifida and individuals with typical development. Data and comprehensive simulation results are available with the rest of this repository at SimTK (TODO). The intent of this repository is to accompany our research article.

# Installation
Currently, this repo can only be built from source. To do so,

## Clone the repo:
```git clone https://github.com/stanfordnmbl/independently-ambulatory-spina-bifida
cd independently-ambulatory-spina-bifida
```

## Download the prerequisite packages:
```pip install -r requirements.txt```

# Experimental Data
Experimental data can be found on SimTK (TODO). The data types are:

- .trc Marker Data: unfiltered marker data in the OpenSim global reference frame (x: forward, y: up, z: left). Marker data largely follow the Plug-in-Gait model with thigh wand markers replaced by a single patella marker.
- .mot Ground Reaction Forces (GRFs): raw and lowpass filtered at 6 Hz.
- .xlsx EMG Data: raw EMG data.

Marker and ground reaction force data are available for static and overground walking trials. EMG data are available for 6 independently ambulatory children with spina bifida and 1 child with typical development.

# Musculoskeletal Model
The generic model is described in [Uhlrich et al. 2022](https://doi.org/10.1038/s41598-022-13386-9) and is based on the model described in [Rajagopal et al. 2016](https://doi.org/10.1109/TBME.2016.2586891). This model is available at SimTK: [https://simtk.org/projects/fbmodpassivecal/](https://simtk.org/projects/fbmodpassivecal/).

# Musculoskeletal Simulations
The musculoskeletal simulation pipeline is automated in MATLAB scripts. These scripts are largely based on previous scripts and code written by Apoorva Rajagopal and Scott Uhlrich. Their scripts can be viewed [here](https://simtk.org/projects/full_body) and [here](https://simtk.org/projects/coordretraining/).

To run the entire simulation pipeline for a subject, execute [code/run_all.m](https://github.com/marissalee20/independently-ambulatory-spina-bifida/blob/main/code/run_all.m). This script reads in motion data from the data directory and makes simulations for all walks associated with the subject. Simulation setup files and results are saved in the simulation directory.

Mean and standard deviation IK errors and RRA residuals can be calculated by executing [calculate_ik_errors.m](https://github.com/marissalee20/independently-ambulatory-spina-bifida/blob/main/code/calculate_ik_errors.m) and [calculate_rra_residuals](https://github.com/marissalee20/independently-ambulatory-spina-bifida/blob/main/code/calculate_rra_residuals.m), respectively.

# Group Post-Processed Data
All post-simulation analyses were performed in Python and can be found in [code/postprocessing_analysis.ipynb](https://github.com/marissalee20/independently-ambulatory-spina-bifida/blob/main/code/postprocessing_analysis.ipynb). Helper code can be found in [code/postprocessing_helpers](https://github.com/marissalee20/independently-ambulatory-spina-bifida/tree/main/code/postprocessing_helpers).

# Citation
Please cite our paper in your publications if our repository helps your research.

TODO

# Contact
Please feel welcome to reach out to Marissa Lee ([marissalee@stanford.edu](marissalee@stanford.edu)) with any questions.

# License
Distributed under the BSD 3-clause License. See LICENSE for more information.
