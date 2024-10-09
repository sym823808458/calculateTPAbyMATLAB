# calculateTPAbyMATLAB
 A workflow to calculate excited states using Gaussian, extract dipole moment information with Multiwfn, and process the data with MATLAB to compute Two-Photon Absorption (TPA).
# Excited-State Calculation and Two-Photon Absorption (TPA) Analysis

This repository contains a workflow that facilitates excited-state calculations using Gaussian, extracts dipole moment information using Multiwfn, and processes the extracted data with MATLAB to compute Two-Photon Absorption (TPA) properties.

## Overview

This workflow involves three main steps:

1. **Gaussian Excited-State Calculation**: Perform excited-state calculations using Gaussian and generate a formatted checkpoint (`.fchk`) file.
2. **Dipole Moment Extraction with Multiwfn**: Use Multiwfn to extract dipole moment information from the `.fchk` file.
3. **Data Processing with MATLAB**: Process the output `.txt` files from Multiwfn using MATLAB to extract relevant dipole moment information and calculate TPA.

## Requirements

- **Gaussian**: For performing excited-state calculations.
- **Multiwfn**: For extracting transition dipole moments. http://sobereva.com/multiwfn/
- **MATLAB**: For processing the `.txt` files and computing the TPA values.

## Step-by-Step Guide

### 1. Gaussian Excited-State Calculation

#### Step 1: Prepare the Gaussian Input File

Create a Gaussian input file (`.com`) that includes the necessary settings to perform excited-state calculations. Below is an example configuration:

#P TD(nstates=20) IOp(9/40=4) [other settings like functional and basis set]
TD(nstates=20): Specifies that 20 excited states will be computed.
IOp(9/40=4): Ensures detailed excited-state information is written to the output.
Functional and Basis Set: Customize as per your requirements (e.g., B3LYP/6-31G(d)).
Once the .com file is prepared, run the calculation in Gaussian to generate the .fchk (formatted checkpoint) file.

#### Step 2: Run the Gaussian Calculation
Submit the job to Gaussian, which will generate the necessary output files, including the .fchk file required for the next steps.

### 2. Extract Dipole Moment Information with Multiwfn
Once the Gaussian calculation is complete, use Multiwfn to extract the dipole moment data from the .fchk file.

#### Step 1: Open .fchk File in Multiwfn
Launch Multiwfn.exe.
Input the .fchk file generated from Gaussian.
Follow these steps:  
Input 18 (Excited state analysis).  
Input 5 (Analysis of transition dipole moments).  
Press Space to proceed log file.  
Input 2 to generate a detailed .txt file containing information about the dipole moments.  
After completing the steps above, Multiwfn will generate a .txt file containing all the relevant dipole moment information. This file will be used for further analysis in MATLAB.

### 3. Process the Generated .txt Files Using MATLAB
Once you have the .txt files extracted from Multiwfn, use the provided MATLAB script to process all .txt files in the current folder and extract the relevant information.

Step 1: Run the MATLAB Script
Ensure all the .txt files are located in the same folder as the MATLAB script. The script will:
Read each .txt file.
Extract dipole moment information.
Compute Two-Photon Absorption (TPA) values.

### Citation
If you use this workflow or the code provided in this repository for your research, please cite the following papers:

Metal-Organic Layers with an Enhanced Two-Photon Absorption Cross-Section and UP-Converted Emission.DOI:10.1021/acs.chemmater.0c03457  
Hexa-Branched Nanographenes with Large Two-Photon Absorption. DOI:10.1021/jacs.3c05662  
Interpretable Machine Learning of Two-Photon Absorption. DOI:10.1002/advs.202204902 
In this study, we developed models in MLatom and XACS that can directly predict two-photon absorption cross-sections using machine learning techniques.
### License
This project is licensed under the MIT License. See the LICENSE file for details.

### Contact
If you have any questions or need further assistance, please reach out to:
Dr.Yuming Su  
Email: 823808458@qq.com
