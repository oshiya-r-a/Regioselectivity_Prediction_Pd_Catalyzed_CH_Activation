# Machine Learning Predicts Regioselectivity in Pd-Catalyzed Directing Group-Assisted C–H Activation
![Image](https://github.com/user-attachments/assets/14c958e7-d7dd-4835-93fc-b1e34e2248f7)


# Quick Links
  1. [Overview](#overview)
  2. [Repository Structure](#repository-structure)
  3. [Requirements](#requirements)
  4. [Getting Started](#getting-started)
  6. [Citation](#citation)

## Overview
This repository accompanies our study where machine learning models were developed to predict the regioselective outcome of palladium-catalyzed directing group-assisted C–H activation reactions. Using datasets curated from molecular fingerprints derived from the SMILES representations of reactant molecules, we trained and validated multiple models (SVM and Logistic Regression) to capture subtle structural features influencing site selectivity. The pretrained models provided here allow users to make regioselectivity predictions on new molecules using `.mol` file as input.

For full experimental details and data curation methodology, please refer to the original publication [citation](#citation) below.

## Repository Structure
```
├── notebooks/                     # Google Colab notebooks for each model
│   ├── svm_standard.ipynb            
│   ├── svm_large_fragments.ipynb           
│   ├── lr_standard.ipynb           
│   └── lr_large_fragments.ipynb           
│
├── omp_predictor_main/            # Python package directory
|   ├── setup.py 
|   └── omp_predictor/    
│       ├── __init__.py
│       ├── fingerprint.py         # Generates Morgan fingerprints
│       ├── encorder.py            
│       ├── omp.py                 # Core prediction functions
│       ├── svm_model.pkl          # Pretrained ML model
│       ├── lf_svm_model.pkl                
│       ├── lg_model.pkl         
│       └── lf_lg_model.pkl          
│
├── dataset/                
│   ├── omp-data-326.csv           # Dataset used to train the model
│   └── mol_files/                 # Sample .mol files for testing
│
├── requirements.txt              # Required Python packages
└── README.md                     # Project documentation 
```

## Requirements
All required Python packages are listed in the [`requirements.txt`](./requirements.txt) file.  
To install them, run:
```
pip install -r requirements.txt
```

## Getting Started
This section will guide you on how to clone the repository, install the package, and use it to make regioselectivity predictions on your local machine (Windows/Linux) or Google Colab.

### Clone the Repository
To get started, clone the repository using:
```
git clone https://github.com/oshiya-r-a/omp-predictor.git
```
Alternatively, you can download `omp-predictor-main` as a ZIP file and extract it.

### Install the Package (Windows/Linux)
Navigate into the package directory and install it using:
```
cd omp-predictor-main
pip install .
```

Once installed, you can use the `omp` command from your terminal to make regioselectivity predictions on molecules provided in `.mol` format.
Place all the `.mol` files to be tested inside a single folder (e.g., a folder named `mol_files`), and run:
```
omp ./mol_files
```

**Optional Flags:**

```--lf``` : Uses the Large Fragment encoder (by default, the standard Morgan fingerprint is used)

```--logreg``` : Uses the Logistic Regression model (default is SVM)

**Examples:**

```
omp ./mol_files --lf --logreg   # Uses large fragment + logistic regression
omp ./mol_files --logreg        # Uses Morgan fingerprint + logistic regression
omp ./mol_files --lf            # Uses large fragment + svm
```

⚠️ Note: Ensure that all required dependencies are installed. The list of dependencies can be found in the [`requirements.txt`](./requirements.txt) file.

## Citation
The methodology and model implementation used in this repository are described in detail in the following publication: <br /> 
<br />
R. A. Oshiya; Mohari A.; Datta A. Machine Learning Predicts Regioselectivity in Pd-Catalyzed Directing Group-Assisted C–H Activation. *Org. Lett.* 2025, 27, (19), 4909-4914. DOI: https://doi.org/10.1021/acs.orglett.5c01158 <br />
<br />
If you use any part of this code, model, or the associated dataset in your work, please cite the above publication.
