# Fed-GWAS Framework

## Overview
Fed-GWAS is a framework designed for conducting federated genome-wide association studies (GWAS). This project includes the implementation of various functionalities to allow users to analyze genetic data in a federated manner.

## Getting Started

To run the code, please follow the steps below:

1. Ensure you have all necessary dependencies installed. You may use `pip` to install them.
2. Open `autorun.py` in your preferred Python environment.

### Configuration
In the `autorun.py` file, you will need to configure the following parameters before running the program:
- **Number of miners**: Specify how many miners will be involved in the GWAS.
- **Number of SNPs**: Define the number of Single Nucleotide Polymorphisms (SNPs) to analyze.
- **Number of individuals**: Set the total number of individuals participating in the study.
- **Auxiliary variables**: If there are any auxiliary variables that need to be included, make sure to define them.

### Data Source
The data files required for the analysis are located in the `data` folder. This dataset consists of breast cancer data, which has been sourced from the Scikit-learn library.

## Installation

To install the necessary packages, you can use the following command:
```bash
pip install -r requirements.txt
