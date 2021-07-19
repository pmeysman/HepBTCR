# Preexisting memory CD4 T cells in naïve individuals confer robust immunity upon vaccination

Code accompanying the publication "Preexisting memory CD4 T cells in naïve individuals confer robust immunity upon vaccination".

To replicate the entire analysis, the full data set needs to be downloaded and all processing steps are required. Figures and results can be regenerated using the processed files available within this repository (starting from step 3).

## Step 1: Download data

Data set is stored separately in Zenodo:
https://doi.org/10.5281/zenodo.3989144

This data set needs to be extracted and placed in the folder 'data'.

The data folder should then contain the following:
* peptideTCRab folder: A folder containing the raw MiXCR output of peptide-specific TCR clonotypes, named by donor, peptide, replicate, sample id and sequencing lane.
* repTCRb folder: A folder containing the raw immunoSEQ output of CD4+ memory repertoire-level TCR clonotypes, named by donor, time point.
* df_all.tsv: A table with meta and FC data for samples and volunteers
* freqCD154.txt: A table with CD154-based group definitions
* groups.csv: A table with group definitions used for classes
* CD137_CD154_assay_data.tsv: A table with the CD137, CD154 assay data

## Step 2: Run python processing scripts

The 'src' folder contains all scripts necessary to replicate the processing steps described in the manuscript.

The three main python scripts are:
* tcrdata.py: A library containing TCR-specific functions and modules that are used in the other scripts.
* timescatter.py: A script to process the longitudinal data analysis steps. Outputs the processed file 'volunteercounts.txt' in the results folder.
* peptide.py: A script to process the epitope-specific TCR data (in a leave-one-out-cross validation) and generate the feature matrix for further analysis. Outputs the processed file 'loocvrep_ab.txt' in the results folder.

The analyses can be run with the following command on a device that has python3 installed.
> python src/timescatter.py
> python scr/peptide.py

As the output files are already available in the github repository, they will be overwritten by running the python scripts.

## Step 3: Run R visualisation scripts

The 'src' folder contains two R script used to generate the main manuscript images and supplementary materials. They can be run starting from the data that is already provided in the git repository.
* main figures.R: Script to generate the main figures
* supp figures.R: Script to generate the supplementary figures.
