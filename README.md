# Preexisting memory CD4 T cells in naïve individuals confer robust immunity upon vaccination

Code accompanying the publication "Preexisting memory CD4 T cells in naïve individuals confer robust immunity upon vaccination".

To replicate the entire analysis, the full data set needs to be downloaded and all processing steps are required. Figures and results can be regenerated using the processed files available within this repository (starting from step 3).

When run correctly, the code should result in two new folders, namely 'figures' and 'tables'.

## Step 1: Download Sequencing data

The full sequencing data set is stored separately in Zenodo:
https://doi.org/10.5281/zenodo.3989144

This data set needs to be extracted and placed in the folder 'data'.

The data folder should then contain the following:
* peptideTCRab folder: A folder containing the raw MiXCR output of peptide-specific TCR clonotypes, named by donor, peptide, replicate, sample id and sequencing lane.
* repTCRb folder: A folder containing the raw immunoSEQ output of CD4+ memory repertoire-level TCR clonotypes, named by donor, time point.
* df_all.tsv: A table combining meta and FC data for samples and volunteers
* metadata.csv: A table with volunteer metadata
* freqCD154.txt: A table with CD154-based group definitions
* groups.csv: A table with group definitions used for classes
* CD137_CD154_assay_data.tsv: A table with the CD137, CD154 assay data
* all_diversity.tsv: A table
* Antibody and Cytokine Data.csv: A table with FCM antibody and cytokine data
* Cell count.csv: A table with FCM cell counts
* Lymphocytes.csv: A table with FCM lymphocyte data
* Subset_specific_sign.csv: A table with FCM subset data
* Table (CD127 MFI).csv: A table with CD127 MFI
* Table (CD137 MFI).csv: A table with CD137 MFI
* Table (CD154 MFI).csv: A table with CD154 MFI
* Table (CD25 MFI).csv: A table with CD25 MFI
* Table_1.csv: A table with FCM data (set 1)
* Table_2.csv: A table with FCM data (set 2)
* Table_3.csv: A table with FCM data (set 3)
* Table_4.csv: A table with FCM data (set 4)
* Table_5.csv: A table with FCM data (set 5)
* Table_Treg_profile.csv: A table with regulatory T-cell profile data

## Step 2: Run processing scripts

The 'src' folder contains all scripts necessary to replicate the processing steps described in the manuscript.

The three main processing scripts for the sequencing data are:
* tcrdata.py: A library containing TCR-specific functions and modules that are used in the other scripts.
* timescatter.py: A script to process the longitudinal data analysis steps. Outputs the processed file 'volunteercounts.txt' in the results folder.
* peptide.py: A script to process the epitope-specific TCR data (in a leave-one-out-cross validation) and generate the feature matrix for further analysis. Outputs the processed file 'loocvrep_ab.txt' in the results folder.

The analyses can be run with the following command on a device that has python3 installed.

> python src/timescatter.py

> python scr/peptide.py

The processing script to parse the flow cytometry data is:
* FCM_analysis.R: A script to combine all FCM measurements into a single table called 'df_1.tsv'.

It can be run with:

> R src/FCM_analysis.R

As the output files are already available in the github repository, they will be overwritten by running the processing scripts.

## Step 3: Run R visualisation scripts

The 'src' folder contains R scripts used to generate the main manuscript images and supplementary materials. They can be run starting from the data that is already provided in the git repository.
* Fig_#.R: Script to generate the main figures
* Fig_Supp_###.R: Script to generate the supplementary figures (from FCM or TCR data).
* Tables.R: Scrpt to generate tables for the manuscript
