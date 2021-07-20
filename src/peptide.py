#!/usr/local/bin/python3
# Peptide-specific analysis
# End result of the script is a csv file that can be further analyzed with the attached R script (because plotting in ggplot is easier/nicer)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re

from tcrdata import TcrData
from os import listdir
from os.path import isfile, join


datadir = "../data/repTCRb/"
pepdir = "../data/peptideTCRab/"
resdir = "../results/peptide/"
respfile = "../data/groups.csv"
cd154file = "../data/freqCD154.txt"

pp = 'PP'

# Setting the c threshold from the paper at the reviewers request
# Set based on prior Bioinformatics paper
# Higher values than 1 will cause approximations as the hashing will not be perfect
distance_set = 1

# Get the responder type data

responders = pd.read_csv(respfile,sep=",")
responders = responders.set_index('Vaccinee')

# Get CD154 assay data

cd154 = pd.read_csv(cd154file,sep=" ")
cd154 = cd154.set_index('Vaccinee')

# Read the peptide-specific TCR data sets

pepfiles = [p for p in listdir(pepdir) if (isfile(join(pepdir, p))
                                           and ('beta.txt' in p)
                                           and not ('VZV' in p))]
tcrdatatmp = dict()
for p in pepfiles:
    # Regex to match (donor) (peptide) (replicate)
    m = re.match("(H\d+)_(P+\d*)_([abcd])_", p)

    if m:
        print(m)
        volunteer = m.group(1) + m.group(3)
        peptide = m.group(2)
        data = TcrData()
        data.read_mixcr(join(pepdir,p),minFreq=0)

        tcrset = set(data.raw['cdr3'].unique())


        if peptide in tcrdatatmp:
            if volunteer in tcrdatatmp[peptide]:
                tcrdatatmp[peptide][volunteer].update(tcrset)
            else:
                tcrdatatmp[peptide][volunteer] = tcrset
        else:
            tcrdatatmp[peptide] = dict()
            tcrdatatmp[peptide][volunteer] = tcrset

tcrdata = dict()

# Only keep those TCRs that occur in both replicates (a/b or c/d if available)

for peptide in tcrdatatmp:
    for volunteer in tcrdatatmp[peptide]:
        if 'a' in volunteer:

            root = volunteer.replace('a', '')
            c = volunteer.replace('a','c')
            if c in tcrdatatmp[peptide]: #If c/d is available, these are preferred over a/b
                volunteer = c
                b = volunteer.replace('a','d')
            else:
                b = volunteer.replace('a','b')

            if peptide not in tcrdata:
                tcrdata[peptide] = dict()

            tcrdata[peptide][root] = tcrdatatmp[peptide][volunteer].intersection(tcrdatatmp[peptide][b])


# This is the output feature matrix
results = dict()
results[pp] = dict() # Count of the peptide pool TCRs
results['responder'] = dict() # Responder status
results['B0'] = dict() # CD4+ memory breadth at Time 0
results['B60'] = dict() # CD4+ memory breadth at Time 60
results['CD154'] = dict() # CD154 values

# Values used to calculate the Rhbs metric =  numerator / denominator (done in R script)
results['PPnrB0'] = dict() # Leave-one-out non-early-responder peptide pool match with Time 0 (denominator)
results['PPnrB60'] = dict() # Leave-one-out non-early-responder peptide pool match with Time 60 (denominator)
results['PSB0'] = dict() # Leave-one-out peptide-specific TCR matches with Time 0 (numerator)
results['PSB60'] = dict() # Leave-one-out peptide-specific TCR matches with Time 60 (numerator)

# Values used for additional plots
# Not to be used for any predictions as they were not run in a cross validation format!
results['iPPB0'] = dict() # Own peptide pool match with Time 0
results['iPPB60'] = dict() # Own peptide pool match with Time 60
results['IPSB0'] = dict() # Overlap between peptide pool and peptide-specific matches at Time 0
results['IPSB60'] = dict() # Overlap between peptide pool and peptide-specific matches at Time 60


# Do for each repertoire file

files = [f for f in listdir(datadir) if (isfile(join(datadir, f)) and "B0.tsv" in f)]

for f in files:
    root = f[0:-7]
    print(root)
    pre = TcrData()
    pre.read_immunoseq(join(datadir, root + '_B0.tsv'))
    post = TcrData()
    post.read_immunoseq(join(datadir, root + '_B60.tsv'))
    volunteer = root
    results[pp][volunteer] = len(tcrdata[pp][volunteer])

    preOverlap = set()
    postOverlap = set()
    prePP = set()
    postPP = set()
    prePPnr = set()
    postPPnr = set()
    preiPP = set()
    postiPP = set()

    PS = set()
    PP = set()
    PPnr = set()


    for peptide in tcrdata:

        #Leave-one-out cross validation, where each sample is left out exactly once
        for loosample in tcrdata[peptide]:

            #If sample is the sample that is left out, do not do anything with it,
            #except checking the iPP variable (which is just a check for the own PP-specific TCRs)
            if loosample == volunteer:
                if peptide == 'PP':
                    preiPP.update(pre.hammingintersect(tcrdata[peptide][loosample],dist = distance_set))
                    postiPP.update(post.hammingintersect(tcrdata[peptide][loosample],dist = distance_set))

            #This is the main leave-one-out loop
            #All other samples are collected as training data
            #Hamming distances are calculated and compared
            #Sets are kept with unique matching TCRs
            else:
                if peptide == 'PP':
                    #This is the peptide pool data
                    if not responders.loc[loosample]['Status_2'] == 'Early-converter':
                        #This is the non-early converter samples

                        #This builds the denominator for Day 0 (pre) and Day 60 (post)
                        prePPnr.update(pre.hammingintersect(tcrdata[peptide][loosample],dist = distance_set))
                        postPPnr.update(post.hammingintersect(tcrdata[peptide][loosample],dist = distance_set))
                        PPnr.update(tcrdata[peptide][loosample])

                    #Extra data for additional plots
                    prePP.update(pre.hammingintersect(tcrdata[peptide][loosample],dist = distance_set))
                    postPP.update(post.hammingintersect(tcrdata[peptide][loosample],dist = distance_set))
                    PP.update(tcrdata[peptide][loosample])

                else:
                    #This is the epitope specific data
                    #Used to build the numerator
                    prePS = set(pre.hammingintersect(tcrdata[peptide][loosample],dist = distance_set))
                    preOverlap.update(prePS)
                    postPS = set(post.hammingintersect(tcrdata[peptide][loosample],dist = distance_set))
                    postOverlap.update(postPS)
                    PS.update(tcrdata[peptide][loosample])

    #Calculate the Overlap Coefficient for each set.
    #This is saved to a file and passed on to the R scripts

    #Numerator
    results['PSB0'][volunteer] = len(preOverlap) / len(PS)
    results['PSB60'][volunteer] = len(postOverlap) / len(PS)

    #Denominator
    results['PPnrB0'][volunteer] = len(prePPnr) / len(PPnr)
    results['PPnrB60'][volunteer] = len(postPPnr) / len(PPnr)

    #Own matches to make expansion plots
    results['iPPB0'][volunteer] = len(preiPP)
    results['iPPB60'][volunteer] = len(postiPP)

    #PS and PP overlap to make overlap plots
    results['IPSB0'][volunteer] = len(preOverlap.intersection(prePP))
    results['IPSB60'][volunteer] = len(postOverlap.intersection(postPP))

    #Total breadth
    results['B60'][volunteer] = len(post.raw['cdr3'].unique())
    results['B0'][volunteer] = len(pre.raw['cdr3'].unique())

    #Save responder status for easy parsing
    results['responder'][volunteer] = responders.loc[volunteer]['Status_2']

    #Save CD154 value for asy parsing
    if volunteer in cd154.index:
        results['CD154'][volunteer] = cd154.loc[volunteer]['response']
    else:
        results['CD154'][volunteer] = 'NA'


respd = pd.DataFrame.from_dict(results)
respd.to_csv(resdir + 'loocvrep_ab.txt')
