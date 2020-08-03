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



results = dict()
results[pp] = dict()
results['responder'] = dict()
results['PPB0'] = dict()
results['PPB60'] = dict()
results['PPnrB0'] = dict()
results['PPnrB60'] = dict()
results['iPPB0'] = dict()
results['iPPB60'] = dict()
results['PSB0'] = dict()
results['PSB60'] = dict()
results['IPSB0'] = dict()
results['IPSB60'] = dict()
results['B0'] = dict()
results['B60'] = dict()
results['CD154'] = dict()

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
        for loosample in tcrdata[peptide]:
            if loosample == volunteer:
                if peptide == 'PP':
                    preiPP.update(set(pre.raw['cdr3'].unique()).intersection(tcrdata[peptide][loosample]))
                    postiPP.update(set(post.raw['cdr3'].unique()).intersection(tcrdata[peptide][loosample]))
            else:
                if peptide == 'PP':
                    if not responders.loc[loosample]['Status_2'] == 'Early-converter':
                        prePPnr.update(set(pre.raw['cdr3'].unique()).intersection(tcrdata[peptide][loosample]))
                        postPPnr.update(set(post.raw['cdr3'].unique()).intersection(tcrdata[peptide][loosample]))
                        PPnr.update(tcrdata[peptide][loosample])
                    prePP.update(set(pre.raw['cdr3'].unique()).intersection(tcrdata[peptide][loosample]))
                    postPP.update(set(post.raw['cdr3'].unique()).intersection(tcrdata[peptide][loosample]))
                    PP.update(tcrdata[peptide][loosample])

                else:
                    prePS = set(pre.raw['cdr3'].unique()).intersection(tcrdata[peptide][loosample])
                    preOverlap.update(prePS)
                    postPS = set(post.raw['cdr3'].unique()).intersection(tcrdata[peptide][loosample])
                    postOverlap.update(postPS)
                    PS.update(tcrdata[peptide][loosample])

    results['PSB0'][volunteer] = len(preOverlap) / len(PS)
    results['PSB60'][volunteer] = len(postOverlap) / len(PS)

    results['PPB0'][volunteer] = len(prePP) / len(PP)
    results['PPB60'][volunteer] = len(postPP) / len(PP)

    results['PPnrB0'][volunteer] = len(prePPnr) / len(PPnr)
    results['PPnrB60'][volunteer] = len(postPPnr) / len(PPnr)

    results['iPPB0'][volunteer] = len(preiPP)
    results['iPPB60'][volunteer] = len(postiPP)

    results['IPSB0'][volunteer] = len(preOverlap.intersection(prePP))
    results['IPSB60'][volunteer] = len(postOverlap.intersection(postPP))

    results['B60'][volunteer] = len(post.raw['cdr3'].unique())
    results['B0'][volunteer] = len(pre.raw['cdr3'].unique())


    results['responder'][volunteer] = responders.loc[volunteer]['Status_2']

    if volunteer in cd154.index:
        results['CD154'][volunteer] = cd154.loc[volunteer]['response']
    else:
        results['CD154'][volunteer] = 'NA'


respd = pd.DataFrame.from_dict(results)
respd.to_csv(resdir + 'loocvrep_ab.txt')
