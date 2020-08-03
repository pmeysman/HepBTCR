#!/usr/local/bin/python3
# Script to generate scattr plots to compare day 0 and day 60 TCR repertoires
# Colors overlay identified peptide / pepetide pool specific TCRs

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re

from tcrdata import TcrData
from os import listdir
from os.path import isfile, join
from math import log


datadir = "../data/repTCRb/"
pepdir = "../data/peptideTCRab/"
figdir = "../results/timescatter/"
respfile = "../data/groups.csv"

# Get the responder type data

responders = pd.read_csv(respfile,sep=",")
responders = responders.set_index('Vaccinee')

# Get all rep files at time 0 from directory
files = [f for f in listdir(datadir) if (isfile(join(datadir, f)) and "B0.tsv" in f)]


results_list = []


for f in files:
    root = f[0:-7]
    print(root)
    pre = TcrData()
    pre.read_immunoseq(join(datadir, root + '_B0.tsv'))
    post = TcrData()
    post.read_immunoseq(join(datadir, root + '_B60.tsv'))

    # print(pre.df)

    # Calculate overlap

    overlaptcr = list(pre.df.index[np.in1d(pre.df.index, post.df.index, assume_unique=True)])

    print('B0-B60 overlap: ' + str(len(overlaptcr)))

    print(pre.df.loc[overlaptcr]['freq'])

    # Get peptide specific TCRs

    pepfiles = [p for p in listdir(pepdir) if (isfile(join(pepdir, p))
                                               and (root + '_' in p)
                                               and ('beta.txt' in p)
                                               and not ('VZV' in p))]

    pepdata = dict()
    peptotdatatmp = dict()

    for p in pepfiles:

        print(p)

        m = re.match("H\d+_(P+\d*)_([ab])", p)

        if m:

            peptide = m.group(1)
            sample = m.group(2)

            if sample not in peptotdatatmp:
                peptotdatatmp[sample] = dict()

            pep = TcrData()
            pep.read_mixcr(join(pepdir, p))

            if peptide not in peptotdatatmp[sample]:
                peptotdatatmp[sample][peptide] = set()
            peptotdatatmp[sample][peptide].update(set(pep.raw['cdr3'].unique()))

    peptotdata = dict()
    for peptide in peptotdatatmp['a']:
        print(peptide)
        #Only keep those TCRs that occur in both in vitro assays
        peptotdata[peptide] = peptotdatatmp['a'][peptide].intersection(peptotdatatmp['b'][peptide])

        #Calculate overlap between assay TCRs and overlaptcr (those that occur at day 0 and day 60)
        pepoverlap = peptotdata[peptide].intersection(set(overlaptcr))

        if(len(pepoverlap) > 0):
            if peptide not in pepdata:
                pepdata[peptide] = set()
            pepdata[peptide].update(pepoverlap)
            print(pepdata[peptide])
            print(len(pepdata[peptide]))


    # Calculate some statistics

    resdict = {}

    resdict['volunteer'] = root
    resdict['B0total'] = len(pre.df.index)
    resdict['B0entropy'] = pre.entropy()
    resdict['B60total'] = len(post.df.index)
    resdict['B60entropy'] = post.entropy()
    resdict['B0B60overlap'] = len(set(pre.df.index).intersection(set(post.df.index)))
    resdict['B0reads'] = pre.df['count'].sum()
    resdict['B60reads'] = post.df['count'].sum()
    resdict['responder'] = responders.loc[root]['Status_2']
    if 'PP' in peptotdata:
        resdict['PPTotal'] = len(peptotdata['PP'])
        if 'PP' in pepdata:
            resdict['PPIntersect'] = len(pepdata['PP'])
            resdict['PPIntInc'] = len([p for p in pepdata['PP']
                                       if pre.df.loc[p]['freq']
                                       < post.df.loc[p]['freq']])
            resdict['PPIntB0freq'] = sum([pre.df.loc[p]['freq'] for p in pepdata['PP']])
            resdict['PPIntB60freq'] = sum([post.df.loc[p]['freq'] for p in pepdata['PP']])
        else:
            resdict['PPIntersect'] = 0
            resdict['PPIntInc'] = 0
            resdict['PPIntB0freq'] = 0
            resdict['PPIntB60freq'] = 0
        resdict['PPB0'] = len(peptotdata['PP'].intersection(pre.df.index))
        resdict['PPB60'] = len(peptotdata['PP'].intersection(post.df.index))

    results_list.append(resdict)

    # Plot

    plt.scatter(pre.df.loc[overlaptcr]['freq'].apply(log),
                post.df.loc[overlaptcr]['freq'].apply(log),
                label="Aspecific")
    plt.xlabel(root + ' Time 0')
    plt.ylabel(root + ' Time 60')
    ax = plt.gca()
    low_x, high_x = ax.get_xlim()
    low_y, high_y = ax.get_ylim()
    plt.plot([max(low_x, low_y), min(high_x, high_y)], [max(low_x, low_y), min(high_x, high_y)], 'k--')

    colorcode = {'PP': '#FF0000'}
    for i in range(1, 60):
        colorcode['P'+str(i)] = '#' + str(int(100-i)) + str(int(30 + (i**2) % 70)) + '0' + str(int(i % 10))

    for peptide in pepdata:
        if peptide == 'PP':
            plt.scatter(pre.df.loc[pepdata[peptide]]['freq'].apply(log),
                        post.df.loc[pepdata[peptide]]['freq'].apply(log),
                        c=colorcode[peptide], label=peptide)

    lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(figdir + root + '_PP.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close()

    #quit()

resultsdf = pd.DataFrame(results_list)
resultsdf.to_csv(figdir + 'volunteercounts.txt')
