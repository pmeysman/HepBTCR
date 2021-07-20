#Classes for use in the Elias et al. Engerix study

import pandas as pd
import re
import numpy as np

def parse_imgt(gene):
    #Convert TRB1-1*00(xx) into TRB1-01 as per IMGT format
    gene = re.sub('\*00.+$','',gene)


    return gene

class TcrData:
    def __init__(self):
        self.raw = pd.DataFrame(columns=['tcrseq', 'cdr3', 'count', 'freq'])
        self.df = pd.DataFrame(columns=['tcrseq', 'cdr3', 'count', 'freq'])
        self.cdr3hash = dict()

    def read_immunoseq(self, file):
        data = pd.read_csv(file, sep='\t', low_memory=False,
                           usecols=["vGeneName", "aminoAcid", "jGeneName",'count (templates/reads)','frequencyCount (%)'])

        data = data[data['aminoAcid'].str.match(r'C[A-Z]+F', na=False)]

        self.raw['tcrseq'] = data['vGeneName'] + "\t" +  data['aminoAcid'] + "\t" + data['jGeneName']
        self.raw['cdr3'] = data['aminoAcid']
        self.raw['count'] = data['count (templates/reads)']
        self.raw['freq'] = data['frequencyCount (%)']
        self.raw = self.raw.dropna(axis=0)
        # print(self.df)

        self.df = self.raw.groupby(['cdr3']).agg({'count':'sum','freq':'sum','tcrseq':'sum'})

        #print(self.df)

    def read_mixcr(self,file, minFreq = 0):
        data = pd.read_csv(file, sep='\t', low_memory=False)

        data = data[data['aaSeqCDR3'].str.match(r'C[A-Z]+F', na=False)]

        data['vGeneName']=data['allVHitsWithScore'].apply(lambda vgene: parse_imgt(vgene))

        data['jGeneName']=data['allJHitsWithScore'].apply(lambda jgene: parse_imgt(jgene))


        self.raw['tcrseq'] = data['vGeneName'] + "\t" + data['aaSeqCDR3'] + "\t" + data['jGeneName']
        self.raw['cdr3'] = data['aaSeqCDR3']
        self.raw['count'] = data['cloneCount']
        self.raw['freq'] = data['cloneFraction']

        self.raw = self.raw.dropna(axis=0)

        self.raw = self.raw[self.raw['count'] >= minFreq]

        self.df = self.raw.groupby(['cdr3']).agg({'count': 'sum', 'freq': 'sum', 'tcrseq': 'sum'})
    def read_qiaseq(self,file,regexset = ""):
        data = pd.read_csv(file, sep="\t", low_memory=False)

        print("Read "+str(len(data)))

        data = data[data['chain'] == 'TRBC']

        if regexset != "":
            data = data[data['read set'].str.match(regexset, na=False)]


        data['vGeneName'] = 'TRB' + data['V-region'].astype(str)
        data['jGeneName'] = 'TRB' + data['J-region'].astype(str)

        self.raw['tcrseq'] = data['vGeneName'] + "\t" + data['CDR3 amino acid seq'] + "\t" + data['jGeneName']
        self.raw['cdr3'] = data['CDR3 amino acid seq']
        self.raw['count'] = data['reads']
        self.raw['freq'] = data['frequency']

        self.df = self.raw.groupby(['cdr3']).agg({'count': 'sum', 'freq': 'sum', 'tcrseq': 'sum'})

    def buildhash(self):
        for cdr in set(self.raw['cdr3']):
            for hash in (cdr[::2], cdr[1::2]):
                if hash not in self.cdr3hash:
                    self.cdr3hash[hash] = set()
                self.cdr3hash[hash].add(cdr)

    def hammingintersect(self,comparisonset,dist = 1):
        #Get all cdrs that are within a distance of "dist" form the comparison in self tcr set

        #Uses a hashing function to speed up, so will be unreliable at large distances.

        # Returns a set

        if len(self.cdr3hash) == 0:
            #If hash doesn't exist yet, we need to make it first
            self.buildhash()


        results = set()
        for cdr in comparisonset:
            for hash in (cdr[::2], cdr[1::2]):
                if hash in self.cdr3hash:
                    for cdrlist in self.cdr3hash[hash]:
                        if sum(ch1 != ch2 for ch1, ch2 in zip(cdr, cdrlist)) < dist:
                            #print(cdr + ' ' + cdrlist + ':' + str(sum(ch1 != ch2 for ch1, ch2 in zip(cdr, cdrlist))))
                            results.add(cdrlist)

        return results

    def hammingmatch(self,comparisonset,dist = 1):
        # Get all cdrs that are within a distance of "dist" from the self tcr set in the comparison set

        #Uses a hashing function to speed up, so will be unreliable at large distances.

        if len(self.cdr3hash) == 0:
            #If hash doesn't exist yet, we need to make it first
            self.buildhash()


        results = set()
        for cdr in comparisonset:
            for hash in (cdr[::2], cdr[1::2]):
                if hash in self.cdr3hash:
                    for cdrlist in self.cdr3hash[hash]:
                        if sum(ch1 != ch2 for ch1, ch2 in zip(cdr, cdrlist)) < dist:
                            #print(cdr + ' ' + cdrlist + ':' + str(sum(ch1 != ch2 for ch1, ch2 in zip(cdr, cdrlist))))
                            results.add(cdr)

        return results

    def createnetwork(self,dist=1):

        if len(self.cdr3hash) == 0:
            #If hash doesn't exist yet, we need to make it first
            self.buildhash()

        network = set()
        for hash in self.cdr3hash:
            if len(self.cdr3hash[hash]) >= 2:
                for cdr1 in self.cdr3hash[hash]:
                    for cdr2 in self.cdr3hash[hash]:
                        if cdr1 < cdr2:
                            if sum(ch1 != ch2 for ch1, ch2 in zip(cdr1, cdr2)) == dist:
                                network.add(cdr1 + "\t" + cdr2)

        return network

    def entropy(self, equivalent = True):

        def calculate_entropy(row):
            return row['freq'] * np.log(row['freq'])

        entropy = -sum(self.raw.apply(calculate_entropy, axis=1))

        if not equivalent:

            return entropy

        else:

            return entropy / np.log(len(self.raw))
