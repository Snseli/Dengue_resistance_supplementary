# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 10:27:56 2024

@author: selin.sen

This script takes and nucleotide alignment of DENV2 NS4B gene sequences as an input 
and outputs a table showing the percentage of sequences exhibiting one or more of 
the mosnodenvir resistance mutations described in Kaptein, et al. Nature (2021) doi:10.1038/s41586-021-03990-6.

usage python spot_D2_resistance_mutations.py aln.fasta
"""
from __future__ import (absolute_import, division, print_function)
import pandas as pd
import numpy as np
from Bio.Seq import Seq
import argparse # To deal with arguments :
# https://docs.python.org/2/library/argparse.html

def split_seq(seq,size=1):
    """ From sequence return a list of patterns of a given size
    """
    return [seq[i:i+size] for i in range(0, len(seq), size)]

def validNT(nucleotide):
    """ Returns True if the amino acid is not undetermined, not a gap, not a stop 
    """
    return (nucleotide != 'n') and (nucleotide != '-') and (nucleotide != '?')

####################### PARSE ARGUMENTS #######################
# Here we read the argument given and we check them
parser = argparse.ArgumentParser(
    description="This script is used to create a mutation table from a nucleotide alignment (in fasta format)")
parser.add_argument("aln",type=str,
                    help="path : path to the .fasta file to process")

parser.add_argument('-v','--verbose', action='store_true',
    help='Display more information')

args = parser.parse_args()

print("Initializing ...")

# Extract csv filename  and path
nameOfFastaFile = args.aln.split("/")[-1]
nameFastaFile = nameOfFastaFile[:-6]
pathToFastaFile = "/".join(args.aln.split("/")[:-1])
if pathToFastaFile == "":
    pathToFastaFile = "."

# Loading files : sequence alignment
print("Processing fasta file ...")

infile = open(args.aln, 'r+')
content = infile.readlines() #Reading sequence file in txt format to a file named content

# Extracting sequences and names to lists
content = [w.replace('\n', '') for w in content]#takes out \n from file
names=list() #creates empty list named names
sequences=list() #creates empty list named sequences
for i in range(0,len(content),2):
    names.append(content[i])#fills list with sequence names
    sequences.append(content[i+1]) #fills list with sequence

names = [i.replace('>', '') for i in names]

# Building basic name/sequence table and checking alignment is in frame
sequences=pd.Series(sequences) #creating panda series from sequences list
names=pd.Series(names) #creating panda series from names list

Nseq = sequences.count() #Nseq=number of sequences in pd.series
Nnt  = len(sequences[0])
tempList = [[""] * (2) for i in range(Nseq)]

print(f'There are {Nseq} sequences')
print(f'Length of the alignment is {Nnt} nt')

# Creating dataframe
for i in range(0, Nseq):
    tempList[i][0] = names[i] #creating a column with the names of the sequences
    tempList[i][1] = sequences[i] #creating a column with the sequence string

df = pd.DataFrame(tempList)
columns=list()
columns.append('name')
columns.append('sequence')
df.columns=columns

# Translate nucleic acid into amino acid
print('Translation into amino acid ...') 
df['AA']='nan'
df['sequence'].replace('-','n', regex=True, inplace=True)

for i in range (0, df.shape[0]):
    an_dna = Seq(df.iloc[i, 1])   #seq col
    aa_dna = an_dna.translate()
    df.iloc[i, 2] = str(aa_dna)   #aa col

# Checking AA column to find target mutation 
print('Checking target mutation...')
MutPos = {47:'Y', 85:'L', 91:'A', 94:'F', 104:'S', 108:'I', 137:'T', 216:['N','P']} #MutPos -1 car i commence à 0
MutPosTOcol = {'Y':'F47Y','L':'S85L','A':'V91A','F':'L94F','S':'P104S','I':'T108I','T':'A137T','N':'T216N',
               'P':'T216P'}
for f in MutPosTOcol : 
    df[MutPosTOcol[f]]=0

for i in range (0, df.shape[0]) :
    aa = df.loc[i,'AA']
    for k in MutPos :
        if aa[k-1]==MutPos[k]:
            for f in MutPosTOcol :
                if MutPos[k] == f :
                    df.loc[i,MutPosTOcol[f]]=1
                else : 
                    continue 
        else : 
            continue
 
for value in MutPosTOcol.values() : 
    try :
        print('There are '+str(df[value].value_counts()[1])+' '+str(value)+' mutation.')
    except KeyError :
        print('There is no '+str(value)+' mutation.')

# Extract percentage
print('Generating percentage for each mutation...')
table = df[['F47Y', 'S85L', 'V91A', 'L94F', 'P104S','T108I', 'A137T', 'T216N', 'T216P']]
tableT = table.T
l = tableT.shape[1]
tableT['sum'] = tableT.sum(axis=1)
tableT['frequency (%)'] = (tableT['sum']/l)*100
tableT.reset_index(inplace=True)
tableT.rename(columns={'index':'Mutation'}, inplace = True)
cols = list(range(1, Nseq+2))
tableT.drop(tableT.columns[cols],axis=1,inplace=True)
print(tableT)

# Saving csv tables

nametbl = nameFastaFile +'_mutations_percentage.csv'
print ('Creating '+ nametbl+ ' ...')
tableT.to_csv(nametbl, index = False)
print('done !')