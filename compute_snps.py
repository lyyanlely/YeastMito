# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 18:56:56 2017

@author: kuns
"""

from snps import analyze_sample
from joblib import Parallel, delayed
import os
from fnmatch import filter
from Bio import SeqIO
import numpy as np

cutoffAT = 0.7
min_qual = 30

n_jobs = 14
win_len = 200
win = np.ones(win_len) / win_len

data_dir = '/data/kuns/kuns/MitoMutants/alignments/'
save_dir = 'snps/'
record_dir = '/data/kuns/kuns/GenBank/'
ref_dir = '/data/kuns/kuns/ReferenceGenome/Chromosomes/'

chrV_features = SeqIO.read(record_dir + 'chrV.gb', 'genbank').features
chrM_features = SeqIO.read(record_dir + 'chrM.gb', 'genbank').features

def getSeq(chrom):
    seq = str(SeqIO.read(ref_dir + chrom + '.fa', 'fasta').seq)
    return np.fromstring(seq, dtype='S1')

def getAT(seq):
    at = np.logical_or(seq == 'A', seq == 'T')
    return np.convolve(at, win, 'same')

chrV = getSeq('chrV')
chrM = getSeq('chrM')
chrV_AT = getAT(chrV)
chrM_AT = getAT(chrM)



samples = filter(os.listdir(data_dir), '*.bai')
samples = sorted([sample.split('.')[0] for sample in samples])
samples.remove('1914EL19')

#analyze_sample(samples[0], data_dir, save_dir, chrM_features, chrV_features,
#               chrM_AT, chrV_AT, cutoffAT, min_qual, True)

out = Parallel(n_jobs=n_jobs)(
    delayed(analyze_sample)(
    sample, data_dir, save_dir, chrM_features, chrV_features, chrM_AT, chrV_AT, chrM, chrV,
    cutoffAT, min_qual, False)
    for sample in samples)
    
print 'Done analyzing'
