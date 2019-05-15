# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 18:56:56 2017

@author: kuns
"""

from alleles import analyze_sample
from joblib import Parallel, delayed
import os
from fnmatch import filter


n_jobs = 14

data_dir = '/data/kuns/kuns/MitoMutants/alignments/'
save_dir = '/data/kuns/kuns/MitoMutants/alleles/'

samples = filter(os.listdir(data_dir), '*.bai')
samples = sorted([sample.split('.')[0] for sample in samples])
samples.remove('1914EL19')


out = Parallel(n_jobs=n_jobs)(
    delayed(analyze_sample)(sample, data_dir, save_dir, 30)
    for sample in samples)
    
print 'Done analyzing'