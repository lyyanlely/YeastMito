from __future__ import division
import numpy as np
from alleles import GenomePolymorphisms
from snps import GenomeSNPs
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import os
import subprocess
from fnmatch import filter
from natsort import natsorted
from Bio import SeqIO
from collections import Counter
import counter_stats as cs

data_dir = 'snps/'
save_dir = 'plots/'

samples = [fname.split('.')[0] for fname in os.listdir(data_dir)]

names_L = natsorted([sample for sample in filter(samples, '*L*')])
names_M = natsorted([sample for sample in filter(samples, '*M*')])
names_S = natsorted([sample for sample in filter(samples, '*S*')])

samples_L = {sample: GenomeSNPs(sample, data_dir)
             for sample in names_L}
samples_M = {sample: GenomeSNPs(sample, data_dir)
             for sample in names_M}
samples_S = {sample: GenomeSNPs(sample, data_dir)
             for sample in names_S}