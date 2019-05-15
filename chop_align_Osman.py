#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 15:16:32 2019

@author: lyan
"""


from __future__ import division
import os
from joblib import Parallel, delayed
#import subprocess
from fnmatch import filter
from natsort import natsorted
import gzip
from Bio import SeqIO
import numpy as np
import pysam
#import matplotlib.pyplot as plt

data_dir = '/data/lyan/yeast/Illumina/Osman/reads/' #Illumina/GlycerolSelection/reads/' #/nanopore/Dec17/workspace/pass/' #MitoSNPs
save_dir = '/data/lyan/yeast/Illumina/Osman/reads/'
alig_dir = '/data/lyan/yeast/Illumina/Osman/alignments/'
ref_dir  = '/data/lyan/yeast/BWAIndex/'
bwa_dir  = '/data/lyan/tools/bwa/'
fix_dir  = '/data/lyan/tools/bbmap/'
minimap_dir = '/data/lyan/tools/minimap2/'
dtype    = 'fastq'


min_qual = 20
pc_len  = 70
inc_len = 10
    
def alignment(fname,reference='genome.fa',method='bwa'):
    read_dir = data_dir+fname+'/'
    files = os.listdir(read_dir)
    ff = set()
    tt = []
    for file in files:
        ff.add(file.rsplit('_',2)[0])
        tt.append(file.rsplit('_',2)[-1])
    ff  = list(ff)
    if method=='bwa':
        # -x ont2d
        # -p for paired end reads
        c_align = bwa_dir+'bwa mem -t 16 -M '+ref_dir+reference+' '+read_dir+ff[0]+'_R1_'+tt[0]+' '+read_dir+ff[0]+'_R2_'+tt[1]+' > '+alig_dir+fname+'_paired_'+reference.split('.')[0]+'.sam'
    elif method=='minimap':
        # ./minimap2 -ax map-ont ref.fa ont.fq.gz > aln.sam         # Oxford Nanopore genomic reads
        # ./minimap2 -x ava-ont reads.fa reads.fa > overlaps.paf    # Nanopore read overlap
        # Overlap for PacBio reads (or use "-x ava-ont" for nanopore read overlapping)
        # minimap2/minimap2 -x ava-pb -t8 pb-reads.fq pb-reads.fq | gzip -1 > reads.paf.gz
        # Layout
        # miniasm/miniasm -f reads.fq reads.paf.gz > reads.gfa
        # extract fasta from gfa
        # awk '/^S/{print ">"$2"\n"$3}' in.gfa | fold > out.fa
        c_align = minimap_dir+'minimap2 -ax map-ont '+ref_dir+reference+' '+save_dir+fname+'.'+dtype+' > '+alig_dir+fname+'_'+reference.split('.')[0]+'.sam'
    print( 'aligning %s' % fname)
    os.system(c_align)
    fname = fname+'_paired_'+reference.split('.')[0]
    #fname=fname+'_1'
    c_tobam = 'samtools view -b -S -o '+alig_dir+fname+'.bam '+alig_dir+fname+'.sam'
    print( 'to bamfiles %s' % fname)
    os.system(c_tobam)
    # use samtools to sort the aligned reads
    c_sort = 'samtools sort -T /tmp/lyan.sorted -o '+alig_dir+fname+'.sorted.bam '+alig_dir+fname+'.bam'
    os.system(c_sort)
    # use samtools to index the sorted reads
    c_index = 'samtools index '+alig_dir+fname+'.sorted.bam'
    os.system(c_index)
    os.remove(alig_dir+fname+'.sam')
    os.remove(alig_dir+fname+'.bam')
    print( 'Done alignment %s' % fname)

samples = [fname.split('.')[0] for fname in os.listdir(data_dir)]

for sample in samples:  #name in names:  #
    alignment(sample, reference='mitodnalacostrain.fasta')
    #alignment(sample,reference='genome.fa')
    #alignment(name,reference='W303_chrM_mito1.fasta',method='minimap')
    #alignment(name,reference='Mitochondria2.fa',method='minimap')
    
print('Done alignment')



#samples_L = {sample: GenomeSNPs(sample, data_dir)
#             for sample in names}