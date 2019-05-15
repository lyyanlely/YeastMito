# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 15:36:19 2017

@author: kuns lyan
"""

from __future__ import division
import numpy as np
import pysam
from collections import Counter

short_lim = 5000
long_lim = 15000
min_qual = 20

win_len = 1000
win= np.ones(win_len) / win_len

def compute_coverage(bamname,cflag=False,cname=''):
    if cflag and len(cname)==0:
        print 'no recalibration reference!'
        return

    if isinstance(bamname, basestring):
        bamname = [bamname]

    bamfile = pysam.AlignmentFile(bamname[0], 'rb')
    cov = {bamfile.references[ri]: np.zeros(bamfile.lengths[ri])
           for ri in range(len(bamfile.references))}
    cov_short = {bamfile.references[ri]: np.zeros(bamfile.lengths[ri])
                 for ri in range(len(bamfile.references))}
    cov_long = {bamfile.references[ri]: np.zeros(bamfile.lengths[ri])
                for ri in range(len(bamfile.references))}
    cov_med = {bamfile.references[ri]: np.zeros(bamfile.lengths[ri])
               for ri in range(len(bamfile.references))}
    
    nuc_len = Counter()
    mt_len = Counter()
    nuc_qual = Counter()
    mt_qual = Counter()
    mapped = 0
    unmapped = 0
    
    for bi, bamname1 in enumerate(bamname):
        print 'Analyzing %s' % bamname1

        bamfile = pysam.AlignmentFile(bamname1, 'rb')
        
        for ri, read in enumerate(bamfile):
            if ri % 50000 == 0:
                print ri
            if read.is_unmapped:
                unmapped += 1
                continue
            mapped += 1
            chrom_name = bamfile.getrname(read.reference_id)
            if chrom_name == 'chrM':
                mt_len[read.query_length] += 1
                mt_qual[read.mapping_quality] += 1
            else:
                nuc_len[read.query_length] += 1
                nuc_qual[read.mapping_quality] += 1
            #if read.get_tag('ZF') < 0.1:
            #    continue
            if read.mapping_quality < min_qual:
                continue
            cov_pos = np.array(read.get_reference_positions())
            cov[chrom_name][cov_pos] += 1
            if read.query_length <= short_lim:
                cov_short[chrom_name][cov_pos] += 1
            elif read.query_length >= long_lim:
                cov_long[chrom_name][cov_pos] += 1
            else:
                cov_med[chrom_name][cov_pos] += 1
                
    if cflag: 
        recal = np.load(cname)
        rcov  = {}
        winsize = recal['winsize'][()]
        nbins = recal['nbins'][()]
        covrat = recal['covrat'][()]
        inds  = recal['inds'][()]
        chromAT= recal['chromAT'][()]
        for k, v in cov.items():
            rcov.update({k: v*covrat[inds[k]-1]})
            #np.append(np.array([cv*covrat[min(inds[k][ri],nbins)-1] for ri, cv in enumerate(covmat)]).flatten(),v[(nwins*winsize):]*covrat[min(indends[k][()],nbins)-1])})
    else:
        rcov  = {}
        rccov = {}

    print 'convolving coverage'
    ccov = {k: np.convolve(v, win, 'same') for k, v in cov.items()}
    if cflag:
        rccov = {k: np.convolve(v, win, 'same') for k, v in rcov.items()}
    print 'convolving long'
    ccov_long = {k: np.convolve(v, win, 'same') for k, v in cov_long.items()}
    print 'convolving medium'
    ccov_med = {k: np.convolve(v, win, 'same') for k, v in cov_med.items()}
    print 'convolving short'
    ccov_short = {k: np.convolve(v, win, 'same') for k, v in cov_short.items()}
    print 'Done analyzing'
        
    savename = bamname[0].split('.')[0] + '.npz'
    np.savez_compressed(
        savename, cov=cov, rcov=rcov, cov_short=cov_short, cov_med=cov_med,
        cov_long=cov_long, nuc_len=nuc_len, mt_len=mt_len, nuc_qual=nuc_qual,
        mt_qual=mt_qual, mapped=mapped, unmapped=unmapped, short_lim=short_lim,
        long_lim=long_lim, ccov=ccov, rccov=rccov, ccov_long=ccov_long, ccov_med=ccov_med,
        ccov_short=ccov_short, win_len=win_len)
    print 'Done saving'


# use a mapping to align reads to references
# bwa: /home/lyan/bwa/bwa mem -t 10 -M ../BWAIndex/genome.fa sra_data.fasta > alignment_desai_test.sam &
#  -x ont2d nanopore
# use samtools to compress .sam file to .bam files
# samtools view -b -S -o alignment_desai_test.bam alignment_desai_test.sam
# use samtools to sort the aligned reads
# samtools sort -T /tmp/lyan.sorted -o alignment_desai_test.sorted.bam alignment_desai_test.bam
# use samtools to index the sorted reads
# samtools index alignment_desai_test.sorted.bam
dirc = 'Desai/'
#ids = [6290,7110,7479,7283,7163,7477,7066,7209,7153,7140] #[7142,7045,6210,7403,7130,7184,7211,7087,7114,7398,7152]
#bamnames = []
#for i in ids:
#    bamnames.append(''.join([dirc,'SRR540%d.sorted.bam' % i]))
    
compute_coverage([''.join([dirc,'alignment_desai_test.sorted.bam'])])
#dirc = './'
#compute_coverage(''.join(['alignment_srr5406290.sorted.bam']),True,'Recalibration_Desai.npz')
#compute_coverage([''.join([dirc,'alignment_barcode09_0.sorted.bam']),''.join([dirc,'alignment_barcode09_1.sorted.bam'])])