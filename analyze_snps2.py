# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 10:01:33 2017

@author: kuns
"""

from __future__ import division
import numpy as np
from snps import GenomeSNPs
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import os
import subprocess
from fnmatch import filter
from natsort import natsorted
#from Bio import SeqIO
from collections import Counter
import counter_stats as cs


#plt.ion()

min_freq = 0.2
min_cov = 20

savefigs = True

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
#
#
def get_sample_cdf(samples):
    dtype = 'tot'
    snp_nums = {chrom: Counter() for chrom in ['mt', 'nuc']}
    cov = {chrom: 0 for chrom in ['mt', 'nuc']}
    for sample in samples.values():
        coverage, snp_freq = sample.mt.get_data(dtype)
        cov['mt'] += np.count_nonzero(coverage >= min_cov)
        snps = sample.mt.get_snps(min_freq, min_cov, dtype)
        inds = np.nonzero(snps)[0]
        snp_nums['mt'] += Counter(snps[inds])
        
        coverage, snp_freq = sample.nuc.get_data(dtype)
        cov['nuc'] += np.count_nonzero(coverage >= min_cov)
        snps = sample.nuc.get_snps(min_freq, min_cov, dtype)
        inds = np.nonzero(snps)[0]
        snp_nums['nuc'] += Counter(snps[inds])
    x_mit, cdf_mit = cs.cdf(snp_nums['mt'])
    x, frac_mit = cs.cdf(snp_nums['mt'], norm=False)
    x_nuc, cdf_nuc = cs.cdf(snp_nums['nuc'])
    x, frac_nuc = cs.cdf(snp_nums['nuc'], norm=False)
    frac_mit = frac_mit / cov['mt']
    frac_nuc = frac_nuc / cov['nuc']
    mito = {'x': x_mit, 'cdf': cdf_mit, 'frac': frac_mit}
    nuc = {'x': x_nuc, 'cdf': cdf_nuc, 'frac': frac_nuc}
    return mito, nuc


def get_phenotype_cdf(samples):
    sample_names = samples.keys()
    mut = filter(sample_names, '1914E*')
    wt = [name for name in sample_names if name not in mut]
    m, n = get_sample_cdf({k: samples[k] for k in mut})
    mito = {'mut': m}
    nuc = {'mut': n}
    m, n = get_sample_cdf({k: samples[k] for k in wt})
    mito['wt'] = m
    nuc['wt'] = n
    return mito, nuc   
    

mito_L, nuc_L = get_phenotype_cdf(samples_L)
mito_M, nuc_M = get_phenotype_cdf(samples_M)
mito_S, nuc_S = get_phenotype_cdf(samples_S)



fig1 = plt.figure()
fig1.gca().plot(mito_L['wt']['x'], mito_L['wt']['cdf'], 'C0-', label='large')
fig1.gca().plot(mito_L['mut']['x'], mito_L['mut']['cdf'], 'C0--')
#fig1.gca().plot(nuc_L['wt']['x'], nuc_L['wt']['cdf'], 'C0-', alpha=0.5)
#fig1.gca().plot(nuc_L['mut']['x'], nuc_L['mut']['cdf'], 'C0--', alpha=0.5)
fig1.gca().plot(mito_M['wt']['x'], mito_M['wt']['cdf'], 'C1-', label='medium')
fig1.gca().plot(mito_M['mut']['x'], mito_M['mut']['cdf'], 'C1--')
fig1.gca().plot(mito_S['wt']['x'], mito_S['wt']['cdf'], 'C2-', label='small')
fig1.gca().plot(mito_S['mut']['x'], mito_S['mut']['cdf'], 'C2--')
title = r'CDF of SNPs with '
title += r'$\nu_{{\mathrm{{min}}}} = {:0.2f}$,'.format(min_freq)
title += r'  $c_{{\mathrm{{min}}}} = {:d}$'.format(min_cov)
fig1.gca().set_title(title)
fig1.gca().set_xlabel(r'$\nu_{\mathrm{SNP}}$')
handles, labels = fig1.gca().get_legend_handles_labels()
line1 = mlines.Line2D([], [], color='k', linestyle='-')
line2 = mlines.Line2D([], [], color='k', linestyle='--')
handles.extend([line1, line2])
labels.extend(['wt', 'mut'])
fig1.gca().legend(handles, labels)


#fig2 = plt.figure()
#fig2.gca().plot(mito_L['wt']['x'], mito_L['wt']['frac'], 'C0-', label='large')
#fig2.gca().plot(mito_L['mut']['x'], mito_L['mut']['frac'], 'C0--')
#fig2.gca().plot(mito_M['wt']['x'], mito_M['wt']['frac'], 'C1-', label='medium')
#fig2.gca().plot(mito_M['mut']['x'], mito_M['mut']['frac'], 'C1--')
#fig2.gca().plot(mito_S['wt']['x'], mito_S['wt']['frac'], 'C2-', label='small')
#fig2.gca().plot(mito_S['mut']['x'], mito_S['mut']['frac'], 'C2--')
#title = r'Fraction of SNPs with '
#title += r'$\nu_{{\mathrm{{min}}}} = {:0.2f}$,'.format(min_freq)
#title += r'  $c_{{\mathrm{{min}}}} = {:d}$'.format(min_cov)
#fig2.gca().set_title(title)
#fig2.gca().set_xlabel(r'$\nu_{\mathrm{SNP}}$')
#handles, labels = fig2.gca().get_legend_handles_labels()
#line1 = mlines.Line2D([], [], color='k', linestyle='-')
#line2 = mlines.Line2D([], [], color='k', linestyle='--')
#handles.extend([line1, line2])
#labels.extend(['wt', 'mut'])
#fig2.gca().legend(handles, labels)




if savefigs:
    fig1.savefig(save_dir + 'snps_induced_' + str(min_cov) + '.pdf')
#    fnames = []
#    for fi, fig in enumerate([fig1, fig2]):
#        fname = save_dir + str(fi) + '.pdf'
#        fig.savefig(fname)
#        fnames.append(fname)
#    cmd = ['pdftk']
#    cmd.extend(fnames)
#    cmd.extend(['output', save_dir + 'snps_induced_' + str(min_cov) + '.pdf'])
#    status = subprocess.call(cmd)
#    if status == 0:
#        for f in fnames:
#            os.remove(f)
