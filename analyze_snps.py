# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 10:01:33 2017

@author: kuns lyan
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
cutoffAT = 0.7

savefigs = True

data_dir = 'snps/'
save_dir = 'plots/'
samples = [fname.split('.')[0] for fname in os.listdir(data_dir)]

names_L = natsorted([sample for sample in filter(samples, '*L*')])
names_M = natsorted([sample for sample in filter(samples, '*M*')])
names_S = natsorted([sample for sample in filter(samples, '*S*')])

names_LMS = names_L+names_M+names_S

names_14E = natsorted([sample for sample in list(set(filter(samples, '1914E*'))-set(names_LMS))])
names_14 = natsorted([sample for sample in list(set(filter(samples, '1914*'))-set(names_LMS)-set(names_14E))])
names_08 = natsorted([sample for sample in list(set(filter(samples, '1908*'))-set(names_LMS))])
names_03 = natsorted([sample for sample in list(set(filter(samples, 'W303*'))-set(names_LMS))])

samples_14E = {sample: GenomeSNPs(sample, data_dir)
             for sample in names_14E}
samples_14 = {sample: GenomeSNPs(sample, data_dir)
             for sample in names_14}
samples_08 = {sample: GenomeSNPs(sample, data_dir)
             for sample in names_08}
samples_03 = {sample: GenomeSNPs(sample, data_dir)
             for sample in names_03}


def get_sample_cdf(samples, dtype='tot'):
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


def plot_all_samples(samples_14E, samples_14, samples_08, samples_03, category):
    fig_n = plt.figure()
    fig_f = plt.figure()
    plot_samples(samples_14E, category, fig_n, fig_f, '1914E', 'C0')
    plot_samples(samples_14, category, fig_n, fig_f, '1914', 'C1')
    plot_samples(samples_08, category, fig_n, fig_f, '1908', 'C2')
    plot_samples(samples_03, category, fig_n, fig_f, 'W303', 'C3')
    return fig_n, fig_f
    


def plot_samples(samples, category, fig_n, fig_f, size, color):
    plot_type(samples, 'cdf', 'tot', fig_n, color + '-', label=size)
    plot_type(samples, 'frac', 'tot', fig_f, color + '-', label=size)
    if category == 'AT':
        plot_type(samples, 'cdf', 'high', fig_n, color + '--')
        plot_type(samples, 'cdf', 'low', fig_n, color + ':')
        plot_type(samples, 'frac', 'high', fig_f, color + '--')
        plot_type(samples, 'frac', 'low', fig_f, color + ':')
    elif category == 'coding':
        plot_type(samples, 'cdf', 'non', fig_n, color + '--')
        plot_type(samples, 'cdf', 'code', fig_n, color + ':')
        plot_type(samples, 'frac', 'non', fig_f, color + '--')
        plot_type(samples, 'frac', 'code', fig_f, color + ':')


def plot_type(samples, cdf_type, dtype, fig, sty, label=None):
    mito, nuc = get_sample_cdf(samples, dtype)
    fig.gca().plot(mito['x'], mito[cdf_type], sty, label=label)
    fig.gca().plot(nuc['x'], nuc[cdf_type], sty, alpha=0.5)


def plot_labels(fig, cdf_type, dtype):
    if cdf_type == 'frac':
        fig.gca().set_yscale('log')
        title = 'Fraction of SNPs with '
    elif cdf_type == 'cdf':
        title = 'CDF of SNPs with '
    if dtype == 'AT':
        label1 = r'$\nu_{{\mathrm{{AT}}}} > {:0.2f}$'.format(cutoffAT)
        label2 = r'$\nu_{{\mathrm{{AT}}}} < {:0.2f}$'.format(cutoffAT)
    elif dtype == 'coding':
        label1 = 'coding'
        label2 = 'non-coding'
    fig.gca().set_xlabel(r'$\nu_{\mathrm{SNP}}$')
    title += r'$\nu_{{\mathrm{{min}}}} = {:0.2f}$,'.format(min_freq)
    title += r'  $c_{{\mathrm{{min}}}} = {:d}$'.format(min_cov)
    fig.gca().set_title(title)
    handles, labels = fig.gca().get_legend_handles_labels()
    line_solid = mlines.Line2D([], [], color='k', linestyle='-')
    line1 = mlines.Line2D([], [], color='k', linestyle='--')
    line2 = mlines.Line2D([], [], color='k', linestyle=':')
    line3 = mlines.Line2D([], [], color='k', alpha=0.5)
    handles.extend([line_solid, line1, line2, line_solid, line3])
    labels.extend(['tot', label1, label2, 'mit', 'nuc'])
    fig.gca().legend(handles, labels, ncol=2)
    
    
"""
mito_L, nuc_L = get_sample_cdf(samples_L)
mito_M, nuc_M = get_sample_cdf(samples_M)
mito_S, nuc_S = get_sample_cdf(samples_S)
"""

mito_14E, nuc_14E = get_sample_cdf(samples_14E)
mito_14, nuc_14 = get_sample_cdf(samples_14)
mito_08, nuc_08 = get_sample_cdf(samples_08)
mito_03, nuc_03 = get_sample_cdf(samples_03)

figAT_n, figAT_f = plot_all_samples(samples_14E, samples_14, samples_08, samples_03, 'AT')
figCode_n, figCode_f = plot_all_samples(samples_14E, samples_14, samples_08, samples_03,
                                        'coding')

plot_labels(figAT_n, 'cdf', 'AT')
plot_labels(figAT_f, 'frac', 'AT')
plot_labels(figCode_n, 'cdf', 'coding')
plot_labels(figCode_f, 'frac', 'coding')


if savefigs:
    fnames = []
    for fi, fig in enumerate([figAT_n, figAT_f, figCode_n, figCode_f]):
        fname = save_dir + str(fi) + '.pdf'
        fig.savefig(fname)
        fnames.append(fname)
    cmd = ['pdftk']
    cmd.extend(fnames)
    cmd.extend(['output', save_dir + 'snps_1_' + str(min_cov) + '.pdf'])
    status = subprocess.call(cmd)
    if status == 0:
        for f in fnames:
            os.remove(f)
