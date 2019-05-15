# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 10:01:33 2017

@author: kuns
"""

from __future__ import division
import numpy as np
from alleles import GenomePolymorphisms
from snps import GenomeSNPs
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import os
import subprocess
from fnmatch import filter
from natsort import natsorted
from Bio import SeqIO
from collections import Counter
import counter_stats as cs


#plt.ion()

min_freq = 0.2
min_cov = 20
cutoffAT = 0.7

savefigs = True

data_dir = 'snps/'
save_dir = 'plots/'

record_file = 'GenBank/chrM.gb'
records = SeqIO.read(record_file, 'genbank')
features = records.features

ignore = ['AI5_ALPHA', 'AI5_BETA', 'AI4','AI3', 'AI2', 'AI1','BI4', 'BI3',
          'BI2', 'SCEI']

"""
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

samples = {}
samples.update(samples_14E)
samples.update(samples_14)
samples.update(samples_08)
samples.update(samples_03)
"""

def add_features(fig, features=features, ignore_features=ignore): 
    # add genetic features to a gene plot
    for feature in features:
        name = ''
        if feature.type == 'gene':
            color = 'C1'
            try:
                name = feature.qualifiers['gene'][0]
            except KeyError:
                continue
        elif feature.type == 'rep_origin':
            color = 'C2'
            try:
                name = feature.qualifiers['note'][0]
            except KeyError:
                pass
        elif feature.type == 'CDS':
            color = 'C3'
            try:
                name = feature.qualifiers['gene'][0]
            except KeyError:
                continue
        else:
            continue
        if name in ignore_features:
            continue
        if len(feature.location.parts)>1:
            for part in feature.location.parts:
                start = part.start.position
                end  = part.end.position
                fig.gca().axvspan(start, end, alpha=0.25, color=color)
        else:
            start = feature.location.start.position
            end = feature.location.end.position
            fig.gca().axvspan(start, end, alpha=0.25, color=color)
    gene_box = mpatches.Patch(color='C1', alpha=0.25)
    ori_box = mpatches.Patch(color='C2', alpha=0.25)        
    return fig, gene_box, ori_box

def plot_freqs(samples,fig,save_name):
    for sname in natsorted(samples.keys()):
        sample = samples[sname]
        if '08S' in sname:
            ccode = 'C0.'
            label = '1908S'
        if '08L' in sname:
            ccode = 'C1.'
            label = '1908L'
        if '14S' in sname:
            ccode = 'C2.'
            label = '1914S'
        if '14L' in sname:
            ccode = 'C3.'
            label = '1914L'
        if '14ES' in sname:
            ccode = 'C4.'
            label = '1914ES'
        if '14EM' in sname:
            ccode = 'C5.'
            label = '1914EM'
        if '14EL' in sname:
            ccode = 'C6.'
            label = '1914EL'
        if '08' in sname:
            ccode = 'C2.'
            label = '1908'
        if '14' in sname:
            ccode = 'C0.'
            label = '1914'
        if '14E' in sname:
            ccode = 'C1.'
            label = '1914E'
        if '03' in sname:
            ccode = 'C3.'
            label = 'W303'
        freqs = sample.mt.get_snps(min_freq, min_cov, 'tot')
        freqAT = sample.mt.freqAT
        freqAT = freqAT[freqs>0]
        freqs  = freqs[freqs>0]
        fig.gca().plot(freqAT, freqs, ccode, label=label)

def plot_samples(samples, save_name):
    fnames = []
    for sname in natsorted(samples.keys()):
        sample = samples[sname]
        fig = sample.mt.plot_stats(min_freq, min_cov, features, ignore, 'tot')
        title = sname + '\n'
        title += r'$\nu_{{\mathrm{{min}}}} = {:0.2f}$'.format(min_freq)
        title += r'$\qquad c_{{\mathrm{{min}}}} = {:d}$'.format(min_cov)
        fig.gca().set_title(title)
        fname = save_dir + sname + '.pdf'
        fig.savefig(fname)
        fnames.append(fname)
        plt.close(fig)
    #cmd = ['pdftk']
    #cmd.extend(fnames)
    #cmd.extend(['output', save_dir + save_name])
    #status = subprocess.call(cmd)
    plt.close('all')
    #if status == 0:
    #    for fname in fnames:
    #        os.remove(fname)
"""
#  
print 'plotting 1914'
plot_samples(samples_14, '1914_snps.pdf')
#
print 'plotting 1914E'
plot_samples(samples_14E, '1914E_snps.pdf')
#
print 'plotting 1908'
plot_samples(samples_08, '1908_snps.pdf')
#
print 'plotting W303'
plot_samples(samples_03, 'W303_snps.pdf')
#
#print 'plotting medium'
#plot_samples(samples_M, 'medium_snps.pdf')

#print 'plotting small'
#plot_samples(samples_S, 'small_snps.pdf')

"""
def plot_sample_stats(samples, size, dtypes=['tot']):
    nums = {chrom: {dtype: np.zeros(len(samples.keys())) for dtype in dtypes} for chrom in ['mt', 'nuc']}
    freqs = {chrom: {dtype: np.zeros(len(samples.keys())) for dtype in dtypes} for chrom in ['mt', 'nuc']}
    for si, sname in enumerate(sorted(samples.keys())):
        sample = samples[sname]
        for dtype in dtypes:
            s, c = sample.mt.snp_stats(min_freq, min_cov, dtype)
            nums['mt'][dtype][si] = s
            freqs['mt'][dtype][si] = s/(c + 1e-6)

            s, c = sample.nuc.snp_stats(min_freq, min_cov, dtype)
            nums['nuc'][dtype][si] = s
            freqs['nuc'][dtype][si] = s/(c + 1e-6)
        
    fig1 = plt.figure()
    fig1.gca().plot(nums['mt']['tot'], 'C0.', label='mt')
    fig1.gca().plot(nums['nuc']['tot'], 'C1.', label='nuc')
    fig1.gca().legend()
    for dtype in dtypes:
        if dtype == 'code':
            fig1.gca().plot(nums['mt'][dtype], 'C0+')
            fig1.gca().plot(nums['nuc'][dtype], 'C1+')
        if dtype == 'non':
            fig1.gca().plot(nums['mt'][dtype], 'C0x')
            fig1.gca().plot(nums['nuc'][dtype], 'C1x')
    title = size + ' Samples\n'
    title += 'Number of SNPS with   '
    title += r'$\nu_{{\mathrm{{min}}}} = {:0.2f}$'.format(min_freq)
    title += r'$\qquad c_{{\mathrm{{min}}}} = {:d}$'.format(min_cov)
    fig1.gca().set_title(title)
    
    fig2 = plt.figure()
    fig2.gca().semilogy(freqs['mt']['tot'], 'C0.', label='mt')
    fig2.gca().semilogy(freqs['nuc']['tot'], 'C1.', label='nuc')
    fig2.gca().legend()
    for dtype in dtypes:
        if dtype == 'code':
            fig2.gca().plot(freqs['mt'][dtype], 'C0+')
            fig2.gca().plot(freqs['nuc'][dtype], 'C1+')
        if dtype == 'non':
            fig2.gca().plot(freqs['mt'][dtype], 'C0x')
            fig2.gca().plot(freqs['nuc'][dtype], 'C1x')
    title = size + ' Samples\n'
    title += 'Fraction of SNPs with   '
    title += r'$\nu_{{\mathrm{{min}}}} = {:0.2f}$'.format(min_freq)
    title += r'$\qquad c_{{\mathrm{{min}}}} = {:d}$'.format(min_cov)
    fig2.gca().set_title(title)
    return fig1, fig2
#
"""
num_14E, frac_14E = plot_sample_stats(samples_14E, '1914E')
num_14, frac_14 = plot_sample_stats(samples_14, '1914')
num_08, frac_08 = plot_sample_stats(samples_08, '1908')
num_03, frac_03 = plot_sample_stats(samples_03, 'W303')
#
f14E = save_dir + 'num_1914E.pdf'
f14 = save_dir + 'num_1914.pdf'
f08 = save_dir + 'num_1908.pdf'
f03 = save_dir + 'num_W303.pdf'
#
num_14E.savefig(f14E)
num_14.savefig(f14)
num_08.savefig(f08)
num_03.savefig(f03)

f14E = save_dir + 'frac_1914E.pdf'
f14 = save_dir + 'frac_1914.pdf'
f08 = save_dir + 'frac_1908.pdf'
f03 = save_dir + 'frac_W303.pdf'

frac_14E.savefig(f14E)
frac_14.savefig(f14)
frac_08.savefig(f08)
frac_03.savefig(f03)
#

num_L, frac_L = plot_sample_stats(samples_L, 'Large')
num_M, frac_M = plot_sample_stats(samples_M, 'Medium')
num_S, frac_S = plot_sample_stats(samples_S, 'Small')
#
fL = save_dir + 'fracL_code.pdf'
fM = save_dir + 'fracM_code.pdf'
fS = save_dir + 'fracS_code.pdf'
#
num_L.savefig(fL)
num_M.savefig(fM)
num_S.savefig(fS)
#
cmd = ['pdftk']
cmd.extend([fL, fM, fS])
cmd.extend(['output', save_dir + 'num_snps.pdf'])
status = subprocess.call(cmd)
if status == 0:
    for fn in [fL, fM, fS]:
        os.remove(fn)
#        
frac_L.savefig(fL)
frac_M.savefig(fM)
frac_S.savefig(fS)

cmd = ['pdftk']
cmd.extend([fL, fM, fS])
cmd.extend(['output', save_dir + 'frac_snps.pdf'])
status = subprocess.call(cmd)
if status == 0:
    for fn in [fL, fM, fS]:
        os.remove(fn)
"""


def get_sample_cdf(samples):
    snp_nums = {chrom: Counter() for chrom in ['mt', 'nuc']}
    cov = {chrom: 0 for chrom in ['mt', 'nuc']}
    for sample in samples.values():
        cov['mt'] += np.count_nonzero(sample.mt.coverage >= min_cov)
        snps = sample.mt.get_snps(min_freq, min_cov)
        inds = np.nonzero(snps)[0]
        snp_nums['mt'] += Counter(snps[inds])
        
        cov['nuc'] += np.count_nonzero(sample.nuc.coverage >= min_cov)
        snps = sample.nuc.get_snps(min_freq, min_cov)
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



"""
mito_14E, nuc_14E = get_sample_cdf(samples_14E)
mito_M, nuc_M = get_sample_cdf(samples_M)
mito_S, nuc_S = get_sample_cdf(samples_S)


mito_14E, nuc_14E = get_sample_cdf(samples_14E)
mito_14, nuc_14 = get_sample_cdf(samples_14)
mito_08, nuc_08 = get_sample_cdf(samples_08)
mito_03, nuc_03 = get_sample_cdf(samples_03)

plot_labels(figAT_n, 'cdf', 'AT')
plot_labels(figAT_f, 'frac', 'AT')
plot_labels(figCode_n, 'cdf', 'coding')
plot_labels(figCode_f, 'frac', 'coding')

fig1 = plt.figure()
fig1.gca().plot(mito_L['x'], mito_L['cdf'], 'C0-', label='large')
fig1.gca().plot(nuc_L['x'], nuc_L['cdf'], 'C0--')
fig1.gca().plot(mito_M['x'], mito_M['cdf'], 'C1-', label='medium')
fig1.gca().plot(nuc_M['x'], nuc_M['cdf'], 'C1--')
fig1.gca().plot(mito_S['x'], mito_S['cdf'], 'C2-', label='small')
fig1.gca().plot(nuc_S['x'], nuc_S['cdf'], 'C2--')
fig1.gca().set_xlabel(r'$\nu_{\mathrm{SNP}}$')
title = 'CDF of SNPS with   '
title += r'$\nu_{{\mathrm{{min}}}} = {:0.2f}$'.format(min_freq)
title += r'$\qquad c_{{\mathrm{{max}}}} = {:d}$'.format(min_cov)
fig1.gca().set_title(title)
handles, labels = fig1.gca().get_legend_handles_labels()
mito_line = mlines.Line2D([], [], color='k', linestyle='-')
nuc_line = mlines.Line2D([], [], color='k', linestyle='--')
handles.extend([mito_line, nuc_line])
labels.extend(['mit', 'nuc'])
fig1.gca().legend(handles, labels)


fig2 = plt.figure()
fig2.gca().semilogy(mito_L['x'], mito_L['frac'], 'C0-', label='large')
fig2.gca().semilogy(nuc_L['x'], nuc_L['frac'], 'C0--')
fig2.gca().semilogy(mito_M['x'], mito_M['frac'], 'C1-', label='medium')
fig2.gca().semilogy(nuc_M['x'], nuc_M['frac'], 'C1--')
fig2.gca().semilogy(mito_S['x'], mito_S['frac'], 'C2-', label='small')
fig2.gca().semilogy(nuc_S['x'], nuc_S['frac'], 'C2--')
fig2.gca().set_xlabel(r'$\nu_{\mathrm{SNP}}$')
title = 'Fraction of SNPs with   '
title += r'$\nu_{{\mathrm{{min}}}} = {:0.2f}$'.format(min_freq)
title += r'$\qquad c_{{\mathrm{{max}}}} = {:d}$'.format(min_cov)
fig2.gca().set_title(title)
handles, labels = fig2.gca().get_legend_handles_labels()
mito_line = mlines.Line2D([], [], color='k', linestyle='-')
nuc_line = mlines.Line2D([], [], color='k', linestyle='--')
handles.extend([mito_line, nuc_line])
labels.extend(['mit', 'nuc'])
fig2.gca().legend(handles, labels)


fN = save_dir + 'num.pdf'
fF = save_dir + 'frac.pdf'
fig1.savefig(fN)
fig2.savefig(fF)
cmd = ['pdftk']
cmd.extend([fN, fF])
cmd.extend(['output', save_dir + 'snp_distributions40.pdf'])
status = subprocess.call(cmd)
if status == 0:
    for f in [fN, fF]:
        os.remove(f)


mt_locus = []
nuc_locus = []
for sname in natsorted(samples_L.keys()):
    sample = samples_S[sname]
    mt_freq, mt_loc = sample.mt.get_snps(0.2,20)
    nuc_freq,nuc_loc= sample.nuc.get_snps(0.2,20)
    mt_locus.extend(mt_loc.flatten())
    nuc_locus.extend(nuc_loc.flatten())
"""