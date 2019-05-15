# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 10:01:33 2017

@author: kuns
"""

from __future__ import division
import numpy as np
#from alleles import GenomePolymorphisms
from snps import GenomeSNPs
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy.ma as ma
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import os
import subprocess
from fnmatch import filter
from natsort import natsorted
from Bio import SeqIO
from collections import Counter
import counter_stats as cs


plt.ion()

min_freq = 0.1
min_cov = 30

record_file = 'GenBank/chrM.gb'
records = SeqIO.read(record_file, 'genbank')
features = records.features

record_file = 'GenBank/chrIV.gb'
records = SeqIO.read(record_file, 'genbank')
features_nuc = records.features

features_mut = {'rep_origin': [[80416,80687],[52556,52821],[50559,50832],[41188,43892],[31904,32174],[29967,30254],[13224,13529]], 
                'gene': [[14567,24550],[42687,42917],[25505,25651],[70119,70868],[36139,39456],[44914,46128],[26343,27122],[75764,76573],[69379,69817]],
                'exon': [[14567,14735],[17184,17219],[19738,19775],[21291,21767],[22778,23164],[24076,24550],[42687,42917],[25505,25651],[70119,70868],[36139,36894],[38312,38362],[39106,39456],[44914,46128],[26343,27122],[75764,76573],[69379,69817]],
                'intron': [[14736,17183],[17220,19737],[19776,21290],[21768,22777],[23165,24075],[36895,38311],[38363,39105]],
                'rRNA': [[6517,8203],[53987,56695],[57838,58418],[83129,83557]],
                'rRNAintron': [[56696,57837]],
                'tRNA': [[63107,63179],[63265,63339],[66008,66083],[64107,64181],[65692,65767],[67305,67375],[68488,68563],[631,702],[82817,82894],[61859,61943],[60273,60348],[34970,35041],[73939,74013],[62862,62935]]}

ignore = ['AI5_ALPHA', 'AI5_BETA', 'AI4','AI3', 'AI2', 'AI1','BI4', 'BI3',
          'BI2', 'SCEI']

data_dir = 'Osman/snps/'
save_dir = 'Osman/plots/'
names = [fname.split('.')[0] for fname in os.listdir(data_dir)]
# names = natsorted([sample for sample in filter(names, '*chrM*')])

#names_L = natsorted([sample for sample in filter(samples, '*L*')])
#names_M = natsorted([sample for sample in filter(samples, '*M*')])
#names_S = natsorted([sample for sample in filter(samples, '*S*')])

samples = {sample: GenomeSNPs(sample, data_dir)
             for sample in names}


def add_fdict(fig, features=features_mut):
    for k,v in features.items():
        name = ''
        if k == 'gene' or k == 'rRNA':
            color = 'C1'
        elif k == 'rep_origin':
            color = 'C2'
        else:
            continue
        for i, t in enumerate(v):
            start = t[0]
            end = t[1]
            fig.gca().axvspan(start, end, alpha=0.25, color=color)
            gene_box = mpatches.Patch(color='C1', alpha=0.25)
            ori_box = mpatches.Patch(color='C2', alpha=0.25)
    return fig, gene_box, ori_box

def add_features(fig, features, ignore_features=ignore, isW303=False):
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
        else:
            continue
        if name in ignore_features:
            continue
        start = feature.location.start.position
        end = feature.location.end.position
        if isW303:
            if start>30000:
                start =start+1802
                end   =end+1802
        fig.gca().axvspan(start, end, alpha=0.25, color=color)
        gene_box = mpatches.Patch(color='C1', alpha=0.25)
        ori_box = mpatches.Patch(color='C2', alpha=0.25)
    return fig, gene_box, ori_box

def plot_freqs(samples,fig,save_name):
    for sname in natsorted(samples.keys()):
        sample = samples[sname]
        if 'T3' in sname:
            ccode = 'C0.'
            label = 'T3'
        if 'T6' in sname:
            ccode = 'C1.'
            label = 'T6'
        if 'T7' in sname:
            ccode = 'C2.'
            label = 'T7'
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
#        

def plot_samples_ref(samples, save_name, ref_locus, dtype='tot'):
    # ref_freqs, ref_locus = sample_ref.mt.get_snps(min_freq, min_cov, dtype)
    def plot_stats(sample, min_freq, min_cov, features=None, ignore=None, ref_locus=None, dtype='tot', ltype = 'cov'):
        snp_freqs, snp_locus = sample.get_snps(min_freq, min_cov, dtype)
        # (del_freqs, del_locus) , (ins_freqs, ins_locus) = sample.get_indels(min_freq, min_cov, dtype)
        # snp_freqs = ma.masked_where(snp_freqs == 0, snp_freqs)
        # del_freqs = ma.masked_where(del_freqs == 0, del_freqs)
        # ins_freqs = ma.masked_where(ins_freqs == 0, ins_freqs)
        # ref_cov   = np.ones(self.chrom_len, dtype=float)
        fig = plt.figure()
        # if ltype == 'cov': 
        #     win = np.ones(1000)/1000
        #     cov = np.convolve(self.coverage[dtype], win, 'same')
        #     max_cov = np.max(cov)
        #     cov = cov / max_cov
        #     fig.gca().plot(np.array([0,self.chrom_len]),np.array([1,1])*min_cov/max_cov, 'C4--', label=r'$c_{\mathrm{min}}/c_{\mathrm{max}}$')
        # else:
        #     cov = self.freqAT
        # fig.gca().plot(cov, 'C9-', label=r'$c/c_{\mathrm{max}}$', alpha=0.5)
        fig.gca().plot(snp_locus, snp_freqs, 'C0.', label=r'$\nu_{\mathrm{SNP}}$')
        if ref_locus is not None:
            fig.gca().plot(snp_locus[np.isin(snp_locus,ref_locus)], snp_freqs[np.isin(snp_locus,ref_locus)], 'C1.', label=r'$\nu_{\mathrm{SNP far from GC}}$')
        #fig.gca().plot(del_locus, del_freqs, 'C3.', label=r'$\nu_{\mathrm{del}}$')
        #fig.gca().plot(ins_locus, ins_freqs, 'C8.', label=r'$\nu_{\mathrm{ins}}$')
        fig.gca().set_xlabel('position (kbp)')
        title = r'$\nu_{{\mathrm{{min}}}} = {:0.2f}$'.format(min_freq)
        title += r'$\qquad c_{{\mathrm{{min}}}} = {:d}$'.format(min_cov)
        fig.gca().set_title(title)
        fig.gca().set_ylim([0, 1.2])
        fig.gca().get_xaxis().set_major_formatter(
            mpl.ticker.FuncFormatter(lambda x, p: int(x*1e-3)))
        if features is not None:
            fig, gene_box, ori_box = add_fdict(fig, features=features_mut)  # add_features(fig, features, isW303=True)  #add_features(fig, features)
            handles, labels = fig.gca().get_legend_handles_labels()
            phandles = handles[1:]
            plabels = labels[1:]
            phandles.extend([handles[0], gene_box, ori_box])
            plabels.extend([labels[0], 'gene', 'rep. origin'])
            fig.legend(phandles, plabels, ncol=2, loc='best')
        else:
            handles, labels = fig.gca().get_legend_handles_labels()
            fig.legend(handles, labels, ncol=2, loc='best')
        return fig
    fnames = []
    for sname in natsorted(samples.keys()):
        sample = samples[sname]
        fig = plot_stats(sample.mt, min_freq, min_cov, features, ignore, ref_locus, 'tot')
        title = sname + '\n'
        title += r'$\nu_{{\mathrm{{min}}}} = {:0.2f}$'.format(min_freq)
        title += r'$\qquad c_{{\mathrm{{min}}}} = {:d}$'.format(min_cov)
        fig.gca().set_title(title)
        fname = save_dir + sname + '_GC.pdf'
        fig.savefig(fname)
        fnames.append(fname)
        plt.close(fig)
    #cmd = ['pdftk']
    #cmd.extend(fnames)
    #cmd.extend(['output', save_dir + save_name])
    #status = subprocess.call(cmd)
    plt.close('all')

print('plotting large')
plot_samples(samples, 'large_mut_snps.pdf')
#
#print 'plotting medium'
#plot_samples(samples_M, 'medium_snps.pdf')

#print 'plotting small'
#plot_samples(samples_S, 'small_snps.pdf')


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
    keys = samples.keys()
    keys = [name.rsplit('_',4)[0] for name in keys]
        
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
    fig1.gca().set_xticks(np.arange(si+1))
    fig1.gca().set_xticklabels(sorted(keys), rotation=45 )
    
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
    fig2.gca().set_xticks(np.arange(si+1))
    fig2.gca().set_xticklabels(sorted(keys), rotation=45 )
    return fig1, fig2
#
#

num_L, frac_L = plot_sample_stats(samples, 'Large')
#num_M, frac_M = plot_sample_stats(samples_M, 'Medium')
#num_S, frac_S = plot_sample_stats(samples_S, 'Small')
#
fL = save_dir + 'num_tot.pdf'
#fM = save_dir + 'fracM_code.pdf'
#fS = save_dir + 'fracS_code.pdf'
#
num_L.savefig(fL)
#num_M.savefig(fM)
#num_S.savefig(fS)
#
#cmd = ['pdftk']
#cmd.extend(fL) #[fL, fM, fS])
#cmd.extend(['output', save_dir + 'num_snps.pdf'])
#status = subprocess.call(cmd)
#if status == 0:
#    for fn in fL: #[fL, fM, fS]:
#        os.remove(fn)
#        
fL = save_dir + 'fracL_tot.pdf'
frac_L.savefig(fL)
#frac_M.savefig(fM)
#frac_S.savefig(fS)

#cmd = ['pdftk']
#cmd.extend(fL) #[fL, fM, fS])
#cmd.extend(['output', save_dir + 'frac_snps.pdf'])
#status = subprocess.call(cmd)
#if status == 0:
#    for fn in fL: #[fL, fM, fS]:
#        os.remove(fn)
plt.close('all')


def get_sample_cdf(samples):
    snp_nums = {chrom: Counter() for chrom in ['mt', 'nuc']}
    cov = {chrom: 0 for chrom in ['mt', 'nuc']}
    for sample in samples.values():
        cov['mt'] += np.count_nonzero(sample.mt.coverage['tot'] >= min_cov)
        snps = sample.mt.get_snps(min_freq, min_cov)
        #inds = np.nonzero(snps)[0]
        snp_nums['mt'] += Counter(snps[0])
        
        cov['nuc'] += np.count_nonzero(sample.nuc.coverage['tot'] >= min_cov)
        snps = sample.nuc.get_snps(min_freq, min_cov)
        #inds = np.nonzero(snps)[0]
        snp_nums['nuc'] += Counter(snps[0])
    x_mit, cdf_mit = cs.cdf(snp_nums['mt'])
    x, frac_mit = cs.cdf(snp_nums['mt'], norm=False)
    x_nuc, cdf_nuc = cs.cdf(snp_nums['nuc'])
    x, frac_nuc = cs.cdf(snp_nums['nuc'], norm=False)
    frac_mit = frac_mit / cov['mt']
    frac_nuc = frac_nuc / cov['nuc']
    mito = {'x': x_mit, 'cdf': cdf_mit, 'frac': frac_mit}
    nuc = {'x': x_nuc, 'cdf': cdf_nuc, 'frac': frac_nuc}
    return mito, nuc

mito_L, nuc_L = get_sample_cdf(samples)
fig1 = plt.figure()
fig1.gca().plot(mito_L['x'], mito_L['cdf'], 'C0-', label='mitochondrial')
fig1.gca().plot(nuc_L['x'], nuc_L['cdf'], 'C0--', label='nuclear')
fig1.gca().set_xlabel(r'$\nu_{\mathrm{SNP}}$')
fig1.gca().legend()
fig1.savefig(save_dir+'snps_freq.pdf')
plt.close(fig1)

"""
trimer_nuc = {}
for i in np.arange(nuc_len-2):
    trimer = nuc_seq[i:i+3]
    if trimer not in trimer_nuc:
        trimer_nuc.update({trimer:[i+1]})
    else:
        trimer_nuc[trimer].append(i+1)

trimer_mt = {}
trimer = mt_seq[-1]+mt_seq[0:2]
if trimer not in trimer_mt:
    trimer_mt.update({trimer:[0]})
else:
    trimer_mt[trimer].append(0)
for i in np.arange(mt_len-2):
    trimer = mt_seq[i:i+3]
    if trimer not in trimer_mt:
        trimer_mt.update({trimer:[i+1]})
    else:
        trimer_mt[trimer].append(i+1)
trimer = mt_seq[-2:]+mt_seq[0]
if trimer not in trimer_mt:
    trimer_mt.update({trimer:[mt_len-1]})
else:
    trimer_mt[trimer].append(mt_len-1)


for i, k in enumerate(np.sort(trimer_nuc.keys())):
    trimer_count[i] = len(trimer_nuc[k])

for ri, key in enumerate(samples.keys()):
    name = key.rsplit('_',4)[0]
    sample = samples[key]
    inds = np.logical_and(sample.mt.snp_freq['tot'][:,0] >= min_freq, sample.mt.coverage['tot'] >= min_cov)
    lc_mt = np.argwhere(inds)
    inds = np.logical_and(sample.nuc.snp_freq['tot'][:,0] >= min_freq, sample.nuc.coverage['tot'] >= min_cov)
    lc_nuc = np.argwhere(inds)
    mt_snpfrac = np.zeros(64,)
    nuc_snpfrac = np.zeros(64,)
    for i, k in enumerate(np.sort(trimer_mt.keys())):
        mt_snpfrac[i]  = np.mean(np.isin(trimer_mt[k],lc_mt))
        nuc_snpfrac[i] = np.mean(np.isin(trimer_nuc[k],lc_nuc))

    fig = plt.figure(figsize=(9,8))
    fig.gca().bar(np.arange(64,), mt_snpfrac)
    fig.gca().set_xticks(np.arange(64,))
    fig.gca().set_xticklabels(np.sort(trimer_mt.keys()),rotation=90)
    fig.gca().set_ylabel('snp fraction')
    fig.gca().set_title(name)
    plt.savefig('mt_trimer_snpfrac_'+name+'_LacO_chrM.pdf')
    plt.close(fig)

    fig = plt.figure(figsize=(9,8))
    fig.gca().bar(np.arange(64,), nuc_snpfrac)
    fig.gca().set_xticks(np.arange(64,))
    fig.gca().set_xticklabels(np.sort(trimer_nuc.keys()),rotation=90)
    fig.gca().set_ylabel('snp fraction')
    fig.gca().set_title(name)
    plt.savefig('nuc_trimer_snpfrac_'+name+'_S288c_chrIV.pdf')
    plt.close(fig)

mito_M, nuc_M = get_sample_cdf(samples_M)
mito_S, nuc_S = get_sample_cdf(samples_S)


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