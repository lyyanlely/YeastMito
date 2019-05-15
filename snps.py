# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 17:12:38 2017

@author: kuns lyan
"""

from __future__ import division
import numpy as np
import pysam
from collections import Counter
import os
import counter_stats as cs
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches
import numpy.ma as ma


covkeys = ['code', 'non', 'low', 'high', 'tot']
ignore = ['AI5_ALPHA', 'AI5_BETA', 'AI4','AI3', 'AI2', 'AI1','BI4', 'BI3',
          'BI2', 'SCEI']


def analyze_sample(sample_name, data_dir, save_dir, mtFeatures, nucFeatures,
                   mtAT, nucAT, mtref, nucref, cutoffAT=0.7,
                   min_qual=30, VERBOSE=False):
    print 'Analyzing', sample_name
    bamname = data_dir + sample_name + '.sorted.bam'
    genome = GenomeSNPs(sample_name)
    genome.mt.analyze_bam(bamname, mtFeatures, mtAT, mtref, cutoffAT=cutoffAT,
                          min_qual=min_qual, VERBOSE=VERBOSE)
    genome.nuc.analyze_bam(bamname, nucFeatures, nucAT, nucref, cutoffAT=cutoffAT,
                           min_qual=min_qual, VERBOSE=VERBOSE)
    genome.save_data(save_dir)
    print 'Done saving', sample_name


def origin_regions(features, chrom_len):
    origin = np.zeros(chrom_len, dtype=int)
    for feature in features:
        if feature.type == 'rep_origin':
            start = feature.location.start.position
            end = feature.location.end.position
            origin[start:end] = True
    return origin

def coding_regions(features, chrom_len):
    coding = np.zeros(chrom_len, dtype=int)
    for feature in features:
        if feature.type == 'gene':
            start = feature.location.start.position
            end = feature.location.end.position
            coding[start:end] = True
    return coding


def exon_regions(features, chrom_len, ignore_features=ignore):
    exon = np.zeros(chrom_len, dtype=int)
    for feature in features:
        if feature.type == 'CDS':
            try:
                name = feature.qualifiers['gene'][0]
            except KeyError:
                continue
        else:
            continue
        if name in ignore_features:
            continue
        for part in feature.location.parts:
            start = part.start.position
            end  = part.end.position
            exon[start:end] = True
    return exon

def RNA_regions(features, chrom_len, ntype, ignore_features=ignore):
    rna = np.zeros(chrom_len, dtype=int)
    for feature in features:
        if feature.type == ntype:
            try:
                name = feature.qualifiers['note'][0]
            except KeyError:
                pass
        else:
            continue
        if name in ignore_features:
            continue
        for part in feature.location.parts:
            start = part.start.position
            end  = part.end.position
            rna[start:end] = True
    return rna


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

def add_fdict(fig, features):
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

class ChromosomeSNPs:
    def __init__(self, chrom_name, data_dir=None):
        self.name = chrom_name
        self.chrom_len = None
        self.coverage = None
        self.snp_freq = None
        self.snps = None
        self.deletions = None
        self.insertions = None
        self.consensus = None
        self.min_qual = None
        self.freqAT = None
        self.cutoffAT = None
        if data_dir is not None:
            self.load_data(data_dir)
    
    def analyze_bam(self, bamname, features, freqAT, seqref, cutoffAT=0.7, min_qual=30,
                    VERBOSE=False):
        print 'Analyzing', self.name
        bamfile = pysam.AlignmentFile(bamname, 'rb')
        self.chrom_len = bamfile.lengths[bamfile.gettid(self.name)]
        if self.name=='chrM':
            self.chrom_len = int(self.chrom_len)
        self.coverage = {key: np.zeros(self.chrom_len, dtype=int)
                         for key in covkeys}
        self.snp_freq = {key: np.zeros((self.chrom_len,3), dtype=float)
                         for key in covkeys}
        self.snps = {key: np.chararray((self.chrom_len,3))
                         for key in covkeys}
        self.deletions = {key: np.zeros(self.chrom_len, dtype=int)
                         for key in covkeys}
        self.insertions = {key: np.zeros(self.chrom_len, dtype=int)
                         for key in covkeys}
        self.consensus = {key: np.zeros(self.chrom_len, dtype='S1')
                         for key in covkeys}
        self.min_qual = min_qual
        self.freqAT  = freqAT
        self.cutoffAT = cutoffAT
        coding = coding_regions(features, self.chrom_len)
        for pi, pileup_column in enumerate(bamfile.pileup(self.name)):
            if VERBOSE and pi % 10000 == 0:
                print pi*1e-3
            ref_pos = pileup_column.reference_pos
            if self.name=='chrM':
                ref_pos = ref_pos % self.chrom_len
            if freqAT[ref_pos] >= cutoffAT:
                keyAT = 'high'
            else:
                keyAT = 'low'
            if coding[ref_pos]:
                keyCoding = 'code'
            else:
                keyCoding = 'non'
                
            alleles = Counter()
            for pileup_read in pileup_column.pileups:
                read = pileup_read.alignment
                if read.mapping_quality < min_qual:
                    continue
                if pileup_read.is_refskip:
                    continue
                if pileup_read.indel > 0:
                    self.insertions[keyCoding][ref_pos] += 1
                    self.insertions[keyAT][ref_pos] += 1
                    self.insertions['tot'][ref_pos] += 1
                    continue
                if pileup_read.is_del:
                    self.deletions[keyCoding][ref_pos] += 1
                    self.deletions[keyAT][ref_pos] += 1
                    self.deletions['tot'][ref_pos] += 1
                    continue
                read_pos = pileup_read.query_position
                if read.query_qualities[read_pos] < min_qual:
                    continue
                alleles[read.query_sequence[read_pos]] += 1
                
            cov = np.sum(alleles.values())  # alleles['A']+alleles['C']+alleles['G']+alleles['T']
            self.coverage[keyCoding][ref_pos] += cov
            self.coverage[keyAT][ref_pos] += cov
            self.coverage['tot'][ref_pos] += cov
            if cov > 0:
                counts = alleles.most_common()
                #if len(counts) > 1:
                #freq = counts[1][1] / cov
                self.consensus[keyCoding][ref_pos] = counts[0][0]
                self.consensus[keyAT][ref_pos] = counts[0][0]
                self.consensus['tot'][ref_pos] = counts[0][0]
                ordc = 0
                for (k,v) in counts:
                    if k!=seqref[ref_pos] and k in ['A','C','G','T']:
                        self.snp_freq[keyCoding][ref_pos][ordc] = v/cov
                        self.snp_freq[keyAT][ref_pos][ordc] = v/cov
                        self.snp_freq['tot'][ref_pos][ordc] = v/cov
                        self.snps[keyCoding][ref_pos][ordc] = k
                        self.snps[keyAT][ref_pos][ordc] = k
                        self.snps['tot'][ref_pos][ordc] = k
                        ordc += 1
        bamfile.close()
    
    def get_snps(self, min_freq, min_cov, max_freq=1, dtype='tot'):
        #coverage, snp_freq = self.get_data(dtype)
        #freqs = np.zeros((self.chrom_len,), dtype=float)
        #for i in range(3): 
        i = 0
        inds = np.logical_and(self.snp_freq[dtype][:,i] >= min_freq, self.coverage[dtype] >= min_cov)
        freqs = self.snp_freq[dtype][inds,i]
        locus = np.argwhere(inds)
        return freqs, locus.flatten()
    
    def snp_stats(self, min_freq, min_cov, dtype='tot'):
        #coverage, snp_freq = self.get_data(dtype)
        cov_inds = self.coverage[dtype] >= min_cov
        snp_inds = np.logical_and(cov_inds, self.snp_freq[dtype][:,0] >= min_freq)
        return np.count_nonzero(snp_inds), np.count_nonzero(cov_inds)
    
    def snp_cdf(self, min_freq, min_cov, dtype='tot', frac=True, norm=False):
        #coverage, snp_freq = self.get_data(dtype)
        cov_ind = np.count_nonzero(self.coverage[dtype] >= min_cov)
        freqs, locus = self.get_snps(min_freq, min_cov, dtype=dtype)
        #inds = np.nonzero(freqs)[0]
        counts = Counter(freqs)
        x, c =  cs.cdf(counts, norm=norm)
        if frac:
            c = c / cov_ind
        return x, c        
        
    def get_indels(self, min_freq, min_cov, dtype='tot'):
        cov = self.coverage[dtype]
        deletions =  self.deletions[dtype] #/ (cov+1e-8)
        insertions =  self.insertions[dtype] #/ (cov+1e-8)
        #dfreqs = np.zeros(self.chrom_len, dtype=float)
        #ifreqs = np.zeros(self.chrom_len, dtype=float)
        d_inds = np.logical_and(
            deletions >= min_freq*cov, cov >= min_cov)
        i_inds = np.logical_and(
            insertions >= min_freq*cov, cov >= min_cov)
        dfreqs = deletions[d_inds] / cov[d_inds]
        ifreqs = insertions[i_inds] / cov[i_inds]
        dlocus = np.argwhere(d_inds)
        ilocus = np.argwhere(i_inds)
        return (dfreqs, dlocus.flatten()), (ifreqs, ilocus.flatten())
        
    def plot_stats(self, min_freq, min_cov, features=None, ignore=None, dtype='tot', ltype = 'cov'):
        snp_freqs, snp_locus = self.get_snps(min_freq, min_cov, dtype)
        (del_freqs, del_locus) , (ins_freqs, ins_locus) = self.get_indels(min_freq, min_cov, dtype)
        #snp_freqs = ma.masked_where(snp_freqs == 0, snp_freqs)
        #del_freqs = ma.masked_where(del_freqs == 0, del_freqs)
        #ins_freqs = ma.masked_where(ins_freqs == 0, ins_freqs)
        #ref_cov   = np.ones(self.chrom_len, dtype=float)
        fig = plt.figure()
        if ltype == 'cov': 
            win = np.ones(1000)/1000
            cov = np.convolve(self.coverage[dtype], win, 'same')
            max_cov = np.max(cov)
            cov = cov / max_cov
            fig.gca().plot(np.array([0,self.chrom_len]),np.array([1,1])*min_cov/max_cov, 'C4--', label=r'$c_{\mathrm{min}}/c_{\mathrm{max}}$')
        else:
            cov = self.freqAT
        fig.gca().plot(cov, 'C9-', label=r'$c/c_{\mathrm{max}}$', alpha=0.5)
        fig.gca().plot(snp_locus, snp_freqs, 'C0.', label=r'$\nu_{\mathrm{SNP}}$')
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
            fig, gene_box, ori_box = add_features(fig, features, isW303=True)  #add_features(fig, features)
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
        
    def save_data(self, data_dir):
        save_name = data_dir + self.name + '.npz'
        np.savez_compressed(
            save_name, coverage=self.coverage, snp_freq=self.snp_freq, snps=self.snps,
            deletions=self.deletions, insertions=self.insertions, consensus=self.consensus,
            min_qual=self.min_qual, freqAT=self.freqAT, cutoffAT=self.cutoffAT)
    
    def load_data(self, data_dir):
        datafile = np.load(data_dir + self.name + '.npz')
        self.coverage = datafile['coverage'][()]
        self.snp_freq = datafile['snp_freq'][()]
        self.snps = datafile['snps'][()]
        self.insertions = datafile['insertions'][()]
        self.deletions = datafile['deletions'][()]
        self.consensus = datafile['consensus'][()]
        self.min_qual = datafile['min_qual'][()]
        self.freqAT  = datafile['freqAT'][()]
        self.cutoffAT = datafile['cutoffAT'][()]
        self.chrom_len = len(self.coverage['code'])


class GenomeSNPs:
    def __init__(self, samp_name, data_dir=None):
        self.name = samp_name
        self.mt = ChromosomeSNPs('chrM')
        self.nuc = ChromosomeSNPs('chrIV')
        if data_dir is not None:
            self.load_data(data_dir)
    
    def snp_stats(self, min_freq, min_cov):
        mt_snps, mt_cov = self.mt.snp_stats(min_freq, min_cov)
        nuc_snps, nuc_cov = self.nuc.snp_stats(min_freq, min_cov)
        print 'Chrom\t #SNPs\t #cov\t frac'
        print 'chrM\t {:d}\t {:d}\t {:0.2e}'.format(
            mt_snps, mt_cov, mt_snps/mt_cov)
        print 'chrIV\t {:d}\t {:d}\t {:0.2e}'.format(
            nuc_snps, nuc_cov, nuc_snps/mt_cov)
    
    def save_data(self, data_dir):
        save_dir = data_dir + self.name + '/'
        try:
            os.mkdir(save_dir)
        except OSError:
            pass
        self.mt.save_data(save_dir)
        self.nuc.save_data(save_dir)
    
    def load_data(self, data_dir):
        load_dir = data_dir + self.name + '/'
        self.mt.load_data(load_dir)
        self.nuc.load_data(load_dir)
