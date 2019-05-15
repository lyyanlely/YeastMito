# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 17:12:38 2017

@author: kuns
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


def analyze_sample(sample_name, data_dir, save_dir, min_qual=30):
    print 'Analyzing', sample_name
    bamname = data_dir + sample_name + '.sorted.bam'
    genome = GenomePolymorphisms(sample_name)
    genome.mt.analyze_bam(bamname, min_qual)
    genome.nuc.analyze_bam(bamname, min_qual)
    genome.save_data(save_dir)
    print 'Done saving', sample_name


def add_features(fig, features, ignore_features):
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
        fig.gca().axvspan(start, end, alpha=0.25, color=color)
        gene_box = mpatches.Patch(color='C1', alpha=0.25)
        ori_box = mpatches.Patch(color='C2', alpha=0.25)
    return fig, gene_box, ori_box


class ChromosomePolymorphisms:
    def __init__(self, chrom_name, data_dir=None):
        self.name = chrom_name
        self.chrom_len = None
        self.coverage = None
        self.deletions = None
        self.insertions = None
        self.snp_freq = None
        self.consensus = None
        self.snps = None
        self.min_qual = None
        if data_dir is not None:
            self.load_data(data_dir)
    
    def analyze_bam(self, bamname, min_qual=30, VERBOSE=False):
        print 'Analyzing', self.name
        bamfile = pysam.AlignmentFile(bamname, 'rb')
        self.chrom_len = bamfile.lengths[bamfile.gettid(self.name)]
        self.coverage = np.zeros(self.chrom_len, dtype=int)
        self.snp_freq = np.zeros(self.chrom_len, dtype=float)
        self.deletions = np.zeros(self.chrom_len, dtype=int)
        self.insertions = np.zeros(self.chrom_len, dtype=int)
        self.consensus = np.zeros(self.chrom_len, dtype='S1')
        self.snps = np.zeros(self.chrom_len, dtype='S1')
        self.min_qual = min_qual
        for pi, pileup_column in enumerate(bamfile.pileup(self.name)):
            if VERBOSE and pi % 10000 == 0:
                print pi*1e-3
            ref_pos = pileup_column.reference_pos
            alleles = Counter()
            for pileup_read in pileup_column.pileups:
                read = pileup_read.alignment
                if read.mapping_quality < min_qual:
                    continue
                self.coverage[ref_pos] += 1
                if pileup_read.indel > 0:
                    self.insertions[ref_pos] += 1
                if pileup_read.is_refskip:
                    continue
                if pileup_read.is_del:
                    self.deletions[ref_pos] += 1
                    continue
                read_pos = pileup_read.query_position
                if read.query_qualities[read_pos] < min_qual:
                    continue
                alleles[read.query_sequence[read_pos]] += 1
            cov = np.sum(alleles.values())
            #self.coverage[ref_pos] = cov
            if cov > 0:
                counts = alleles.most_common()
                self.consensus[ref_pos] = counts[0][0]
                if len(counts) > 1:
                    self.snps[ref_pos] = counts[1][0]
                    self.snp_freq[ref_pos] = counts[1][1] / cov
        bamfile.close()
    
    def get_snps(self, min_freq, min_cov):
        freqs = np.zeros(self.chrom_len, dtype=float)
        inds = np.logical_and(
            self.snp_freq >= min_freq, self.coverage >= min_cov)
        freqs[inds] = self.snp_freq[inds]
        return freqs
    
    def get_indels(self, min_freq, min_cov, norm=False):
        cov = self.coverage
        deletions =  self.deletions / (cov + 1e-6)
        insertions =  self.insertions / (cov + 1e-6)
        dfreqs = np.zeros(self.chrom_len, dtype=float)
        ifreqs = np.zeros(self.chrom_len, dtype=float)
        d_inds = np.logical_and(
            deletions >= min_freq, cov >= min_cov)
        i_inds = np.logical_and(
            insertions >= min_freq, cov >= min_cov)
        dfreqs[d_inds] = deletions[d_inds]
        ifreqs[i_inds] = insertions[i_inds]
        return dfreqs, ifreqs
    
    def snp_stats(self, min_freq, min_cov):
        cov_inds = self.coverage >= min_cov
        snp_inds = np.logical_and(cov_inds, self.snp_freq >= min_freq)
        return np.count_nonzero(snp_inds), np.count_nonzero(cov_inds)
    
    def snp_cdf(self, min_freq, min_cov, frac=True, norm=False):
        cov_ind = np.count_nonzero(self.coverage >= min_cov)
        freqs = self.get_snps(min_freq, min_cov)
        inds = np.nonzero(freqs)[0]
        counts = Counter(freqs[inds])
        x, c =  cs.cdf(counts, norm=norm)
        if frac:
            c = c / cov_ind
        return x, c

    def plot_stats(self, min_freq, min_cov, features=None, ignore=None):
        snp_freqs = self.get_snps(min_freq, min_cov)
        del_freqs , ins_freqs = self.get_indels(min_freq, min_cov, norm=True)
        snp_freqs = ma.masked_where(snp_freqs == 0, snp_freqs)
        del_freqs = ma.masked_where(del_freqs == 0, del_freqs)
        ins_freqs = ma.masked_where(ins_freqs == 0, ins_freqs)
        win = np.ones(1000)/1000
        cov = np.convolve(self.coverage, win, 'same')
        cov = cov / (np.max(cov) + 1e-6)
        fig = plt.figure()
        fig.gca().plot(cov, 'C9-', label=r'$c/c_{\mathrm{max}}$', alpha=0.5)
        fig.gca().plot(snp_freqs, 'C0.', label=r'$\nu_{\mathrm{SNP}}$')
        fig.gca().plot(del_freqs, 'C3.', label=r'$\nu_{\mathrm{del}}$')
        fig.gca().plot(ins_freqs, 'C8.', label=r'$\nu_{\mathrm{ins}}$')
        fig.gca().set_xlabel('position (kbp)')
        title = r'$\nu_{{\mathrm{{min}}}} = {:0.2f}$'.format(min_freq)
        title += r'$\qquad c_{{\mathrm{{max}}}} = {:d}$'.format(min_cov)
        fig.gca().set_title(title)
        fig.gca().get_xaxis().set_major_formatter(
            mpl.ticker.FuncFormatter(lambda x, p: int(x*1e-3)))
        if features is not None:
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
                if name in ignore:
                    continue
                start = feature.location.start.position
                end = feature.location.end.position
                fig.gca().axvspan(start, end, alpha=0.25, color=color)
            gene_box = mpatches.Patch(color='C1', alpha=0.25)
            ori_box = mpatches.Patch(color='C2', alpha=0.25)
            handles, labels = fig.gca().get_legend_handles_labels()
            phandles = handles[1:]
            plabels = labels[1:]
            phandles.extend([handles[0], gene_box, ori_box])
            plabels.extend([labels[0], 'gene', 'rep. origin'])
            fig.legend(phandles, plabels, ncol=2, loc='upper right')
        else:
            fig.legend()
        return fig
        
    
    def save_data(self, data_dir):
        save_name = data_dir + self.name + '.npz'
        np.savez_compressed(
            save_name, coverage=self.coverage, snp_freq=self.snp_freq,
            consensus=self.consensus, snps=self.snps, min_qual=self.min_qual,
            deletions=self.deletions, insertions=self.insertions)
    
    def load_data(self, data_dir):
        datafile = np.load(data_dir + self.name + '.npz')
        self.coverage = datafile['coverage'][()]
        self.snp_freq = datafile['snp_freq'][()]
        self.consensus = datafile['consensus'][()]
        self.snps = datafile['snps'][()]
        self.min_qual = datafile['min_qual'][()]
        self.chrom_len = len(self.coverage)
        self.deletions = datafile['deletions'][()]
        self.insertions = datafile['insertions'][()]


class GenomePolymorphisms:
    def __init__(self, samp_name, data_dir=None):
        self.name = samp_name
        self.mt = ChromosomePolymorphisms('chrM')
        self.nuc = ChromosomePolymorphisms('chrIV')
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