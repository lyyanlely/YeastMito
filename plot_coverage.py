# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 23:19:52 2017

@author: kuns
"""

from __future__ import division
import numpy as np
import counter_stats as cs
from Bio import SeqIO
import subprocess
import os
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib as mpl

plt.ion()

mpl.rcParams.update({'text.usetex': True, 'font.family': 'serif',
                     'font.serif': 'Computer Modern Roman',
                     'font.size': 14,
                     'figure.autolayout': True})

fnames = []
savefigs = True

data11 = np.load('Desai/alignment_srr5406290.npz') #Yeast_Nanopore_Aug30_barcode07/alignment_barcode07.npz')
#data12 = np.load('Yeast_Nanopore_Aug30_barcode09/alignment_barcode09.npz')


x11_nuc, cdf11_nuc = cs.cdf(data11['nuc_len'][()])
x11_mt, cdf11_mt = cs.cdf(data11['mt_len'][()])
#x12_nuc, cdf12_nuc = cs.cdf(data12['nuc_len'][()])
#x12_mt, cdf12_mt = cs.cdf(data12['mt_len'][()])
#
fig_len = plt.figure()
fig_len.gca().plot(
    x11_mt, 1 - cdf11_mt, 'C0-', label='mitochondria')
fig_len.gca().plot(
    x11_nuc, 1 - cdf11_nuc, 'C1-', label='nuclear')
#fig_len.gca().semilogy(x12_mt * 1e-3, 1 - cdf12_mt, 'C0--')
#fig_len.gca().semilogy(x12_nuc * 1e-3, 1 - cdf12_nuc, 'C1--')
fig_len.gca().set_xlabel(r'$r$ (bp)')
##fig_len.gca().set_title('1 - CDF of read length')
fig_len.gca().set_title(r'$\mathrm{Prob}(\mathrm{read\; length} > r)$',
                        fontsize=14)
handles, labels = fig_len.gca().get_legend_handles_labels()
line11 = mlines.Line2D([], [], color='k', linestyle='-')
line12 = mlines.Line2D([], [], color='k', linestyle='--')
handles.extend([line11, line12])
labels.extend(['wildtype', 'petite mutant'])
fig_len.gca().legend(handles, labels)
if savefigs:
    fig_len.savefig('lengthcdf_desai_srr5406290.pdf')
    plt.close(fig_len)

short_lim = int(data11['short_lim'])
long_lim = int(data11['long_lim'])
norm_chrom = 'chrIV'
norm11 = np.mean(data11['ccov'][()][norm_chrom])
#norm12 = np.mean(data12['ccov'][()][norm_chrom])


def plot_chrom(chrom_name):
    cov11 = data11['ccov'][()][chrom_name] / norm11
    cov_short11 = data11['ccov_short'][()][chrom_name] / norm11
    cov_med11 = data11['ccov_med'][()][chrom_name] / norm11
    cov_long11 = data11['ccov_long'][()][chrom_name] / norm11
    #cov12 = data12['ccov'][()][chrom_name] / norm12
    #cov_short12 = data12['ccov_short'][()][chrom_name] / norm12
    #cov_med12 = data12['ccov_med'][()][chrom_name] / norm12
    #cov_long12 = data12['ccov_long'][()][chrom_name] / norm12
        
    fig = plt.figure()
    fig.gca().plot(cov11, 'C0-', label='total')
    lab_long = r'$\mathrm{{reads}}\geq {:d}\;\mathrm{{kbp}}$'.format(
        int(long_lim*1e-3))
    lab_medium = r'$\mathrm{{reads}}\in [{:d}, {:d}]\;\mathrm{{kbp}}$'.format(
        int(short_lim*1e-3), int(long_lim*1e-3))
    lab_short = r'$\mathrm{{reads}}\leq {:d}\;\mathrm{{kbp}}$'.format(
        int(short_lim*1e-3))
    fig.gca().plot(cov_long11, 'C1-', label=lab_long)
    fig.gca().plot(cov_med11, 'C2-', label=lab_medium)
    fig.gca().plot(cov_short11, 'C4-', label=lab_short)
    
    #fig.gca().plot(cov12, 'C0--')
    #fig.gca().plot(cov_long12, 'C1--')
    #fig.gca().plot(cov_med12, 'C2--')
    #fig.gca().plot(cov_short12, 'C4--')
    
    fig.gca().set_title(chrom_name + ' Coverage')
    fig.gca().set_xlabel('position (kbp)')
    fig.gca().set_ylabel('normalized coverage')
    fig.gca().get_xaxis().set_major_formatter(
        mpl.ticker.FuncFormatter(lambda x, p: int(x*1e-3)))
    handles, labels = fig.gca().get_legend_handles_labels()
    line11 = mlines.Line2D([], [], color='k', linestyle='-')
    line12 = mlines.Line2D([], [], color='k', linestyle='--')
    handles.extend([line11, line12])
    labels.extend(['wildtype', 'mutant'])
    fig.gca().legend(handles, labels, fontsize=8)
    return fig

#chroms = ['chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrVI', 'chrVII',
#          'chrVIII', 'chrIX', 'chrX', 'chrXI', 'chrXII', 'chrXIII',
#          'chrXIV', 'chrXV', 'chrXVI', 'chrM']
#
#for chrom in chroms:
#    print chrom
#    fig = plot_chrom(chrom)
#    if savefigs:
#        sname = chrom + '.pdf'
#        fig.savefig(sname)
#        fnames.append(sname)
#        plt.close(fig)
#
#
#if savefigs:
#    cmd = ['pdftk']
#    cmd.extend(fnames)
#    cmd.extend(['output', 'coverage_qual.pdf'])
#    subprocess.call(cmd)
#    for fname in fnames:
#        os.remove(fname)


record_file = '/data/kuns/kuns/GenBank/chrM.gb'
records = SeqIO.read(record_file, 'genbank')
features = records.features

ignore = ['AI5_ALPHA', 'AI5_BETA', 'AI4','AI3', 'AI2', 'AI1','BI4', 'BI3',
          'BI2', 'SCEI']

fig_mt = plot_chrom('chrM')
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
    fig_mt.gca().axvspan(start, end, alpha=0.25, color=color)

handles, labels = fig_mt.gca().get_legend_handles_labels()
line11 = mlines.Line2D([], [], color='k', linestyle='-')
line12 = mlines.Line2D([], [], color='k', linestyle='--')
handles.extend([line11, line12])
labels.extend(['wildtype', 'petite mutant'])
gene_box = mpatches.Patch(color='C1', alpha=0.25)
ori_box = mpatches.Patch(color='C2', alpha=0.25)
handles = [handles[0], handles[1], line11, gene_box, handles[2],
           handles[3], line12, ori_box]
labels = [labels[0], labels[1], 'wildtype', 'gene', labels[2], labels[3],
          'petite mutant', 'rep.\ origin']
fig_mt.gca().legend(handles, labels, ncol=2, fontsize=10, loc='upper left')
fig_mt.gca().set_title('Mitochondria Coverage')
fig_mt.gca().get_xaxis().set_major_formatter(
    mpl.ticker.FuncFormatter(lambda x, p: int(x*1e-3)))
fig_mt.gca().set_xlabel('position (kbp)')

if savefigs:
    fig_mt.savefig('mitoch_desai_srr5406290.pdf')
    plt.close(fig_mt)