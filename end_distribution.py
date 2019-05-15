# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 11:07:30 2017

@author: kuns
"""

from __future__ import division
import numpy as np
import pysam
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.lines as mlines
import subprocess
import os

#plt.ion()

bamname = 'align11.sorted.bam'
cov_name = 'align11.npz'

data = np.load(cov_name)

bamfile = pysam.AlignmentFile(bamname, 'rb')

ref_start = {bamfile.references[ri]: np.zeros(bamfile.lengths[ri])
             for ri in range(len(bamfile.references))}
ref_end = {bamfile.references[ri]: np.zeros(bamfile.lengths[ri])
           for ri in range(len(bamfile.references))}   

ref_start_l = {bamfile.references[ri]: np.zeros(bamfile.lengths[ri])
               for ri in range(len(bamfile.references))}
ref_end_l = {bamfile.references[ri]: np.zeros(bamfile.lengths[ri])
             for ri in range(len(bamfile.references))} 
              
for ri, read in enumerate(bamfile):
    if ri % 10000 == 0:
        print ri
    if read.is_unmapped or read.is_secondary:
        continue
    chrom_name = bamfile.getrname(read.reference_id)
    ref_start[chrom_name][read.reference_start] += 1
    # "end" coordinate is 1 to the right of the last base of the mapping!!!
    ref_end[chrom_name][read.reference_end - 1] += 1
    if read.query_length > 1500:
        ref_start_l[chrom_name][read.reference_start] += 1
        ref_end_l[chrom_name][read.reference_end - 1] += 1

win = np.ones(1000)/1000

def plot_chrom(chrom_name):
    start_pos = np.convolve(ref_start[chrom_name], win, 'same')
    end_pos = np.convolve(ref_end[chrom_name], win, 'same')
    start_pos_l = np.convolve(ref_start_l[chrom_name], win, 'same')
    end_pos_l = np.convolve(ref_end_l[chrom_name], win, 'same')
    cov = data['ccov'][()][chrom_name]
    norm = np.max([np.max(start_pos), np.max(end_pos)]) / np.max(cov)
        
    fig = plt.figure()
#    fig.gca().plot(start_pos, 'C0-', label='start position')
#    fig.gca().plot(end_pos, 'C1-', label='end')
    fig.gca().plot(start_pos, 'C0-', label='total')
    fig.gca().plot(start_pos_l, 'C1-', label='long')
    fig.gca().plot(end_pos, 'C0--')
    fig.gca().plot(end_pos_l, 'C1--')
    fig.gca().plot(cov * norm, 'C2-', label='coverage', alpha=0.5)
    
    fig.gca().set_title(chrom_name + ' Mapping Positions')
    fig.gca().set_xlabel('position (kbp)')
    fig.gca().set_ylabel('normalized coverage')
    fig.gca().get_xaxis().set_major_formatter(
        mpl.ticker.FuncFormatter(lambda x, p: int(x*1e-3)))
    handles, labels = fig.gca().get_legend_handles_labels()
    line_start = mlines.Line2D([], [], color='k', linestyle='-')
    line_end = mlines.Line2D([], [], color='k', linestyle='--')
    handles.extend([line_start, line_end])
    labels.extend(['start position', 'end_position'])
    fig.gca().legend(handles, labels)
    return fig


chroms = ['chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrVI', 'chrVII',
          'chrVIII', 'chrIX', 'chrX', 'chrXI', 'chrXII', 'chrXIII',
          'chrXIV', 'chrXV', 'chrXVI', 'chrM']

fnames = []

for chrom in chroms:
    print chrom
    fig = plot_chrom(chrom)
    sname = chrom + '.pdf'
    fnames.append(sname)
    fig.savefig(sname)
    plt.close(fig)

cmd = ['pdftk']
cmd.extend(fnames)
cmd.extend(['output', 'end_distribution_long.pdf'])
subprocess.call(cmd)
for fname in fnames:
    os.remove(fname)
    