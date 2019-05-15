from __future__ import division
from Bio import SeqIO
import Bio
import numpy as np
import matplotlib.pyplot as plt
#from math import floor
from collections import Counter

winsize = 100
mt_name = 'chrM'

win_len = 30
win= np.ones(win_len) / win_len

bins = np.arange(0,1.01,0.02)
xx = np.arange(0.01,1.,0.02)

nuc_mcov = np.zeros(xx.shape)
mt_mcov  = np.zeros(xx.shape)

def getAT(seq):
	#seq = str(SeqIO.read(ref_dir + chrom + '.fa', 'fasta').seq)
	#seq = np.fromstring(str(seq), dtype='S1')
	at = np.logical_or(seq == 'A', seq == 'T')
	return np.convolve(at, win, 'same')

def atfrac(seq):
	if isinstance(seq,Bio.Seq.Seq):
		return (seq.count('A')+seq.count('T'))/len(seq)
	return ((seq=='A').sum()+(seq=='T').sum())/len(seq) #

nuc_cov = []
mt_cov  = []
nuc_val = []
mt_val  = []

genome = list(SeqIO.parse('BWAIndex/genome.fa','fasta'))
data = np.load('Desai/alignment_desai_test.npz')

inds = {}
#indends = {}
chromAT = {}

for chrom in genome:
	print chrom.id
	nwins  = int(len(chrom)/winsize)
	#covmat = np.reshape(data['cov'][()][chrom.id][0:(winsize*nwins)],(nwins,winsize))
	#seqmat = np.reshape(chrom.seq[0:(winsize*nwins)],(nwins,winsize))
	#val = [atfrac(seq) for seq in seqmat] [cov.mean() for cov in covmat]
	cov = data['cov'][()][chrom.id]
	val = getAT(chrom.seq)
	if chrom.id != mt_name:
		nuc_cov.extend(cov[win_len:-win_len])
		nuc_val.extend(val[win_len:-win_len])
	else:
		mt_gen = chrom
		mt_cov.extend(cov[win_len:-win_len])
		mt_val.extend(val[win_len:-win_len])
	indx = np.digitize(val,bins)
	indx[indx>len(xx)] = len(xx)
	inds.update({chrom.id: indx})
	#indends.update({chrom.id: np.digitize(atfrac(chrom.seq[(winsize*nwins):]),bins)})
	chromAT.update({chrom.id: val})

nuc_cov = np.array(nuc_cov)
mt_cov  = np.array(mt_cov)
nuc_inds = np.digitize(nuc_val,bins)
mt_inds = np.digitize(mt_val,bins)
for ii in range(len(xx)):
	nuc_mcov[ii] = np.mean(nuc_cov[nuc_inds==(ii+1)])
	mt_mcov[ii]  = np.mean(mt_cov[mt_inds==(ii+1)])

for ii in range(len(xx)):
	if ii<len(xx)-1 and nuc_mcov[ii]==0:
		nuc_mcov[ii] = (nuc_mcov[ii-1]+nuc_mcov[ii+1])/2
	else:
		if np.isnan(nuc_mcov[ii]):
			nuc_mcov[ii] = nuc_mcov[ii-1]

refcov = np.mean(nuc_cov[np.digitize(nuc_val,[0.4,0.6])==1])
covrat = refcov/nuc_mcov
covrat[np.isnan(covrat)] = 1 

#nwins  = int(len(data['cov'][()][mt_name])/winsize)
#mt_covmat = np.reshape(data['cov'][()][mt_name][0:(winsize*nwins)],(nwins,winsize))
#mt_endat  = atfrac((mt_gen.seq[(winsize*nwins):]))
#mt_covend = data['cov'][()][mt_name][(winsize*nwins):]*covrat[np.digitize(mt_endat,bins)-1]
#mt_covcor = np.append(np.array([cov*covrat[min(mt_inds[ri],len(xx))-1] for ri, cov in enumerate(covmat)]).flatten(),...
#					data['cov'][()][mt_name][(winsize*nwins):]*covrat[np.digitize(atfrac((mt_gen.seq[(winsize*nwins):])),bins)-1])
#mt_ccov = np.convolve(mt_covcor, win, 'same')

savename = 'Desai/Recalibration_Desai.npz'
np.savez_compressed(savename, 
	winsize=winsize, nbins=len(xx), bins=bins, xx=xx, covrat=covrat,
	inds=inds, chromAT=chromAT)

