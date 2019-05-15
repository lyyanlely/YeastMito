from __future__ import division
import os
import subprocess
import pdfkit
from fnmatch import filter
from natsort import natsorted
from collections import Counter
from snps import GenomeSNPs
from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq
import re
import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt
from Bio import Phylo
import pylab
from snps import origin_regions
from snps import coding_regions
from snps import exon_regions
from snps import RNA_regions


min_freq = 0.1
min_cov  = 30
min_qual = 20
min_len  = 14
def_dist = 10  # default distance to GC cluster
def_stem = 20  # default stem size
bon_alig = 1   # bonus of alignment
pen_mis  = 0  # penalty of mismatch
pen_gap  = -1.5  # penalty for making gap
pen_exgp = -.5  # penalty for extending gap
min_scor = 0.8
min_pscor= 0.7
min_q    = 0.1 # minimal overlap
w = 0.2

win_len = 30
win= np.ones(win_len) / win_len

#totlen =34

save_dir= 'seqs/'
ref_dir = 'Chromosomes/'
record_file = 'GenBank/chrM.gb'
records = SeqIO.read(record_file, 'genbank')
features = records.features
record_file = 'GenBank/chrIV.gb'
records = SeqIO.read(record_file, 'genbank')
features_nuc = records.features

ignore = ['AI5_ALPHA', 'AI5_BETA', 'AI4','AI3', 'AI2', 'AI1','BI4', 'BI3',
			'BI2', 'SCEI']

GCcluster = {'M1': Seq('TTCCGGGGCCCGGCCACGGGAGCCGGAACCCCGAAAGGAG'), 
			'M1p': Seq('TTCCCGCTTCGCGGGAACCCCGTAAGGAG'),
			'M2': Seq('TCCGGCCCGCCCCGCGGGGCGGACCCCGAAGGAG'),
			'M2p': Seq('TCCCCGCCCCGGCGGGGACCCCGAAGGAG'),
			'M2pp': Seq('TCCGGCCGAAGGAG'),
			'M3': Seq('TGAGGGACCCCTCCCTATACTATGGGAGGGGGACCGAACCCTTTAAAGAAGAG'),
			'M3p': Seq('TGAGGGACCCCTCCCTATGGGAGGGGGACCGAACCCCGAAGGAG'),
			'M4': Seq('TGAACACCTTTATTTAATTATAAAGGTGTGAACCAATCCGCAAGGCAAG'),
			'G': Seq('CCCGGTTTCTTACGAAACCGGGACCTCGGAGACGT'),
			'V': Seq('CCCCGCGGGCGCCAATCCGGTTGTTCACCGGATTGGTCCCGCGGGG')}

#TGAGGGACCCCCTCCCGTTA------GGGAGGGGGACCGAACCCCGAAGGAG
#TGAGGGACCCCCTCCC-TAACCTAATGGGAGGGGGACCGAACCCCGAAGGAG
#GGGGTTCGGTCCCCCTCCC-----TA-ACGGGAGGGGGTCCCTCA
#GGGGTTCGGTCCCCCTCCCGTTAGTACACGGGAGGGGGTCTCTCA
#CTAACCGCCCCCGCGGGGGCGGTTTATATAATTTAATAATTAATATATTAATAATAATTATAATAGTCCGGCCCGCCCCCTGCGGGGCGGACCCCGAAGGAG
#TGAGGGACCCCCTCCCTAACCTAATGGGAGGGGGACCGAACCCCGAAGGAG

GCstemloop = {'M1': {'Probability':.95, 'stem': [[0,10],[17,27]], 'loop': [11,16]},
			'M1p': {'Probability':.9, 'stem': [[0,5],[11,16]], 'loop': [6,10]},
			'M2': {'Probability':.99, 'stem': [[0,8],[14,22]], 'loop': [9,13]},
			'M2p': {'Probability':.8, 'stem': [[0,5],[12,17]], 'loop': [6,11]},
			'M2pp': {'Probability':.01, 'stem': None, 'loop': None},
			'M3': {'Probability':.5, 'stem': [[4,13],[24,34]], 'loop': [14,23]},
			'M3p': {'Probability':.5, 'stem': [[4,13],[24,34]], 'loop': [14,23]},
			'M4': {'Probability':.6, 'stem': [[3,11],[20,28]], 'loop': [12,19]},
			'G': {'Probability':.99, 'stem': [[0,7],[14,21]], 'loop': [8,13]},
			'V': {'Probability': .99, 'stem': [[0,19],[27,45]], 'loop': [20,26]}}

def get_label(leaf):
	return leaf.name

def getAT(seq):
	#seq = str(SeqIO.read(ref_dir + chrom + '.fa', 'fasta').seq)
	at = np.logical_or(seq == 'A', seq == 'T')
	return np.convolve(at, win, 'same')

def getSeq(chrom):
	seq = str(SeqIO.read(ref_dir + chrom + '.fa', 'fasta').seq)
	return np.frombuffer(seq, dtype='S1')

def gap0_function(x, y):  # x is gap position in seq, y is gap length
	#global totlen
	if y == 0:  # No gap
		return 0
	elif y == 1:  # Gap open penalty
		return pen_gap
	return pen_gap+(y-1)*pen_exgp

def gap_function(x, y):  # x is gap position in seq, y is gap length
	#global totlen
	if x>0 and x<totlen: 
		if y == 0:  # No gap
			return 0
		elif y == 1:  # Gap open penalty
			return pen_gap
		return pen_gap+(y-1)*pen_exgp
	return 0

def getOverlap(a, b):
	return max(0, min(a[1], b[1]) - max(a[0], b[0]))

def matchrange(align,lp=-1):
	if lp==-1:
		lp = len(align)
	startid = -1
	m_for = re.search(r'[^-]', align)
	if m_for:
		startid = m_for.start()
	endid = -1
	m_back = re.search(r'[^-]', align[::-1])
	if m_back:
		endid = lp-m_back.start()-1
	return [startid, endid]

def idPalin(seq):
	lp = len(seq)
	seq_f = seq
	seq_r = seq.reverse_complement()
	apalin= pairwise2.align.globalms(seq_f, seq_r, bon_alig, pen_mis, pen_gap, pen_exgp)[0]
	scor  = apalin[2]/apalin[4]
	atemp = apalin
	lptemp= lp
	ind_f = {i for i, x in enumerate(atemp[0]) if x == '-'}
	ind_r = {i for i, x in enumerate(atemp[1]) if x == '-'}
	ind = ind_f.union(ind_r)
	mrang = [0,lp]
	flag= 1
	while len(ind)>0 and flag==1:
		# no mismatch
		flag = 0
		mid1 = min(ind)
		if mid1 < atemp[4]/2: # on the left
			mid2 = mid1+1
			while mid2 in ind:
				mid2 += 1
			mrang[0] = mid2-sum(i < mid2 for i in ind_f)
			mrang[1] = lptemp-(mid2-sum(i < mid2 for i in ind_r))
			stemp_f = seq_f[mrang[0]:mrang[1]]
			stemp_r = stemp_f.reverse_complement()
			lptemp  = len(stemp_f)
			atemp   = pairwise2.align.globalms(stemp_f, stemp_r, bon_alig, pen_mis, pen_gap, pen_exgp)[0]
			if atemp[2]/atemp[4] > scor:
				flag  = 1
				scor  = atemp[2]/atemp[4]
				seq_f = stemp_f
				seq_r = stemp_r
				ind_f = {i for i, x in enumerate(atemp[0]) if x == '-'}
				ind_r = {i for i, x in enumerate(atemp[1]) if x == '-'}
				ind = ind_f.union(ind_r)

	mrang[0] = str(seq).index(str(seq_f))
	mrang[1] = mrang[0]+lptemp-1
	return (seq_f,mrang), scor


def idGC(seq, refAT=0.6, lmin=min_len, dist=def_dist):
	# identify the GC clusters and high GC regions in a long sequences
	seqAT = getAT(seq)
	seqlen = len(seq)
	mapped = np.zeros(seq.shape, dtype=int)  # mapped patch
	mapid = 0
	identifiedGC = []
	(GCpatches, npatch) = ndimage.label(seqAT<refAT)
	for i in range(npatch):
		(GCid,) = np.where(GCpatches==(i+1))
		#if 2558 in GCid:
		#	break
		lp = len(GCid)
		if lp >= lmin:
			seq_test = Seq(''.join(mt_chrom[max(0,(GCid[0]-dist)):min(seqlen,(GCid[-1]+dist+1))]))
			#print((GCid[0]-dist),(GCid[-1]+dist+1))
			#ifor = re.search(r'[^T]', str(seq_test)).start()
			#irev = re.search(r'[^T]', str(seq_test.reverse_complement())).start()
			#seq_test = seq_test[ifor:(len(seq_test)-irev)]
			GCidp = np.arange((GCid[0]-dist),(GCid[-1]+dist+1))
			#(seq_opt, prang), pscor = idPalin(seq_test)
			#if prang[1]-prang[0]+1 >= lmin and pscor > min_pscor:
			#else:
			#	palin = True
			#	pscor = 0
			#	prang = []

			lpp = lp+2*dist #-ifor-irev
			for refname, ref_clust in GCcluster.items():
				lref = len(ref_clust)
				global totlen
				totlen = lref
				align_f = pairwise2.align.globalmc(seq_test,ref_clust,bon_alig,pen_mis,gap_function,gap_function)[0]
				align_r = pairwise2.align.globalmc(seq_test,ref_clust.reverse_complement(),bon_alig,pen_mis,gap_function,gap_function)[0]
				if align_f[2]>align_r[2]:
					align = align_f
					direction = 'Forward'
				else:
					align = align_r
					direction = 'Reverse'
				alrang = matchrange(align[1],lpp)
				scor = align[2]/lref
				if scor > min_scor:
					if not any(mapped[GCidp[alrang[0]]:(GCidp[alrang[1]]+1)]):  # not mapped before
						mapid += 1
						mapped[GCidp[alrang[0]]:(GCidp[alrang[1]]+1)] = mapid
						identifiedGC.append({'GCcluster':refname,'range':GCidp[alrang],'score':scor,'direction':direction})
					else:
						mids = set(mapped[GCidp[alrang[0]]:(GCidp[alrang[1]]+1)]) - set([0])
						remap = []
						for mid in mids:
							if scor>identifiedGC[mid-1]['score']:
								rang = identifiedGC[mid-1]['range']
								if getOverlap(GCidp[alrang],rang)/(rang[1]-rang[0]+1)>min_q:
									remap.append(mid)
						if len(remap) > 0:  # needs to be remapped
							for mid in remap:
								rang = identifiedGC[mid-1]['range']
								mapped[rang[0]:(rang[1]+1)] = 0
							mapid += 1
							mapped[GCidp[alrang[0]]:(GCidp[alrang[1]]+1)] = mapid
							identifiedGC.append({'GCcluster':refname,'range':GCidp[alrang],'score':scor,'direction':direction})
						"""
							if len(remap) == 1: # only one set to remap
								mid = remap[0]
								rang = identifiedGC[mid]['range']
								mapped[rang[0]:(rang[1]+1)] = 0
								identifiedGC[mid] = {'GCcluster':refname,'range':GCidp[alrang],'score':scor}
								mapped[GCidp[alrang[0]]:(GCidp[alrang[1]]+1)] = mid
							else:
								mid = min(remap)
								remset = set(remap)-set([mid])
								rang = identifiedGC[mid]['range']
								mapped[rang[0]:(rang[1]+1)] = 0
								identifiedGC[mid] = {'GCcluster':refname,'range':GCidp[alrang],'score':scor}
								mapped[GCidp[alrang[0]]:(GCidp[alrang[1]]+1)] = mid
						"""

			if all(mapped[GCid]==0):
				mapid += 1
				mapped[GCid] = mapid
				identifiedGC.append({'GCcluster':'high GC','range':[GCid[0],GCid[-1]],'score': 0,'direction': None})



	return identifiedGC, mapped

def checkGC(loc, dist=def_dist):
	scor = 0
	for refname, ref_clust in GCcluster.items():
		lref  = len(ref_clust)  # length of reference
		seql  = Seq(''.join(mt_chrom[(loc+1):(loc+dist+lref+1)]))
		seqr  = Seq(''.join(mt_chrom[(loc-dist-lref):loc]))
		seqc  = Seq(''.join(mt_chrom[(loc-int((dist+lref)/2)):(loc+int((dist+lref)/2))]))
		lalign = pairwise2.align.globalms(seql,ref_clust,bon_alig,pen_mis,pen_gap,pen_exgp)[0]
		ralign = pairwise2.align.globalms(seqr,ref_clust,bon_alig,pen_mis,pen_gap,pen_exgp)[0]
		calign = pairwise2.align.globalms(seqc,ref_clust,bon_alig,pen_mis,pen_gap,pen_exgp)[0]
		lalign_r = pairwise2.align.globalms(seqr,ref_clust.reverse_complement(),bon_alig,pen_mis,pen_gap,pen_exgp)[0]
		ralign_r = pairwise2.align.globalms(seql,ref_clust.reverse_complement(),bon_alig,pen_mis,pen_gap,pen_exgp)[0]
		calign_r = pairwise2.align.globalms(seqc,ref_clust.reverse_complement(),bon_alig,pen_mis,pen_gap,pen_exgp)[0]
		ss = [lalign[2],ralign[2],calign[2],lalign_r[2],ralign_r[2],calign_r[2]]
		ms = max(ss)
		if ms/lref>scor:
			scor = ms/lref
			idx = ss.index(ms)
			if idx == 0:
				feat = [refname, 'Left', 'Forward']
			if idx == 1:
				feat = [refname, 'Right', 'Forward']
			if idx == 2:
				feat = [refname, 'Center', 'Forward']
			if idx == 3:
				feat = [refname, 'Left', 'Reverse']
			if idx == 4:
				feat = [refname, 'Right', 'Reverse']
			if idx == 5:
				feat = [refname, 'Center', 'Reverse']
	if scor > min_scor:
		return feat
	else:
		if getAT(mt_chrom[(loc-int(win_len/2)):(loc+int(win_len/2))])[int(win_len/2)]<0.5:
			return ['highGC']
		return None


def checkParlindrome(loc, stemsize=def_stem):
	seqr = Seq(''.join(mt_chrom[(loc+1):(loc+stemsize+1)]))
	seql = Seq(''.join(mt_chrom[(loc-stemsize):loc]))
	align = pairwise2.align.globalms(seql,seqr.reverse_complement(),bon_alig,pen_mis,pen_gap,pen_exgp)[0]
	if align[2]/stemsize > min_scor*0.5:
		return ['Loop']
	else:
		return None

"""
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

samples = {}
samples.update(samples_14E)
samples.update(samples_14)
samples.update(samples_08)
samples.update(samples_03)

samples_L = {sample: GenomeSNPs(sample, data_dir)
			for sample in names_L}
samples_M = {sample: GenomeSNPs(sample, data_dir)
			for sample in names_M}
samples_S = {sample: GenomeSNPs(sample, data_dir)
			for sample in names_S}
"""

mt_chrom  = np.frombuffer(str(SeqIO.read('DelATP21/W303_chrM_mito0.fasta', 'fasta').seq.upper()),dtype='S1') #getSeq('chrM') # Osman/mitodnalacostrain
mt_len    = len(mt_chrom)
nuc_chrom = getSeq('chrIV')
nuc_len   = len(nuc_chrom)



GCdict, GCmap = idGC(mt_chrom)

"""
hotspots = []
for sname in natsorted(samples.keys()):
	sample = samples[sname]
	freqs1, locus1 = sample.mt.get_snps(min_freq, min_cov)
	locus1 = locus1.flatten()
	snps1 = sample.mt.snps['tot'][locus1,0]
	hotspots.extend([(loc1,snps1[i]) for i,loc1 in enumerate(locus1)])
hotspots = list(set(hotspots))
spotlocs = set([loc[0] for loc in hotspots])

hs_cat = Counter()
locfeat = {}
for loc in natsorted(spotlocs):
	ft1 = checkParlindrome(loc)
	ft2 = checkGC(loc, dist=50)
	locfeat.update({loc: [ft1, ft2]})
	if ft1 is not None:
		if ft2 is not None:
			hs_cat['highGC & Parlindrome'] += 1
		else:
			hs_cat['Parlindrome'] += 1
	else:
		if ft2 is not None:
			hs_cat['highGC'] += 1
		else:
			hs_cat['Not Categorized'] += 1

GCloc = {}
for gcname in GCcluster.keys():
	loc = []
	dirc = []
	scor = []
	for i, v in enumerate(GCdict):
		if v['GCcluster']==gcname:
			loc.append(v['range'])
			dirc.append(v['direction'])
			scor.append(v['score'])
	GCloc.update({gcname: {'locations': loc,'directions':dirc, 'scores':scor}})

loc=[]
for i, v in enumerate(GCdict):
	if v['GCcluster']=='high GC':
		loc.append(v['range']) 
GCloc.update({'high GC':{'locations':loc}})

spotfeat = {}
for loc in natsorted(spotlocs):
	if GCmap[loc]==0:
		if GCmap[loc+15]==0:
			if GCmap[loc-15]==0:
				spotfeat.update({loc: None})
			else:
				clust = GCdict[GCmap[loc-15]-1]
				if clust['direction']=='Forward':
					rloc = loc-clust['range'][0]
				else:
					rloc = clust['range'][1]-loc
				spotfeat.update({loc:[None,[clust['GCcluster'],rloc,clust['direction']]]})
		else:
			clust = GCdict[GCmap[loc+15]-1]
			if clust['GCcluster'] != 'high GC':
				if clust['direction']=='Forward':
					rloc = loc-clust['range'][0]
				else:
					rloc = clust['range'][1]-loc
				spotfeat.update({loc:[None,[clust['GCcluster'],rloc,clust['direction']]]})
			else:
				spotfeat.update({loc:[None,[clust['GCcluster'], loc-clust['range'][0]]]})
	else:
		clust = GCdict[GCmap[loc]-1]
		if clust['GCcluster'] != 'high GC':
			if clust['direction']=='Forward':
				rloc = loc-clust['range'][0]
			else:
				rloc = clust['range'][1]-loc
			if GCstemloop[clust['GCcluster']]['Probability']>0.8:
				if rloc in range(GCstemloop[clust['GCcluster']]['loop'][0],GCstemloop[clust['GCcluster']]['loop'][1]+1):
					ft1 = 'loop'
				else:
					if rloc in range(GCstemloop[clust['GCcluster']]['stem'][0][0],GCstemloop[clust['GCcluster']]['stem'][0][1]+1) or rloc in range(GCstemloop[clust['GCcluster']]['stem'][1][0],GCstemloop[clust['GCcluster']]['stem'][1][1]+1):
						ft1 = 'stem'
					else:
						ft1 = None
			else:
				ft1 = None
			ft2 = [clust['GCcluster'], rloc, clust['direction']]
		else:
			ft1 = None
			ft2 = [clust['GCcluster'], loc-clust['range'][0]]
		spotfeat.update({loc:[ft1,ft2]})

for k,spot in enumerate(hotspots):
	loc = spot[0]
	snp = spot[1]
	if spotfeat[loc] is not None:
		gcname = spotfeat[loc][1][0]
		l0 = spotfeat[loc][1][1]
		mloc = []
		if gcname!='high GC':
			d0 = spotfeat[loc][1][2]
			locs = GCloc[gcname]['locations']
			dirc = GCloc[gcname]['directions']
			for i, l in enumerate(locs):
				if d0 == dirc[i]:   # same direction
					if d0=='Forward':
						l1 = l[0]+l0
					else:
						l1 = l[1]-l0
					if mt_chrom[l1]==snp:
						mloc.append(l1)
				else:
					if d0=='Forward':
						l1 = l[1]-l0
					else:
						l1 = l[0]+l0
					if mt_chrom[l1]==str(Seq(snp).complement()):
						mloc.append(l1)
		hotspots[k] = (loc,snp,mloc)

keys = GCcluster.keys()
keys.extend(['Uncat','high GC'])

gc_cat = Counter()
for sname in natsorted(samples_03.keys()):
	sample = samples_03[sname]
	freqs1, locus1 = sample.mt.get_snps(min_freq, min_cov)
	locus1 = locus1.flatten()
	snps1 = sample.mt.snps['tot'][locus1,0]
	for i, loc1 in enumerate(locus1):
		if spotfeat[loc1]:
			gc_cat[spotfeat[loc1][1][0]] += 1
		else:
			gc_cat['Uncat'] += 1

for key in keys:
	if key not in gc_cat.keys():
		gc_cat[key]=0

# common snp locations
spotcomm = {}  
for sname in natsorted(samples.keys()):
	sample = samples[sname]
	freqs1, locus1 = sample.mt.get_snps(min_freq, min_cov)
	locus1 = locus1.flatten()
	snps1 = sample.mt.snps['tot'][locus1,0]
	for i, loc1 in enumerate(locus1):
		if loc1 not in spotcomm.keys():
			spotcomm.update({loc1: [sname]})
		else:
			spotcomm[loc1].append(sname)

locs = spotcomm.keys()
bss  = mt_chrom[locs]

spotfreq = {loc1:{} for loc1 in locs}
for sname in natsorted(samples.keys()):
	sample = samples[sname]
	freqs1 = sample.mt.snp_freq['tot'][locs,:]
	snps1 = sample.mt.snps['tot'][locs,:]
	for i,loc1 in enumerate(locs):
		fa = 0
		fc = 0
		fg = 0
		ft = 0
		ftot = 0
		j0 = 0
		while j0<3 and freqs1[i,j0]>0:
			b1 = snps1[i,j0]
			if b1 == 'A':
				fa = freqs1[i,j0]
			elif b1 == 'C':
				fc = freqs1[i,j0]
			elif b1 == 'G':
				fg = freqs1[i,j0]
			else:
				ft = freqs1[i,j0]
			ftot += freqs1[i,j0]
			j0 += 1
		if bss[i] == 'A':
			fa = 1-ftot
		elif bss[i] == 'C':
			fc = 1-ftot
		elif bss[i] == 'G':
			fg = 1-ftot
		else:
			ft = 1-ftot
		spotfreq[loc1].update({sname:{'A':fa,'C':fc,'G':fg,'T':ft}})

def plotspot(loc):
	if loc in spotfreq.keys():
		spot = spotfreq[loc]
		fig = plt.figure()
		for i,k in enumerate(natsorted(spot.keys())):
			fig.gca().bar(i,spot[k]['A'],color='r',width=w,align='center')
			fig.gca().bar(i,spot[k]['C'],color='y',width=w,align='center',bottom=spot[k]['A'])
			fig.gca().bar(i,spot[k]['G'],color='g',width=w,align='center',bottom=spot[k]['A']+spot[k]['C'])
			fig.gca().bar(i,spot[k]['T'],color='b',width=w,align='center',bottom=spot[k]['A']+spot[k]['C']+spot[k]['G'])
		plt.xticks(range(len(spot)), natsorted(spot.keys()),rotation='vertical')
		plt.ylabel('frequency')
		return fig
	else:
		print 'Error! Not a SNP site'
		return None

freqc = {}
name1 = ['1914E11','1914E12','1914E14']
for loc in locs:
	for sname in name1:
		sample = samples[sname]
		if sample.mt.coverage['tot'][loc] > min_cov:
			f1 = sample.mt.snp_freq['tot'][loc,0]
			b1 = sample.mt.snps['tot'][loc,0]
			if f1 not in freqc.keys():
				freqc.update({f1:[]})
			name2 = filter(names_14E,sname+'*')
			name2.remove(sname)
			for s2 in name2:
				sample2 = samples[s2]
				if sample2.mt.coverage['tot'][loc] > min_cov and sample2.mt.snps['tot'][loc,0]==b1:
					f2 = sample2.mt.snp_freq['tot'][loc,0]
					freqc[f1].append(f2)

records = []
for sname in natsorted(samples.keys()):
	sample = samples[sname]
	freqs1, locus1 = sample.mt.get_snps(min_freq, min_cov,dtype='non')
	locus1 = locus1.flatten()
	snps1 = sample.mt.snps['non'][locus1,0]
	recs = SeqIO.parse('Mitochondria1.fa','fasta')
	r = recs.next()
	seq = np.fromstring(str(r.seq), dtype='S1')
	#for i,loc in enumerate(locus1):
	#	if spotfeat[loc]:
	#		seq[loc] = snps1[i]
	seq[locus1] = snps1
	r.seq = Seq(seq.tostring())
	r.id = sname
	r.name = sname
	records.append(r)

"""

origin = origin_regions(features, mt_len)
coding = coding_regions(features, mt_len)
exon = exon_regions(features, mt_len)
tRNA = RNA_regions(features, mt_len,'tRNA')
rRNA = RNA_regions(features, mt_len,'rRNA')
RNA = np.sign(tRNA+rRNA)
# exon
"""
ex_cat = Counter()
ncat = 0
for sname in natsorted(samples_08.keys()):
	sample = samples_08[sname]
	mcov = sample.mt.coverage['tot'].mean()
	if mcov > min_cov:
		ncat += 1
		freqs1, locus1 = sample.mt.get_snps(min_freq, min_cov)
		locus1 = locus1.flatten()
		snps1 = sample.mt.snps['tot'][locus1,0]
		for i, loc1 in enumerate(locus1):
			if coding[loc1]:
				if exon[loc1]:
					ex_cat['exon'] += 1
				else:
					ex_cat['intron'] += 1
			else:
				ex_cat['noncode'] += 1

fixed = [  796, 20933, 33742, 39516, 82506]
# codon cat
cd_cat = {'exon+GC':0,'exon+AT':0,'RNA+GC':0,'RNA+AT':0,'intron+GC':0,'intron+AT':0,'noncode+GC':0,'noncode+AT':0}
ncat = 0
for sname in natsorted(samples_08.keys()):
	sample = samples_08[sname]
	mcov = sample.mt.coverage['tot'].mean()
	if mcov > min_cov:
		ncat += 1
		freqs1, locus1 = sample.mt.get_snps(min_freq, min_cov)
		snps1 = sample.mt.snps['tot'][locus1,0]
		for i, loc1 in enumerate(locus1):
			if loc1 not in fixed:
				if coding[loc1]:
					if exon[loc1]:
						if spotfeat[loc1]:
							cd_cat['exon+GC'] += 1
						else:
							cd_cat['exon+AT'] += 1
					elif RNA[loc1]:
						if spotfeat[loc1]:
							cd_cat['RNA+GC'] += 1
						else:
							cd_cat['RNA+AT'] += 1
					else:
						if spotfeat[loc1]:
							cd_cat['intron+GC'] += 1
						else:
							cd_cat['intron+AT'] += 1
				else:
					if spotfeat[loc1]:
						cd_cat['noncode+GC'] += 1
					else:
						cd_cat['noncode+AT'] += 1
"""

def nearGC(loc,dist=15):
	if loc<dist:
		if GCmap[loc] or GCmap[loc+dist] or GCmap[loc-dist+len(mt_chrom)]:
			return True
		else:
			return False
	elif len(mt_chrom)-loc<=dist:
		if GCmap[loc] or GCmap[loc+dist-len(mt_chrom)] or GCmap[loc-dist]:
			return True
		else:
			return False
	else:
		if GCmap[loc] or GCmap[loc+dist] or GCmap[loc-dist]:
			return True
		else:
			return False

"""
nncat = np.zeros(8)
for i in range(len(mt_chrom)):
	if coding[i]:
		if exon[i]:
			if nearGC(i):
				nncat[5] += 1
			else:
				nncat[2] += 1
		elif RNA[i]:
			if nearGC(i):
				nncat[0] += 1
			else:
				nncat[1] += 1
		else:
			if nearGC(i):
				nncat[3] += 1
			else:
				nncat[7] += 1
	else:
		if nearGC(i):
			nncat[4] += 1
		else:
			nncat[6] += 1


f = 'path/to/my/file'
tree = Phylo.read(f, 'newick')
tree.ladderize()
Phylo.draw(tree, label_func=get_label, do_show=False)
pylab.axis('off')
pylab.savefig('tree2.svg',format='svg', bbox_inches='tight', dpi=300)

hs_numb = Counter()
for loc in natsorted(hotspots):
	if len(loc)>2:
		hs_numb[len(loc[2])] += 1

gcen_GC = {k:np.zeros(len(GCloc[k]['locations'])) for k in GCloc.keys()}
for k,v in GCloc.items():
    for i in range(len(v['locations'])):
        gcen_GC[k][i] = np.median(v['locations'][i])

ROmap = np.zeros(mt_len)
for fi in [0,1,2,3,5,6,7]:
	f  = origs[fi]
	fl = f.location.start.position
	fr = f.location.end.position
	ROmap[fl:fr] = fi

GEmap = np.zeros(mt_len)
for f in features:
	if f.type=='gene':
		fl = f.location.start.position
		fr = f.location.end.position
		GEmap[fl:fr] = 1

mt_tot = 286 # kuhn length=300
shif   = 128
rep_ori = [13,41,100,107,181,188,274]
M1 = [0,6,8,9,14,27,29,31,34,35,36,37,40,99,103,111,113,116,118,143,146,154,157,159,166,174,175,179,189,209,214,218,222,226,229,231,237,243,251,253,254,258,267,270,276,281]
# [-0,6,-8,9,-14,27,29,-31,34,35,36,37,-40,99,-103,-111,113,116,118,143,146,154,-157,159,166,-174,175,179,-189,209,214,218,222,-226,229,231,-237,243,-251,253,254,258,267,270,276,281]
M1p = [142, 234, 259, 266]  #  [142, 234, 259, -266]
M2  = [2,8,20,40,82,93,134,167,168,177,178,184,187,192,216,218,219,222,226,228,242,243,256,270,275,280,282,284]
# [-2,-8,-20,40,82,93,134,-167,-168,177,178,184,-187,192,-216,218,219,222,-226,-228,-242,-243,256,-270,275,-280,-282,284]
M2p = [9, 13, 113, 149, 188, 228, 243]  # [9, -13, 113, 149, -188, -228, 243]
M2pp= [20, 192, 216, 242, 282] # [-20, 192, -216, -242, -282]
M3  = [8, 113, 134, 285] # [-8, 113, 134, 285]
G = [14, 42, 100, 149, 189]  #  [-14, -42, 100, 149, -189]
V = [158, 163]  #  [158, 163]
gene = [[21, 26],[45, 88],[91, 92],[94, 97],[121, 145],[155, 156],[162, 166],[192, 207],[245, 252],[263, 266]]
# [[2, 2],[21, 26],[30, 31],[45, 88],[45, 76],[45, 72],[45, 66],[45, 62],[45, 53],[80, 83],[91, 92],[94, 97],[117, 117],[121, 145],[121, 140],[121, 133],[121, 128],[155, 156],[160, 160],[162, 166],[192, 207],[202, 205],[212, 212],[214, 214],[214, 215],[219, 220],[220, 220],[223, 223],[223, 224],[224, 224],[227, 227],[230, 230],[230, 230],[232, 232],[233, 233],[235, 235],[237, 237],[241, 241],[245, 247],[247, 252],[257, 257],[259, 260],[261, 261],[263, 266],[283, 283],[283, 285]]

entr = np.zeros(300)
for i in range(300):
	ids = []
	olc = []
	for fi in range(8):
	    f = origs[fi]
	    fl = f.location.start.position
	    fr = f.location.end.position
	    olc.append([fl,fr])
    	mtid = np.zeros(mt_len)
	    mtid[fl:fr] = 1
	    mtl = fr-fl
	    lid = int((fl-i)/300)
	    rid = int(ceil((fr-i)/300))
	    mid = int((fl/2+fr/2-i)/300)
	    ids.append([lid,rid,mid])
        for j in np.arange(lid,rid):
            rgid = np.zeros(mt_len)
            rgid[300*j+i:300*(j+1)+i] = 1
            p = np.sum(rgid*mtid)/mtl
            if p>0:
                entr[i] += -p*np.log(p)

entr = np.zeros(30)
for i in range(30):
	for k,v in GCloc.items():
		if k!='high GC':
			for fi in range(len(v['locations'])):
			    f = v['locations'][fi]
			    fl = f[0]
			    fr = f[1]
			    mtid = np.zeros(mt_len)
			    mtid[fl:fr] = 1
			    mtl = fr-fl
			    lid = int((fl-i)/30)
			    rid = int(ceil((fr-i)/30))
			    for j in np.arange(lid,rid):
			        rgid = np.zeros(mt_len)
			        rgid[30*j+i:30*(j+1)+i] = 1
			        p = np.sum(rgid*mtid)/mtl
			        if p>0:
			            entr[i] += -p*np.log(p)

GCpoly = {}
for k,v in GCloc.items():
	gcloc = v['locations']
	if k!='high GC':
		loc = []
		for fi in range(len(gcloc)):
			if v['directions'][fi]=='Reverse':
				sgn = -1
			else:
				sgn = 1
			loc.append(sgn*int((gcloc[fi][0]/2+gcloc[fi][1]/2-i)/300))
		GCpoly.update({k:loc})
        
dists_GC = {}
for k,v in gcen_GC.items():
	if len(v)>0:
		for ii in range(40):
			di = np.zeros(int(ni[ii]))
			for i in range(int(ni[ii])):
				di[i] = min(abs(v-np.median(np.where(pt[ii]==i+1))))
			dists_GC[k].update({ii:di})

gcdist_GC = {k:np.zeros(len(v),) for k,v in gcen_GC.items()}
for k,v in gcen_GC.items():
	gc1 = np.sort(v)
	for i in range(len(v)):
		if i==0:                               
			d1 = abs(gc1[i-1]-mt_len-gc1[i])
		else:
			d1 = abs(gc1[i-1]-gc1[i])
		if i==len(v)-1:
			d2 = abs(gc1[0]+mt_len-gc1[i])
		else:
			d2 = abs(gc1[i+1]-gc1[i])
		gcdist_GC[k][i] = min([d1,d2])

recbs = np.zeros((mt_len,))
d1 = 50
d2 = 75
for i in np.arange(len(loc1)):
	i1 = loc1[i]
	i2 = loc2[i]
	# [i1,i2] = data['pairs'][i,:]
	i1 = int(i1)
	i2 = int(i2)
	if i1>0:
		if i1+d2>mt_len:
			if i1+d1>mt_len:
				recbs[(i1+d1-mt_len):(i1+d2-mt_len)] += 1
			else:
				recbs[(i1+d1):] += 1
				recbs[:(i1+d2-mt_len)] += 1
		else:
			recbs[(i1+d1):(i1+d2)] += 1
	else:
		i1 = abs(i1)
		if i1<d2:
			if i1<d1:
				recbs[(i1-d2):(i1-d1)] += 1
			else:
				recbs[(i1-d2):] += 1
				recbs[:(i1-d1)] += 1
		else:
			recbs[(i1-d2):(i1-d1)] += 1
	if i2>0:
		if i2+d2>mt_len:
			if i2+d1>mt_len:
				recbs[(i2+d1-mt_len):(i2+d2-mt_len)] += 1
			else:
				recbs[(i2+d1):] += 1
				recbs[:(i2+d2-mt_len)] += 1
		else:
			recbs[(i2+d1):(i2+d2)] += 1
	else:
		i2 = abs(i2)
		if i2<d2:
			if i2<d1:
				recbs[(i2-d2):(i2-d1)] += 1
			else:
				recbs[(i2-d2):] += 1
				recbs[:(i2-d1)] += 1
		else:
			recbs[(i2-d2):(i2-d1)] += 1

paired = {}
for ri, read in enumerate(bf):
	if ri%1000000 == 0:
		print(ri)
	if read.is_secondary or read.is_unmapped:
		continue
	if bf.references[read.reference_id]=='chrM':
		if paired.has_key(read.query_name):
			paired[read.query_name].append(read)
		else:
			paired.update({read.query_name:[read]})

pileup = np.zeros(mt_len)  

for k,v in paired.items():
	if len(v) == 2:
		r1 = v[0]
		r2 = v[1]
		p1 = min(r1.get_reference_positions()[0],r2.get_reference_positions()[0])
		p2 = max(r1.get_reference_positions()[-1],r2.get_reference_positions()[-1])
		pileup[p1] += 1
		pileup[p2] += 1


npair = 0
nrevp = 0
pairs = []
plen  = []
mlen  = []
mqual = []
rqual = []
for k,v in paired.items():
	if len(v)>1:
		r1 = v[0]
		r2 = v[1]
		if r1.is_secondary or r2.is_secondary:
			continue
		#if r1.mapping_quality<min_qual or r2.mapping_quality<min_qual:
		#	continue
		npair += 1
		r1pos = r1.get_reference_positions()
		r2pos = r2.get_reference_positions()
		if np.logical_xor(r1.is_reverse, r2.is_reverse):
			pairs.append([np.mean(r1pos),np.mean(r2pos)])
		else:
			pairs.append([np.mean(r1pos),-np.mean(r2pos)])
			nrevp += 1
		mqual.append([r1.mapping_quality,r2.mapping_quality])
		rqual.append([np.mean(r1.query_qualities),np.mean(r2.query_qualities)])
		plen.append(min(abs(np.mean(r2pos)-np.mean(r1pos)),mt_len-abs(np.mean(r2pos)-np.mean(r1pos)))+len(r1pos)/2+len(r2pos)/2)
		mlen.append([len(r1pos),len(r2pos)])
pairs = np.array(pairs)
plen  = np.array(plen)
mlen  = np.array(mlen)
mqual = np.array(mqual)
rqual = np.array(rqual)

"""
#plt.bar(range(len(hs_cat)), hs_cat.values(), align='center')
#plt.xticks(range(len(hs_cat)), hs_cat.keys())
