from __future__ import division
import os
import subprocess
import pdfkit
from fnmatch import filter
from natsort import natsorted
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import colors as clr
from matplotlib import patches
import matplotlib.pyplot as plt

min_freq = 0.15
min_cov  = 50
plt.ion()

save_dir= 'seqs/'
ref_dir = 'Chromosomes/'
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

class htmclr:
	# Text attributes
	START = '<span '
	STYLE = 'style="'
	ENDSTYLE = '">'
	END   = '</span>'
	BOLD = 'font-weight:bold;'
	UNDERLINE = 'text-decoration:underline;'

	# Foreground colors
	FG_BLACK = 'color:black;'
	FG_RED = 'color:red;'
	FG_GREEN = 'color:green;'
	FG_YELLOW = 'color:olive;'
	FG_BLUE = 'color:blue;'
	FG_MAGENTA = 'color:magenta;'
	FG_CYAN = 'color:cyan;'
	FG_WHITE = 'color:white;'

	# Background colors
	BG_BLACK = 'background-color:black;'
	BG_RED = 'background-color:red;'
	BG_GREEN = 'background-color:green;'
	BG_YELLOW = 'background-color:olive;'
	BG_BLUE = 'background-color:blue;'
	BG_MAGENTA = 'background-color:magenta;'
	BG_CYAN = 'background-color:cyan;'
	BG_WHITE = 'background-color:white;'

def getSeq(chrom):
	seq = str(SeqIO.read(ref_dir + chrom + '.fa', 'fasta').seq)
	return np.fromstring(seq, dtype='S1')

def Seq2Bin(seq):
	# convert a ATCG sequence to a binary seq
	slen = len(seq)
	bseq = np.zeros(2*slen)
	bseq[np.where(seq=='A')[0]] = 1
	bseq[np.where(seq=='A')[0]+slen] = 1
	bseq[np.where(seq=='T')[0]] = 1
	bseq[np.where(seq=='T')[0]+slen] = -1
	bseq[np.where(seq=='G')[0]] = -1
	bseq[np.where(seq=='G')[0]+slen] = 1
	bseq[np.where(seq=='C')[0]] = -1
	bseq[np.where(seq=='C')[0]+slen] = -1
	return bseq

def colorseq(seq, locus=[], freq=[]):
	cseq = ''
	for si, bp in enumerate(seq):
		if si in locus:
			if len(freq)>0:
				f = freq[np.argwhere(locus==si)]
			else:
				f = 0
			cseq += clr.FG_WHITE
			if bp == 'A':
				cseq += clr.BG_RED
			if bp == 'T':
				cseq += clr.BG_BLUE
			if bp == 'C':
				cseq += clr.BG_YELLOW
			if bp == 'G':
				cseq += clr.BG_GREEN
			if f>0.5:
				cseq += clr.BOLD
		else:
			#cseq += clr.BG_WHITE
			if bp == 'A':
				cseq += clr.FG_RED
			if bp == 'T':
				cseq += clr.FG_BLUE
			if bp == 'C':
				cseq += clr.FG_YELLOW
			if bp == 'G':
				cseq += clr.FG_GREEN
		cseq += bp+clr.ALL_OFF
	return cseq

def colorseqhtml(seq, locus=[], freq=[]):
	cseq = ''
	for si, bp in enumerate(seq):
		cseq += htmclr.START + htmclr.STYLE
		if si in locus:
			if len(freq)>0:
				f = freq[np.argwhere(locus==si)]
			else:
				f = 0
			cseq += htmclr.FG_WHITE
			if bp == 'A':
				cseq += htmclr.BG_RED
			if bp == 'T':
				cseq += htmclr.BG_BLUE
			if bp == 'C':
				cseq += htmclr.BG_YELLOW
			if bp == 'G':
				cseq += htmclr.BG_GREEN
			if f>0.5:
				cseq += htmclr.BOLD
		else:
			#cseq += htmclr.BG_BLACK
			if bp == 'A':
				cseq += htmclr.FG_RED
			if bp == 'T':
				cseq += htmclr.FG_BLUE
			if bp == 'C':
				cseq += htmclr.FG_YELLOW
			if bp == 'G':
				cseq += htmclr.FG_GREEN
		cseq += htmclr.ENDSTYLE+bp+htmclr.END
	return cseq

def printseq(samples, ref_seq, srange=[], width=100, htmlflag=False, dtype='tot'):
	locus = {}
	freqs = {}
	seq   = {}
	if len(srange) == 1:
		srange = [0,srange]
	if len(srange) == 0:
		srange = [0,len(ref_seq)]
	out_name = 'Mito_SNP_%d_%d' % (srange[0], srange[1])
	cseq = ''
	if htmlflag:
		cseq +='<html>\n<body>\n<pre>\n'

	for sname in natsorted(samples.keys()):
		sample = samples[sname]
		sname = sname.rsplit('_',4)[0]
		freqs1, locus1 = sample.mt.get_snps(min_freq, min_cov, dtype)
		locus1 = locus1.flatten()
		snps1 = sample.mt.snps[dtype][locus1,0]
		seq1  = np.copy(ref_seq)
		seq1[locus1] = snps1
		locus.update({sname: locus1})
		freqs.update({sname: freqs1})
		seq.update({sname: seq1})

	nline = int((srange[1]-srange[0])/width)
	nends = (srange[1]-srange[0]) - width*nline

	for i in range(nline):
		cseq += ('%-10.10s' % '') + ('%-10d' % (srange[0]+i*width)) + (repr((srange[0]+(i+1)*width)).rjust(width-10)) + '\n'
		for sname in natsorted(samples.keys()):
			sname = sname.rsplit('_',4)[0]
			locus1 = locus[sname] - srange[0] - i*width
			freqs1 = freqs[sname]
			seq1 = seq[sname][(srange[0]+i*width):(srange[0]+(i+1)*width)]
			if htmlflag:
				cseq += ('%-10.10s' % sname) + colorseqhtml(seq1, locus1, freqs1) + '\n'
			else:
				cseq += ('%-10.10s' % sname) + colorseq(seq1, locus1, freqs1) + '\n'
		cseq += '....................\n'
	if nends>0:
		if nends>20:
			cseq += ('%-10.10s' % '') + ('%-10d' % (srange[0]+nline*width)) + (repr(srange[1]).rjust(nends-10)) + '\n'
		else:
			cseq += ('%-10.10s' % '') + (repr(srange[1]).rjust(nends)) + '\n'
		for sname in natsorted(samples.keys()):
			sname = sname.rsplit('_',4)[0]
			locus1 = locus[sname] - srange[0] - nline*width
			freqs1 = freqs[sname]
			seq1 = seq[sname][(srange[0]+nline*width):srange[1]]
			if htmlflag:
				cseq += ('%-10.10s' % sname) + colorseqhtml(seq1, locus1, freqs1) + '\n'
			else:
				cseq += ('%-10.10s' % sname) + colorseq(seq1, locus1, freqs1) + '\n'

	if htmlflag:
		cseq +='</pre>\n</body>\n</html>\n'
		htmfile = open(save_dir+out_name+'.html', 'w')
		htmfile.write(cseq)
		htmfile.close()
		#p1 = subprocess.Popen(['echo', cseq], stdout=subprocess.PIPE)
		#p2 = subprocess.Popen(['aha > '+out_name+'.htm'], stdin=p1.stdout, stdout=subprocess.PIPE, shell=True)
		#pdfkit.from_file(out_name+'.htm', out_name+'.pdf')
		#print p2.communicate()
		#subprocess.call(cmd)
	else:
		with clr.pretty_output() as out:
			out.write(cseq)

def kmerFreq(seq, k=10, periodic=0):
	kmer_dic = {}
	slen = len(seq)
	if periodic:
		kmerlen = slen
	else:
		kmerlen = slen-k+1

	for i in np.arange(kmerlen):
		if i>slen-k and periodic:
			kmer = np.concatenate((seq[i:slen],seq[:i+k-slen])).tostring()
		else:
			kmer = seq[i:(i+k)].tostring()
		if kmer_dic.has_key(kmer):
			kmer_dic[kmer] += 1
		else:
			kmer_dic.update({kmer:1})
	return kmer_dic

def kmerCov(seq, kmer_dic, periodic=0):
	slen = len(seq)
	k = len(kmer_dic.keys()[0])

	kcov = np.zeros(slen)
	for i in np.arange(slen):
		"""
		if i<k-1:
			if periodic:
				kcov[i] = kmer_dic[]
			else:
				#kmers = [kmer_dic[seq[ii:ii+k].tostring()] for ii in np.arange(i+1)]
				kcov[i] = kmer_dic[seq[:k].tostring()]
		"""
		if i>slen-k:
			if periodic:
				kmer = np.concatenate((seq[i:],seq[:i+k-slen])).tostring()
			else:
				#kmers = [kmer_dic[seq[ii:ii+k].tostring()] for ii in np.arange(i-k+1,slen-k+1)]
				kmer = seq[-k:].tostring()
		else:
			kmer = seq[i:(i+k)].tostring()
			
		if kmer_dic.has_key(kmer):
			kcov[i] = kmer_dic[kmer]
	return kcov

def compKmer(kmer_ref, seqs, mincov=1, periodic=0):

	k = len(kmer_ref.keys()[0])

	kmer_dic = {}
	for key,vl in kmer_ref.items():
		if vl>=mincov:
			kmer_dic.update({key:0})

	for seq in seqs:
		seq_rc = str(Seq(seq.tostring()).reverse_complement())
		slen = len(seq)
		if periodic:
			kmerlen = slen
		else:
			kmerlen = slen-k+1

		for i in np.arange(kmerlen):
			if i>slen-k and periodic:
				kmer = np.concatenate((seq[i:slen],seq[:i+k-slen])).tostring()
				kmer_rc = seq_rc[2*slen-i-k:]+seq_rc[:slen-i]
			else:
				kmer = seq[i:(i+k)].tostring()
				kmer_rc = seq_rc[slen-(i+k):slen-i]
			if kmer_dic.has_key(kmer):
				kmer_dic[kmer] += 1
			elif kmer_dic.has_key(kmer_rc):
				kmer_dic[kmer_rc] += 1

	return kmer_dic


mt_chrom  = np.fromstring(str(SeqIO.read('Osman/mitodnalacostrain.fasta', 'fasta').seq.upper()),dtype='S1')  # getSeq('chrM')
mt_len = len(mt_chrom)
# nuc_chrom = getSeq('chrIV')

def plotCircle(fig,rpairs,totlen=mt_len,shif=0.2,alpha=0.2,color=['C6','C9'],ellips=1,clockwise=1,wt=[]):
	ax = fig.gca() #fig.add_subplot(111, aspect='auto')
	ac0 = patches.Arc((0,0),2,2,angle=0,color='black',linewidth=2)
	ax.add_patch(ac0)
	for i,r in enumerate(rpairs):
		r1 = abs(r[0])
		r2 = abs(r[1])
		if np.sign(r[0]*r[1])>0:
			cid = 0
		else:
			cid = 1
		t1 = r1/(totlen+1)*2*np.pi
		t2 = r2/(totlen+1)*2*np.pi
		if clockwise:
			t1 = np.pi/2-t1
			t2 = np.pi/2-t2
		x1 = np.cos(t1)
		y1 = np.sin(t1)
		x2 = np.cos(t2)
		y2 = np.sin(t2)
		x12 = (x1+x2)/2
		y12 = (y1+y2)/2
		r12 = np.sqrt(x12**2+y12**2)
		x0 = x12/r12*(1+shif)
		y0 = y12/r12*(1+shif)
		angle = 180*np.arctan(y0/x0)/np.pi
		x01 = x1-x0
		y01 = y1-y0
		xp  = (x01*x0+y01*y0)/(1+shif)
		yp2 = max(x01**2+y01**2-xp**2,1e-5)
		rp  = np.sqrt(xp**2+yp2)
		if ellips:
			a = 2*shif+(1-shif)*(rp-shif)/(np.sqrt(1+(1+shif)**2)-shif)
			b = np.sqrt(yp2/(1-xp**2/a**2))
		else:
			if rp<(1+shif):
				a = rp
				b = rp
			else:
				a = 1+shif
				b = np.sqrt(yp2/(1-xp**2/a**2))
		theta = 180*np.arctan(-np.sqrt(yp2)/abs(xp)*a/b)/np.pi
		if x0>0:
			theta1 = 180+theta
			theta2 = 180-theta
		else:
			theta1 = theta
			theta2 = -theta
		if len(wt)==0:
			ac1 = patches.Arc((x0,y0),2*a,2*b,angle=angle,theta1=theta1,theta2=theta2,alpha=alpha,color=color[cid],linewidth=2)
		else:
			ac1 = patches.Arc((x0,y0),2*a,2*b,angle=angle,theta1=theta1,theta2=theta2,alpha=wt[i]/max(wt),color=color[cid],linewidth=2)
		ax.add_patch(ac1)
	ax.axis('square')
	ax.set_xlim(-1.2,1.2)
	ax.set_ylim(-1.2,1.2)
	return fig

def featureCircle(fig,features,totlen=mt_len,radi=1.1,colorset={'rep_origin':'green','rRNA':'blue','tRNA':'magenta','CDS':'orange'},clockwise=1):
	ax = fig.gca()
	codons = []
	for f in features:
		dx = 0
		if totlen > 85779:
			if f.location.start.position>30000:
				dx = totlen-85779

		if f.type in colorset.keys():
			t1 = (f.location.start.position+dx)/(totlen+1)*360
			t2 = (f.location.end.position+dx)/(totlen+1)*360
			if clockwise:
				tp1 = t1
				t1 = 90-t2
				t2 = 90-tp1
			ac1 = patches.Arc((0,0),2*radi,2*radi,theta1=t1,theta2=t2,color=colorset[f.type],linewidth=2)
			ax.add_patch(ac1)
			if f.type=='CDS':
				for p in f.location.parts:
					pos = [p.start.position+dx,p.end.position+dx]
					if pos not in codons:
						codons.append(pos)
				#if f.qualifiers['gene'][0]=='COX1':
	for pos in codons:
		t1 = pos[0]/(totlen+1)*360
		t2 = pos[1]/(totlen+1)*360
		if clockwise:
			tp1 = t1
			t1 = 90-t2
			t2 = 90-tp1
		ac1 = patches.Arc((0,0),2*radi,2*radi,theta1=t1,theta2=t2,color='red',linewidth=2)
		ax.add_patch(ac1)
	return fig

def fdictCircle(fig,features,totlen=mt_len,radi=1.1,colorset={'rep_origin':'green','rRNA':'blue','rRNAintron':'blue','tRNA':'magenta','intron':'orange','gene':'orange','exon':'red'},clockwise=1):
	ax = fig.gca()
	codons = []
	for k,v in features.items():

		if k in colorset.keys():
			for i, t in enumerate(v):
				t1 = t[0]/(totlen+1)*360
				t2 = t[1]/(totlen+1)*360
				if clockwise:
					tp1 = t1
					t1 = 90-t2
					t2 = 90-tp1
				ac1 = patches.Arc((0,0),2*radi,2*radi,theta1=t1,theta2=t2,color=colorset[k],linewidth=2)
				ax.add_patch(ac1)
	return fig

def GCCircle(fig,GCdict,totlen=mt_len,radi=1.,colorset={'M1':'green','M1p':'green','M2':'red','M2p':'red','M2pp':'red','M3':'blue','M4':'magenta','G':'orange','V':'cyan'},clockwise=1):
	ax = fig.gca()
	for g in GCdict:
		dx = 0
		#if totlen>85779:
		#	if g['range'][0]>30000:
		#		dx = totlen-85779
		if g['GCcluster'] in colorset.keys():
			t1 = (g['range'][0]+dx)/(totlen+1)*2*np.pi
			t2 = (g['range'][1]+dx)/(totlen+1)*2*np.pi
			t0 = (t1+t2)/2
			if clockwise:
				#tp1 = t1
				#t1 = 90-t2
				#t2 = 90-tp1
				t0 = np.pi/2-t0
			ax.plot(radi*np.cos(t0),radi*np.sin(t0),'o',color=colorset[g['GCcluster']],markersize=3)
			#ax.add_patch(ac1)
	return fig

def plotcoverage(fig,cov,totlen=mt_len,shif=0.2,alpha=0.2,color='C1',range=[1.15,1.5],clockwise=1):
	ax = fig.gca() #fig.add_subplot(111, aspect='auto')
	if len(cov)!=totlen:
		totlen = len(cov)
	tt = np.arange(totlen)/totlen*2*np.pi
	if clockwise:
		tt = np.pi/2-tt
	mcov = max(cov)
	rr  = range[0]+(range[1]-range[0])*cov/mcov
	ax.plot(rr*np.cos(tt),rr*np.sin(tt),color=color)
	ax.axis('square')
	ax.set_xlim(-(range[1]+.1),(range[1]+.1))
	ax.set_ylim(-(range[1]+.1),(range[1]+.1))
	return fig

"""
ax = fig.gca()
for s in list(spotlocs):
	theta = s/(mt_len+1)*2*np.pi
	ax.plot([1*np.cos(theta),1.1*np.cos(theta)],[1*np.sin(theta),1.1*np.sin(theta)],'-C0',linewidth=1)

for feature in features:
	name = ''
	if feature.type == 'rep_origin':
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
	printseq(samples_L, mt_chrom, [start,end])
"""
