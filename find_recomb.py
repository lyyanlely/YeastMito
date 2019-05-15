# -*- coding: utf-8 -*-
"""
Created on Sun Nov  5 12:40:34 2017

@author: lyan
"""

from __future__ import division
from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq
import re
import os
import numpy as np
from numpy.fft import fft, ifft
import pysam
#from joblib import Parallel, delayed
from scipy import ndimage

min_freq = 0.2
min_cov  = 20
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

win_len = 30
win= np.ones(win_len) / win_len

choplen = 1000
cpsize  = 70
stepsize= 10
min_step= 10
nmax = 2

totlen =34

min_qual = 20

ovlp_th = 0.5 # threshold of overlap
rc_th   = 2000 

n_jobs = 10

data_dir = '/data/lyan/yeast/nanopore/Dec17/reads/'  #Illumina/MitoSNPs/reads/' #
alig_dir = '/data/lyan/yeast/nanopore/Dec17/alignments/' #nanopore/Dec17/alignments/'
save_dir= 'recomb/' #nanopore/Dec17/recombinations/'
ref_dir = 'Chromosomes/'
record_file = 'GenBank/chrM.gb'
records = SeqIO.read(record_file, 'genbank')
features = records.features

def getAT(seq):
	#seq = str(SeqIO.read(ref_dir + chrom + '.fa', 'fasta').seq)
	at = np.logical_or(seq == 'A', seq == 'T')
	return np.convolve(at, win, 'same')

def ATfrac(seq):
	at = seq.count('A')+seq.count('T')
	return at/len(seq)

def getSeq(chrom):
	seq = str(SeqIO.read(ref_dir + chrom + '.fa', 'fasta').seq)
	return seq #np.fromstring(seq, dtype='S1')

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
 
def movedist(align):
    l0 = stepsize
    dl = align[0][:l0].count('-')
    while l0-dl<stepsize:
        l0 += dl
        dl = align[0][:l0].count('-')
    return l0-align[1][:l0].count('-')
 
def alignscore(seq1,ref1,seq2,ref2):
    a1 = pairwise2.align.globalms(seq1,ref1,bon_alig,pen_mis,pen_gap,pen_exgp)[0]
    a2 = pairwise2.align.globalms(seq2,ref2,bon_alig,pen_mis,pen_gap,pen_exgp)[0]
    #d1 = movedist(a1)
    #d2 = movedist(a2)
    return a1[2], a2[2] #, d1, d2
    
def periodicref(ref,start,end):
    if (start>=0 and end<=mt_len) or (start<0 and end<0):
        return ref[start:end]
    else:
        if start<0:
            return ''.join([ref[start:],ref[:end]])
        else:
            if start>=mt_len:
                return ref[(start-mt_len):(end-mt_len)]
            else:
                return ''.join([ref[start:],ref[:(end-mt_len)]])
    
def forwardloc(seq,start1,start2,dirc=True,rflg=False):
    if dirc:
        seq1 = seq[0:stepsize]
    else:
        seq1 = seq[-stepsize:]
    if rflg:
        seq2 = str(Seq(seq1).reverse_complement())
    else:
        seq2 = seq1
    if dirc:
        end1 = start1+stepsize
        if rflg:
            end2 = start2
            start2 = end2-stepsize
        else:
            end2 = start2+stepsize
    else:
        end1 = start1
        start1 = end1-stepsize
        if rflg:
            end2 = start2+stepsize
        else:
            end2 = start2
            start2 = end2-stepsize
    ref1 = periodicref(mt_chrom,start1,end1)
    ref2 = periodicref(mt_chrom,start2,end2)
    sc1, sc2 = alignscore(seq1,ref1,seq2,ref2)
    dsc1 = sc1
    dsc2 = sc2
    i = 1
    while dsc1<dsc2 and i*stepsize<len(seq):
        if i>=nmax:
            if dirc:
                seq1 = seq[(i-nmax+1)*stepsize:i*stepsize]
            else:
                seq1 = seq[-i*stepsize:-(i-nmax+1)*stepsize]
            if rflg:
                seq2 = str(Seq(seq1).reverse_complement())
            else:
                seq2 = seq1
            if dirc:
                start1 += stepsize # stepsize
                #end1 += stepsize
                if rflg:
                    #start2 -= stepsize
                    end2 -= stepsize #
                else:
                    start2 += stepsize #
                    #end2 += stepsize
            else:
                #start1 -= stepsize
                end1 -= stepsize #
                if rflg:
                    start2 += stepsize #
                    #end2 += stepsize
                else:
                    #start2 -= stepsize
                    end2 -= stepsize #
            ref1 = periodicref(mt_chrom,start1,end1)
            ref2 = periodicref(mt_chrom,start2,end2)
            #print len(ref1), len(ref2)
            sctmp1, sctmp2 = alignscore(seq1,ref1,seq2,ref2)
        else:
            #d1 = stepsize
            #d2 = stepsize
            sctmp1 = sc1
            sctmp2 = sc2
        i+=1
        if dirc:
            if i>nmax:
                seq1 = seq[(i-nmax)*stepsize:i*stepsize]
            else:
                seq1 = seq[:i*stepsize] #(i-1)*stepsize
        else:
            if i>nmax:
                seq1 = seq[-i*stepsize:-(i-nmax)*stepsize]
            else:
                seq1 = seq[-i*stepsize:] #(i-1)*stepsize
        if rflg:
            seq2 = str(Seq(seq1).reverse_complement())
        else:
            seq2 = seq1
        if dirc:
            #start1 += stepsize #
            end1 += stepsize
            if rflg:
                start2 -= stepsize
                #end2 -= stepsize #
            else:
                #start2 += stepsize #
                end2 += stepsize
        else:
            start1 -= stepsize
            #end1 -= stepsize #
            if rflg:
                #start2 += stepsize #
                end2 += stepsize
            else:
                start2 -= stepsize
                #end2 -= stepsize #
        ref1 = periodicref(mt_chrom,start1,end1)
        ref2 = periodicref(mt_chrom,start2,end2)
        #print len(ref1), len(ref2)
        sc1, sc2 = alignscore(seq1,ref1,seq2,ref2)
        dsc1 = sc1-sctmp1
        dsc2 = sc2-sctmp2
    #print i,dsc1,dsc2
        #if i==5:
        #    print sc1
    return i
    
def forwscore(seq,start,dirc=True,rflg=False):
    slen = int((len(seq)-cpsize)/stepsize)+1
    sc = np.zeros(slen)
    
    for i in range(slen):
        if dirc:
            seq1 = seq[i*stepsize:(i*stepsize+cpsize)]
        else:
            if i==0:
                seq1 = seq[-cpsize:]
            else:
                seq1 = seq[-(i*stepsize+cpsize):-i*stepsize]
        if rflg:
            seq1 = str(Seq(seq1).reverse_complement())
        if dirc:
            if rflg:
                end1 = start-i*stepsize
                start1 = end1-cpsize
            else:
                start1 = start+i*stepsize
                end1 = start1+cpsize
        else:
            if rflg:
                start1 = start+i*stepsize
                end1 = start1+cpsize
            else:
                end1 = start-i*stepsize
                start1 = end1-cpsize
        ref1 = periodicref(mt_chrom,start1,end1)
        a1 = pairwise2.align.globalms(seq1,ref1,bon_alig,pen_mis,pen_gap,pen_exgp)[0]
        sc[i] = a1[2]

    return sc
    
def findrecomb(read,end1,end2):
    if read.seq[0:choplen]==end2.seq or read.seq[0:choplen]==str(Seq(end2.seq).reverse_complement()):
        etmp = end1
        end1 = end2
        end2 = etmp
    if read.seq[0:choplen]==end1.seq:
        if read.reference_start != end1.reference_start:
            cg1 = end1.cigartuples[0]
            if cg1[0]==4:
                i1 = forwardloc(read.seq[cg1[1]:],read.reference_start,end1.reference_start)
            else:
                i1 = forwardloc(read.seq,read.reference_start,end1.reference_start)
        else:
            i1 = -1
    else:
        i1 = forwardloc(read.seq,read.reference_start,end1.reference_end,rflg=True)
    if read.seq[-choplen:]==end2.seq:
        if read.reference_end != end2.reference_end:
            i2 = forwardloc(read.seq,read.reference_end,end2.reference_end,dirc=False)
        else:
            i2 = -1
    else:
        i2 = forwardloc(read.seq,read.reference_end,end2.reference_start,dirc=False,rflg=True)
    return i1, i2
    
def readscores(read,end1,end2):
    if read.seq[0:choplen]==end2.seq or read.seq[0:choplen]==str(Seq(end2.seq).reverse_complement()):
        etmp = end1
        end1 = end2
        end2 = etmp
    if read.seq[0:choplen]==end1.seq:
        cg1 = end1.cigartuples[0]
        if cg1[0]==4:
            sc1 = forwscore(read.seq[cg1[1]:],end1.reference_start)
        else:
            sc1 = forwscore(read.seq,end1.reference_start)
    else:
        cg1 = end1.cigartuples[-1]
        if cg1[0]==4:
            sc1 = forwscore(read.seq[cg1[1]:],end1.reference_end,rflg=True)
        else:
            sc1 = forwscore(read.seq,end1.reference_end,rflg=True)
    if read.seq[-choplen:]==end2.seq:
        cg2 = end2.cigartuples[-1]
        if cg2[0]==4:
            sc2 = forwscore(read.seq[:-cg2[1]],end2.reference_end,dirc=False)
        else:
            sc2 = forwscore(read.seq,end2.reference_end,dirc=False)
    else:
        cg2 = end2.cigartuples[0]
        if cg2[0]==4:
            sc2 = forwscore(read.seq[:-cg2[1]],end2.reference_start,dirc=False,rflg=True)
        else:
            sc2 = forwscore(read.seq,end2.reference_start,dirc=False,rflg=True)
    return sc1,sc2
        
    
mt_chrom = getSeq('chrM')
mt_len = len(mt_chrom)

nuc_chrom = getSeq('chrIV')
nuc_len = len(nuc_chrom)

def analyzerecomb(it): # np.arange(12):
    print('analyzing barcode%02d' % it)
    bamfile = pysam.AlignmentFile(alig_dir+'barcode%02d_1.sorted.bam' % (it+1), 'rb')
    bf = pysam.AlignmentFile(alig_dir+'barcode%02d_1000.sorted.bam' % (it+1), 'rb')
    
    data = np.load(save_dir+'barcode%02d.npz' % (it+1))
    
    #mt_len = bf.lengths[bf.references.index('chrM')]
    #nuc_len = bf.lengths[bf.references.index('chrIV')]
    
    mtreads = []
    nucreads = []
    for read in bf:
        if bf.references[read.reference_id]=='chrM':
            if read.mapping_quality>min_qual and not read.is_secondary:
                mtreads.append(read)
        if bf.references[read.reference_id]=='chrIV':
            if read.mapping_quality>min_qual and not read.is_secondary:
                nucreads.append(read)
    
    names_uniq = data['names_uniq'][()]
    #if len(names_uniq):
    #    names_uniq = []
    nnames_uniq = data['nnames_uniq'][()]
    
    readid = data['readid'][()]
    nreadid = data['nreadid'][()]
            
        
    potreads = []
    potnreads = []
    for read in bamfile:
        if read.query_name in names_uniq:
            potreads.append(read)
        if read.query_name in nnames_uniq:
            potnreads.append(read)
    
    mtn=0
    potreads1 = []
    for read in potreads:
        if bamfile.references[read.reference_id]=='chrM' and not read.is_secondary:
            potreads1.append(read)
            mtn +=1
    nucn=0
    potnreads1 = []
    for read in potnreads:
        if bamfile.references[read.reference_id]=='chrIV' and not read.is_secondary:
            potnreads1.append(read)
            nucn +=1
            
    dist = data['dist'][()]
    distn = data['distn'][()]
    rid  = data['rid'][()]
    nrid = data['nrid'][()]
    
    scores = {}
    for ri,i in enumerate(rid):
        if ri*10 % len(rid)==0:
            print(ri/len(rid))
        read = potreads1[i]
        idx  = readid[read.query_name]
        end1 = mtreads[idx[0]]
        end2 = mtreads[idx[1]]
        sc1,sc2 = readscores(read,end1,end2)
        scores.update({i:np.concatenate((sc1,sc2),axis=0)})
        
    np.savez_compressed(save_dir+'barcode%02d' % (it+1),names_uniq=names_uniq,nnames_uniq=nnames_uniq,readid=readid,nreadid=nreadid,dist=dist,distn=distn,rid=rid,nrid=nrid,scores=scores)
    
def readaligns(fname,ref_name=''):
    bf = pysam.AlignmentFile(alig_dir+fname+'_%d-%d' % (cpsize,stepsize)+ref_name+'.sorted.bam','rb')
    #data = np.load(save_dir+'barcode%02d.npz' % it)
    #names_uniq = data['names_uniq'][()]
    reads = {} #{name: {} for name in names_uniq}
    for read in bf:
        rname = read.query_name
        if read.is_unmapped:
            continue
        r0 = rname.split('.')[0]
        r1 = rname.split('.')[1]
        if r0 not in reads.keys():
            reads.update({r0:{int(r1.split('_')[0]):read}})
        else:
            reads[r0].update({int(r1.split('_')[0]):read})
    return reads

def cutAT(reads1, reads2, chrom='chrM'):
    readAT = []
    mapAT  = []
    refAT  = []
    for k, r in reads1.items():
        ni = max(r.keys())
        flg1 = 0
        flg2 = 0
        i = 0
        if ni in r.keys():
            readAT.append(ATfrac(r[ni].seq[-stepsize:]))
        if i in r.keys():
            readAT.append(ATfrac(r[i].seq[0:stepsize]))
        elif i in reads2[k].keys():
            readAT.append(ATfrac(reads2[k][i].seq[0:stepsize]))
        while flg1*flg2==0 and i<int(ni/2):
            if flg1==0:
                if i in r.keys() and not r[i].is_unmapped and not r[i].is_secondary and r[i].reference_name==chrom: # mapped to mt_chrom
                    cg = r[i].cigartuples
                    m0 = 0
                    ii = 0
                    while ii<len(cg) and cg[ii][0]>0 :
                        m0 += cg[ii][1]
                        ii += 1
                    if r[i].is_reverse:
                        r0 = r[i].reference_end
                    else:
                        r0 = r[i].reference_start
                    mapAT.append(ATfrac(r[i].seq[max(0,m0-stepsize):min(cpsize,m0+stepsize)]))
                    if r0-stepsize<0:
                        seq = mt_chrom[r0-stepsize:]
                        seq += mt_chrom[:r0+stepsize]
                    elif r0+stepsize>mt_len:
                        seq = mt_chrom[r0-stepsize:]
                        seq += mt_chrom[:r0+stepsize-mt_len]
                    else:
                        seq = mt_chrom[r0-stepsize:r0+stepsize]
                    refAT.append(ATfrac(seq))
                    flg1 += 1
                else:
                    if i in reads2[k].keys() and not reads2[k][i].is_unmapped and not reads2[k][i].is_secondary and reads2[k][i].reference_name==chrom:
                        cg = reads2[k][i].cigartuples
                        m0 = 0
                        ii = 0
                        while ii<len(cg) and cg[ii][0]>0 :
                            m0 += cg[ii][1]
                            ii += 1
                        if reads2[k][i].is_reverse:
                            r0 = reads2[k][i].reference_end
                        else:
                            r0 = reads2[k][i].reference_start
                        r0 = (r0+int(mt_len/2))%mt_len
                        mapAT.append(ATfrac(reads2[k][i].seq[max(0,m0-stepsize):min(cpsize,m0+stepsize)]))
                        if r0-stepsize<0:
                            seq = mt_chrom[r0-stepsize:]
                            seq += mt_chrom[:r0+stepsize]
                        elif r0+stepsize>mt_len:
                            seq = mt_chrom[r0-stepsize:]
                            seq += mt_chrom[:r0+stepsize-mt_len]
                        else:
                            seq = mt_chrom[r0-stepsize:r0+stepsize]
                        refAT.append(ATfrac(seq))
                        flg1 += 1
            if flg2==0:
                if ni-i in r.keys() and not r[ni-i].is_unmapped and not r[ni-i].is_secondary and r[ni-i].reference_name==chrom: # mapped to mt_chrom
                    cg = r[ni-i].cigartuples
                    m0 = 0
                    ii = 0
                    while ii<len(cg) and cg[ii][0]>0 :
                        m0 += cg[ii][1]
                        ii += 1
                    if r[ni-i].is_reverse:
                        r0 = r[ni-i].reference_end
                    else:
                        r0 = r[ni-i].reference_start
                    mapAT.append(ATfrac(r[ni-i].seq[max(0,m0-stepsize):min(cpsize,m0+stepsize)]))
                    if r0-stepsize<0:
                        seq = mt_chrom[r0-stepsize:]
                        seq += mt_chrom[:r0+stepsize]
                    elif r0+stepsize>mt_len:
                        seq = mt_chrom[r0-stepsize:]
                        seq += mt_chrom[:r0+stepsize-mt_len]
                    else:
                        seq = mt_chrom[r0-stepsize:r0+stepsize]
                    refAT.append(ATfrac(seq))
                    flg2 += 1
                else:
                    if ni-i in reads2[k].keys() and not reads2[k][ni-i].is_unmapped and not reads2[k][ni-i].is_secondary and reads2[k][ni-i].reference_name==chrom:
                        cg = reads2[k][ni-i].cigartuples
                        m0 = 0
                        ii = 0
                        while ii<len(cg) and cg[ii][0]>0 :
                            m0 += cg[ii][1]
                            ii += 1
                        if reads2[k][ni-i].is_reverse:
                            r0 = reads2[k][ni-i].reference_end
                        else:
                            r0 = reads2[k][ni-i].reference_start
                        r0 = (r0+int(mt_len/2))%mt_len
                        mapAT.append(ATfrac(reads2[k][ni-i].seq[max(0,m0-stepsize):min(cpsize,m0+stepsize)]))
                        if r0-stepsize<0:
                            seq = mt_chrom[r0-stepsize:]
                            seq += mt_chrom[:r0+stepsize]
                        elif r0+stepsize>mt_len:
                            seq = mt_chrom[r0-stepsize:]
                            seq += mt_chrom[:r0+stepsize-mt_len]
                        else:
                            seq = mt_chrom[r0-stepsize:r0+stepsize]
                        refAT.append(ATfrac(seq))
                        flg2 += 1
            i += 1
    return readAT, mapAT, refAT

def cutATn(reads, chrom='chrIV'):
    readAT = []
    mapAT  = []
    refAT  = []
    for r in reads:
        readAT.append(ATfrac(r.seq[-stepsize:]))
        readAT.append(ATfrac(r.seq[0:stepsize]))
        rl = len(r.seq)
        cg = r.cigartuples
        ni = len(cg)
        m0 = 0
        ii = 0
        while ii<ni and cg[ii][0]>0:
            m0 += cg[ii][1]
            ii += 1
        m1 = 0
        ii = 0
        while ii<ni and cg[ni-ii-1][0]>0:
            m1 += cg[ni-ii-1][1]
            ii += 1
        if r.is_reverse:
            r0 = r.reference_end
            r1 = r.reference_start
        else:
            r0 = r.reference_start
            r1 = r.reference_end
        mapAT.append(ATfrac(r.seq[max(0,m0-stepsize):min(rl,m0+stepsize)]))
        mapAT.append(ATfrac(r.seq[max(0,m1-stepsize):min(rl,m1+stepsize)]))
        seq = nuc_chrom[max(0,r0-stepsize):min(nuc_len,r0+stepsize)]
        refAT.append(ATfrac(seq))
        seq = nuc_chrom[max(0,r1-stepsize):min(nuc_len,r1+stepsize)]
        refAT.append(ATfrac(seq))
    return readAT, mapAT, refAT

def overlaps(reads,chrom='chrM'):
    ovlp = {}
    for k,r in reads.items():
        if len(r)>1:
            #ni = int((150-cpsize)/stepsize)
            ni = max(r.keys())
            v = np.zeros(ni)
            for i in range(ni):  #max(r.keys())
                if i not in r.keys() or r[i].is_unmapped or r[i].mapping_quality<min_qual:
                    continue
                if i+1 not in r.keys() or r[i+1].is_unmapped or r[i+1].mapping_quality<min_qual:
                    continue
                if r[i].reference_name==chrom and r[i+1].reference_name==chrom:
                    if r[i].is_reverse != r[i+1].is_reverse:
                        sgn = -1
                    else:
                        sgn = 1
                    if r[i].reference_start>=r[i+1].reference_start and r[i].reference_start<r[i+1].reference_end:
                        v[i] = sgn*2/(1-stepsize/cpsize)*(r[i+1].reference_end-r[i].reference_start)/(r[i].reference_length+r[i+1].reference_length)
                    elif r[i+1].reference_start>=r[i].reference_start and r[i+1].reference_start<r[i].reference_end:
                        v[i] = sgn*2/(1-stepsize/cpsize)*(r[i].reference_end-r[i+1].reference_start)/(r[i].reference_length+r[i+1].reference_length)
            ovlp.update({k:v})
    return ovlp

def covadd(vec,i0,i1,clen=mt_len):
    if i0<0:
        if i1<0:
            vec[i0:i1] += 1
        else:
            vec[i0:] += 1
            vec[:i1] += 1
    elif i1>=clen:
        if i0>=clen:
            vec[i0-clen:i1-clen] += 1
        else:
            vec[i0:] += 1
            vec[:i1-clen] += 1
    else:
        vec[i0:i1] += 1
    return vec

def findhotspot(fname,chrom='chrM'):
    print('read '+fname) #barcode%02d' % it)
    reads = readaligns(fname,'_Mitochondria1')
    reads2 = readaligns(fname,'_Mitochondria2')
    print('compute score in '+fname) #barcode%02d' % it)
    ovlp = overlaps(reads,chrom)
    ovlp2 = overlaps(reads2,chrom)
    print('find recombination hotspots in '+fname) #barcode%02d' % it)
    recomb = np.zeros(mt_len)
    foldlp = np.zeros(mt_len)
    rcshort= np.zeros(mt_len)
    rclong = np.zeros(mt_len)
    rkeys = []
    lkeys = []
    rpair = []
    lpair = []
    for k,v in ovlp.items():
        if k in ovlp2.keys():
            if any(v<ovlp_th):
                locs = np.zeros(mt_len)
                locs[np.where((v<ovlp_th)*(ovlp2[k]<ovlp_th))[0]] = 1
                (pt,ni) = ndimage.label(locs)
                for lb in range(ni):
                    i = min(np.where(pt==lb+1)[0])
                    j = max(np.where(pt==lb+1)[0])+1
                    lc = (j-i)*stepsize #critirium for adjacent
                    #if (i not in reads[k].keys() and i not in reads2[k].keys()) or (j not in reads[k].keys() and j not in reads2[k].keys()):
                    #    continue
                    if i in reads[k].keys():
                        r01 = reads[k][i]
                    if i in reads2[k].keys():
                        r02 = reads2[k][i]
                    if j in reads[k].keys():
                        r11 = reads[k][j]
                    if j in reads2[k].keys():
                        r12 = reads2[k][j]
                    if i not in reads[k].keys() or (r01.is_unmapped or r01.reference_name!=chrom or r01.mapping_quality<min_qual):
                        if i not in reads2[k].keys() or (r02.is_unmapped or r02.reference_name!=chrom or r02.mapping_quality<min_qual):
                            continue
                        else:
                            r0 = r02
                            l00 = (r0.reference_start+int(mt_len/2))%mt_len
                            l01 = (r0.reference_end+int(mt_len/2))%mt_len
                    else:
                        r0 = r01
                        l00 = r0.reference_start
                        l01 = r0.reference_end
                    if j not in reads[k].keys() or (r11.is_unmapped or r11.reference_name!=chrom or r11.mapping_quality<min_qual):
                        if j not in reads2[k].keys() or (r12.is_unmapped or r12.reference_name!=chrom or r12.mapping_quality<min_qual):
                            continue
                        else:
                            r1 = r12
                            l10 = (r1.reference_start+int(mt_len/2))%mt_len
                            l11 = (r1.reference_end+int(mt_len/2))%mt_len
                    else:
                        r1 = r11
                        l10 = r1.reference_start
                        l11 = r1.reference_end
                    if r0.is_reverse:
                        sgn0 = -1
                    else:
                        sgn0 = 1
                    if r1.is_reverse:
                        sgn1 = -1
                    else:
                        sgn1 = 1
                    if sgn1!=sgn0:
                        dl = min([abs((l10-l00+l11-l01)/2),abs(mt_len+(l10-l00+l11-l01)/2)]) #+stepsize
                    elif sgn1*l11<sgn0*l01 or sgn1*l10<sgn0*l00:
                        if sgn1*l11<sgn0*l01 and sgn1*l10<sgn0*l00:
                            dl = mt_len+(sgn1*l10-sgn0*l00+sgn1*l11-sgn0*l01)/2
                        else:
                            dl = mt_len/2+(sgn1*l10-sgn0*l00+sgn1*l11-sgn0*l01)/2
                    else:
                        dl = (sgn1*l10-sgn0*l00+sgn1*l11-sgn0*l01)/2
                    if abs(dl/lc-1)>0.3:
                        if r0.is_reverse:
                            #if r0.reference_length<0.7*cpsize:
                            #i00 = l00
                            #i01 = l00+stepsize
                            #else:
                            i01 = l01-int(cpsize/2+lc/2-stepsize/2) #int((l00+l01)/2)
                            i00 = i01-stepsize
                        else:
                            #if r0.reference_length<0.7*cpsize:
                            #i01 = l01
                            #i00 = l01-stepsize
                            #else:
                            i00 = l00+int(cpsize/2+lc/2-stepsize/2) #int((l00+l01)/2)
                            i01 = i00+stepsize
                        if r1.is_reverse:
                            #if r1.reference_length<0.7*cpsize:
                            #i11 = l11
                            #i10 = l11-stepsize
                            #else:
                            i10 = l10+int(cpsize/2+lc/2-stepsize/2) #int((l10+l11)/2)
                            i11 = i10+stepsize
                        else:
                            #if r1.reference_length<0.7*cpsize:
                            #i10 = l10
                            #i11 = l10+stepsize
                            #else:
                            i11 = l11-int(cpsize/2+lc/2-stepsize/2) #int((l10+l11)/2)
                            i10 = i11-stepsize
                        i0 = (i00+i01)/2
                        i1 = (i10+i11)/2
                        if sgn1!=sgn0:
                            di = min([abs(i1-i0),abs(mt_len+i1-i0)]) #+stepsize
                        elif sgn1*i1<sgn0*i0:
                            di = mt_len+sgn1*i1-sgn0*i0
                        else:
                            di = sgn1*i1-sgn0*i0
                        if r0.is_reverse!=r1.is_reverse:  # a folded loop
                            foldlp = covadd(foldlp, i00, i01)
                            foldlp = covadd(foldlp, i10, i11)
                            lkeys.append(k)
                            lpair.append([sgn0*i0,sgn1*i1])
                        recomb = covadd(recomb, i00, i01)
                        recomb = covadd(recomb, i10, i11)
                        inser = lc - di
                        if abs(inser)>rc_th:
                            rclong = covadd(rclong, i00, i01)
                            rclong = covadd(rclong, i10, i11)
                        else:
                            rcshort = covadd(rcshort, i00, i01)
                            rcshort = covadd(rcshort, i10, i11)
                        rkeys.append(k)
                        rpair.append([sgn0*i0,sgn1*i1,inser])
    rkeys = list(set(rkeys))
    lkeys = list(set(lkeys))
    np.savez_compressed(save_dir+'recomb_'+fname+'_%d-%d' % (cpsize,stepsize), rc=recomb,lp=foldlp,rkeys=rkeys,lkeys=lkeys,rpair=rpair,lpair=lpair,rclong=rclong,rcshort=rcshort)
  

#names = [fname.split('_')[1] for fname in filter(os.listdir(save_dir),'*nucIV*')]

cpsize = 70
stepsize = 10
def seqpiece(seq,i0,i1,clen=mt_len):
    if i0<0:
        if i1<0:
            s1 = seq[i0:i1]
        else:
            s1 = seq[i0:]
            s1 += seq[:i1]
    elif i1>clen:
        if i0>clen:
            s1 = seq[(i0-clen):(i1-clen)]
        else:
            s1 = seq[i0:]
            s1 += seq[:(i1-clen)]
    else:
        s1 = seq[i0:i1]
    return s1

def compareRsite(rpair,rsize=cpsize,shif=[],ds = int(stepsize/2),chrom='chrM'):
    i0 = int(abs(rpair[0]))
    i1 = int(abs(rpair[1]))
    sgn0 = np.sign(rpair[0])
    sgn1 = np.sign(rpair[1])
    di = abs(i1-i0) #abs(rpair[2]) #
    if len(shif)==0:
        i00 = i0 - int(rsize/2)
        i01 = i0 + int(rsize/2)
        i10 = i1 - int(rsize/2)
        i11 = i1 + int(rsize/2)
        if sgn0==sgn1:
            sq0 = seqpiece(mt_chrom,i00,i01)
            sq1 = seqpiece(mt_chrom,i10,i11)
        else:
            sq0 = seqpiece(mt_chrom,i00,i01)
            sq1 = str(Seq(seqpiece(mt_chrom,i10,i11)).reverse_complement())
        a1 = pairwise2.align.globalxx(sq0,sq1)
        sc = a1[0][2]/a1[0][4]
        flg = 0
        dc = int(np.sign(np.random.rand(1)-0.5)[0])
        shif = 0
        while flg < 2 and abs(shif)<cpsize:
            i00 = i0 + shif + dc*ds - int(rsize/2)
            i01 = i0 + shif + dc*ds + int(rsize/2)
            i10 = i1 + shif + dc*ds - int(rsize/2)
            i11 = i1 + shif + dc*ds + int(rsize/2)
            if sgn0==sgn1:
                sq0 = seqpiece(mt_chrom,i00,i01)
                sq1 = seqpiece(mt_chrom,i10,i11)
            else:
                sq0 = seqpiece(mt_chrom,i00,i01)
                sq1 = str(Seq(seqpiece(mt_chrom,i10,i11)).reverse_complement())
            a1 = pairwise2.align.globalxx(sq0,sq1)
            sc0 = a1[0][2]/a1[0][4]
            if sc0>sc and abs(shif)<cpsize:
                shif += dc*ds
            else:
                flg += 1
                dc *= -1
    else:
        shif = shif[0]
        i00 = i0 + shif - int(rsize/2)
        i01 = i0 + shif + int(rsize/2)
        i10 = i1 + shif - int(rsize/2)
        i11 = i1 + shif + int(rsize/2)
        if sgn0==sgn1:
            sq0 = seqpiece(mt_chrom,i00,i01)
            sq1 = seqpiece(mt_chrom,i10,i11)
        else:
            sq0 = seqpiece(mt_chrom,i00,i01)
            sq1 = str(Seq(seqpiece(mt_chrom,i10,i11)).reverse_complement())
        a1 = pairwise2.align.globalxx(sq0,sq1)
        sc = a1[0][2]/a1[0][4]
    return sc,shif,di


def matchtest(rpair,rsize=cpsize,chrom='chrM'):
    i0 = int(abs(rpair[0]))
    i1 = int(abs(rpair[1]))
    sgn0 = np.sign(rpair[0])
    sgn1 = np.sign(rpair[1])
    i00 = i0 - int(rsize/2)
    i01 = i0 + int(rsize/2)
    i10 = i1 - int(rsize/2)
    i11 = i1 + int(rsize/2)
    if sgn0==sgn1:
        sq0 = seqpiece(mt_chrom,i00,i01)
        sq1 = seqpiece(mt_chrom,i10,i11)
    else:
        sq0 = seqpiece(mt_chrom,i00,i01)
        sq1 = str(Seq(seqpiece(mt_chrom,i10,i11)).reverse_complement())
    sc = np.sum(np.fromstring(sq0,dtype='S1')==np.fromstring(sq1,dtype='S1'))/rsize
    return sc

"""
out = Parallel(n_jobs=n_jobs)(
    delayed(findhotspot)(
    #sample) for sample in samples)
    'barcode%02d' %(it+1))  #
    for it in range(10))
    
print( 'Done analyze')

simsc = np.zeros((int(mt_len/30),mt_len))
for i in np.arange(int(mt_len/30)):
    for j in np.arange(mt_len):
        sc = matchtest([30*i+12,j],rsize=30)
        simsc[i,j] = sc

simscr = np.zeros((int(mt_len/30),mt_len))
for i in np.arange(int(mt_len/30)):
    for j in np.arange(mt_len):
        sc = matchtest([30*i+12,-j],rsize=30)
        simscr[i,j] = sc

simmx = np.zeros((int(mt_len/30),mt_len))
for i in np.arange(2696,int(mt_len/30)):
    for j in np.arange(mt_len):
        sc = []
        for di in range(11):
            sc.append(matchtest([30*i+12-5+di,j],rsize=30))
        simmx[i,j] = max(sc)

rcwt = np.zeros((int(mt_len/30),mt_len))
for i in np.arange(int(mt_len/30)):
    for j in np.arange(mt_len):
        sc = simsc[i,j]
        di = abs(j-30*i-12)
        if di==0:
            sc = 0
        if di<600 or di>mt_len-600:
            di = 600
        if ido[i] or idl[i]:
            pb = max(np.exp(12.1*sc-6.),np.exp(30.7*sc-19.))
        else:
            pb = max(np.exp(12.1*sc-6.),np.exp(30.7*sc-19.)/15)
        rcwt[i,j] = pb*(di*(1-di/mt_len))**(-2.5)

rcwtr = np.zeros((int(mt_len/30),mt_len))
for i in np.arange(int(mt_len/30)):
    for j in np.arange(mt_len):
        sc = simscr[i,j]
        di = abs(j-30*i-12)
        if di==0:
            sc = 0
        if di<600 or di>mt_len-600:
            di = 600
        if ido[i] or idl[i]:
            pb = max(np.exp(12.1*sc-6.),np.exp(30.7*sc-19.))
        else:
            pb = max(np.exp(12.1*sc-6.),np.exp(30.7*sc-19.)/15)
        rcwtr[i,j] = pb*(di*(1-di/mt_len))**(-2.5)

rcmx = np.zeros((int(mt_len/30),mt_len))
for i in np.arange(int(mt_len/30)):
    for j in np.arange(mt_len):
        sc = simmx[i,j]
        di = abs(j-30*i-12)
        if di<6:
            sc = 0
        if di<600 or di>mt_len-600:
            di = 600
        rcmx[i,j] = np.exp(17.13*sc-7.84)*(di*(1-di/mt_len))**(-2.5)

GCid = np.where(GCmap>0)[0]
ROid = np.where(ROmap>0)[0]
GEid = np.where(GEmap>0)[0]
sc1 = []
sf1 = []
di1 = []
id0 = {}
sc2 = []
sf2 = []
di2 = []
id2 = {}
sc3 = []
sf3 = []
di3 = []
id3 = {}
for i in range(len(rpair)):
    c,f,d = compareRsite(rpair[i],rsize=30)
    i0 = (int(abs(rpair[i][0]))+f+mt_len)%mt_len
    i1 = (int(abs(rpair[i][1]))+f+mt_len)%mt_len
    if i0 not in GCid and i1 not in GCid:
        sc2.append(c)
        sf2.append(f)
        di2.append(d)
        #if c>0.8:
        id2.update({i:[c,f,d]})
    if i0 not in GCid or i1 not in GCid:
        sc3.append(c)
        sf3.append(f)
        di3.append(d)
        #if c>0.8:
        id3.update({i:[c,f,d]})
    sc1.append(c)
    sf1.append(f)
    di1.append(d)
    if abs(d)<20:
        id0.update({i:[c,f,d]})

sc0 = []
sf0 = []
di0 = []
for i in range(30000):
    i0 = int(np.random.rand(1)*mt_len)
    i1 = int(np.random.rand(1)*mt_len)
    c,f,d = compareRsite([np.sign(np.random.rand(1)[0]-0.5)*i0,np.sign(np.random.rand(1)[0]-0.5)*i1,abs(i0-i1)],rsize=30)
    sc0.append(c)
    sf0.append(f)
    di0.append(d)

rc  = np.zeros(mt_len)
rcp = np.zeros(mt_len)
rcn = np.zeros(mt_len)
nrp = 0
idp = np.zeros(len(rpair))
nrn = 0
for i in range(len(rpair)):
    sgn0 = np.sign(rpair[i][0])
    sgn1 = np.sign(rpair[i][1])
    if sgn0*sgn1>0:
        idp[i] = 1
    if di1[i]>30:
        sf = sf1[i]
        i0 = int(abs(rpair[i][0]))
        i1 = int(abs(rpair[i][1]))
        l00 = i0+sf-5
        l01 = i0+sf+5
        l10 = i1+sf-5
        l11 = i1+sf+5
        rc = covadd(rc,l00,l01)
        rc = covadd(rc,l10,l11)
        if sgn0*sgn1<0:
            nrn += 1
            rcn = covadd(rcn,l00,l01)
            rcn = covadd(rcn,l10,l11)
        else:
            nrp += 1
            rcp = covadd(rcp,l00,l01)
            rcp = covadd(rcp,l10,l11)

# nanopore data ana
nidp = {}
ndi  = {}
idp = {}
di = {}
for it in range(10):
    idp_tmp = np.zeros(len(rpair[it+1]))
    di_tmp  = np.zeros(len(rpair[it+1]))
    for i in range(len(rpair[it+1])):
        sgn0 = np.sign(rpair[it+1][i][0])
        sgn1 = np.sign(rpair[it+1][i][1])
        if sgn0*sgn1>0:
            idp_tmp[i] = 1
        di_tmp[i] = min(abs(rpair[it+1][i][2]),mt_len-abs(rpair[it+1][i][2]))
    idp.update({it+1:idp_tmp})
    di.update({it+1:di_tmp})

    nidp_tmp = np.zeros(len(nrpair[it+1]))
    ndi_tmp  = np.zeros(len(nrpair[it+1]))
    for i in range(len(nrpair[it+1])):
        sgn0 = np.sign(nrpair[it+1][i][0])
        sgn1 = np.sign(nrpair[it+1][i][1])
        if sgn0*sgn1>0:
            nidp_tmp[i] = 1
        ndi_tmp[i] = min(abs(nrpair[it+1][i][2]),nuc_len-abs(nrpair[it+1][i][2]))
    nidp.update({it+1:nidp_tmp})
    ndi.update({it+1:ndi_tmp})

# merge different recombination events in nanopore
pairs = {}
dists = {}
counts = {}
covs = {}
c2vs = {}
cov2s = {}
rats  = {}
mcovs = {}
lcovs = {}

npairs = {}
ndists = {}
ncounts = {}
ncovs = {}
nc2vs = {}
ncov2s = {}
nrats  = {}
nmcovs = {}
nlcovs = {}
for it in np.arange(1,11): 
    data = np.load('recomb/recomb_barcode%02d_300-50.npz' % it)
    rpair = data['rpair'][()]

    data = np.load('recomb/recomb_barcode%02d_nucIV_300-50.npz' % it)
    nrpair = data['rpair'][()]

    data = np.load('Nanopore/Dec17/alignments/barcode%02d_1.npz' % it)
    coverage = data['cov'][()]['chrM']
    mcovs.update({it:coverage.mean()})
    lcovs.update({it:coverage})

    regpair = np.zeros((len(rpair),2))
    covpair = np.zeros((len(rpair),2))
    dpair = np.zeros((len(rpair),))
    for i in range(len(rpair)):
        sgn = np.sign(rpair[i,0])*np.sign(rpair[i,1])
        i0 = int(abs(rpair[i,0]))
        i1 = int(abs(rpair[i,1]))
        if i0<i1:
            regpair[i] = [i0,sgn*i1]
            covpair[i] = [max(1,coverage[i0]),max(1,coverage[i1%mt_len])]
        else:
            regpair[i] = [i1,sgn*i0]
            covpair[i] = [max(1,coverage[i1]),max(1,coverage[i0%mt_len])]
        dpair[i] = min(rpair[i,2],mt_len-rpair[i,2])

    pset = []
    dset = []
    cvset = []
    cc   = []
    for i in range(len(regpair)):
        if not pset:  # empty
            pset.append(regpair[i])
            cvset.append(covpair[i])
            dset.append(dpair[i])
            cc.append(1)
        else: # not empty
            dp  = pset - regpair[i]
            dp1 = dp[:,0]
            dp2 = dp[:,1]
            dpp = abs(dp1-dp2)
            #id  = np.argmin(dpp)
            id  = np.argmin(abs(dp1))
            if dpp[id]<150 and abs(dp1[id])<150:
                pset[id] = (cc[id]*pset[id]+regpair[i])/(cc[id]+1)
                cvset[id] = (cc[id]*cvset[id]+covpair[i])/(cc[id]+1)
                dset[id] = (cc[id]*dset[id]+dpair[i])/(cc[id]+1)
                cc[id]  += 1
            else:
                pset.append(regpair[i])
                cvset.append(covpair[i])
                dset.append(dpair[i])
                cc.append(1)
    pset = np.array(pset)
    cvset = np.array(cvset)
    dset = np.array(dset)
    cv2 = cvset[:,0]*cvset[:,1]
    c2v = cvset[:,0]+cvset[:,1]
    cc = np.array(cc)

    pairs.update({it:pset})
    dists.update({it:dset})
    covs.update({it:cvset})
    cov2s.update({it:cv2})
    c2vs.update({it:c2v})
    counts.update({it:cc})
    rats.update({it: counts[it]/cov2s[it]})

#   for nuclear chromosome
    coverage = data['cov'][()]['chrIV']
    nmcovs.update({it:coverage.mean()})
    nlcovs.update({it:coverage})

    regpair = np.zeros((len(nrpair),2))
    covpair = np.zeros((len(nrpair),2))
    dpair = np.zeros((len(nrpair),))
    for i in range(len(nrpair)):
        sgn = np.sign(nrpair[i,0])*np.sign(nrpair[i,1])
        i0 = int(abs(nrpair[i,0]))
        i1 = int(abs(nrpair[i,1]))
        if i0<i1:
            regpair[i] = [i0,sgn*i1]
            covpair[i] = [max(1,coverage[i0]),max(1,coverage[i1])]
        else:
            regpair[i] = [i1,sgn*i0]
            covpair[i] = [max(1,coverage[i1]),max(1,coverage[i0])]
        dpair[i] = min(nrpair[i,2],nuc_len-nrpair[i,2])

    pset = []
    dset = []
    cvset = []
    cc   = []
    for i in range(len(regpair)):
        if not pset:  # empty
            pset.append(regpair[i])
            cvset.append(covpair[i])
            dset.append(dpair[i])
            cc.append(1)
        else: # not empty
            dp  = pset - regpair[i]
            dp1 = dp[:,0]
            dp2 = dp[:,1]
            dpp = abs(dp1-dp2)
            #id  = np.argmin(dpp)
            id  = np.argmin(abs(dp1))
            if dpp[id]<150 and abs(dp1[id])<150:
                pset[id] = (cc[id]*pset[id]+regpair[i])/(cc[id]+1)
                cvset[id] = (cc[id]*cvset[id]+covpair[i])/(cc[id]+1)
                dset[id] = (cc[id]*dset[id]+dpair[i])/(cc[id]+1)
                cc[id]  += 1
            else:
                pset.append(regpair[i])
                cvset.append(covpair[i])
                dset.append(dpair[i])
                cc.append(1)
    pset = np.array(pset)
    cvset = np.array(cvset)
    dset = np.array(dset)
    cv2 = cvset[:,0]*cvset[:,1]
    c2v = cvset[:,0]+cvset[:,1]
    cc = np.array(cc)

    npairs.update({it:pset})
    ndists.update({it:dset})
    ncovs.update({it:cvset})
    ncov2s.update({it:cv2})
    nc2vs.update({it:c2v})
    ncounts.update({it:cc})
    nrats.update({it: counts[it]/cov2s[it]})

# compute the similarity scores for recombination pairs
sc1 = []
id0 = []
sc2 = []
id2 = []
sc3 = []
id3 = []
GCid = np.where(GCmap>0)[0]
for i in range(len(rpair)):
    c = matchtest(rpair[i,:],rsize=30)
    i0 = (int(abs(rpair[i,0]))+mt_len)%mt_len
    i1 = (int(abs(rpair[i,1]))+mt_len)%mt_len
    if i0 not in GCid and i1 not in GCid:
        sc2.append(c)
        id2.append(i)
    if i0 not in GCid or i1 not in GCid:
        sc3.append(c)
        id3.append(i)
    sc1.append(c)
    if abs(rpair[i,2])<20:
        id0.append(i)

ida = np.zeros(len(rpair))
for i in range(len(rpair)):
    i0 = (int(abs(rpair[i,0]))+mt_len)%mt_len
    i1 = (int(abs(rpair[i,1]))+mt_len)%mt_len
    if i0 in hgid or i1 in hgid:
        ida[i] = 1

(y,x) = np.histogram(simsc[np.logical_not(idc),:].reshape(((2859-sum(idc))*85779,)),weights=wt.reshape(((2859-sum(idc))*85779,)),bins=np.arange(0.5/30,1.02,1/30))

# merge different recombination events in illumina
pairs = {}
counts = {}
covs = {}
c2vs = {}
cov2s = {}
rats  = {}
mcovs = {}
lcovs = {}
nforw = {}
for name in names: 
    data = np.load('recomb/recomb_'+name+'_70-10.npz')
    rpair = data['rpair'][()]

    data = np.load('snps/'+name+'/chrM.npz')
    coverage = data['coverage'][()]['tot']
    mcovs.update({name:coverage.mean()})
    lcovs.update({name:coverage})

    regpair = np.zeros((len(rpair),2))
    covpair = np.zeros((len(rpair),2))
    for i in range(len(rpair)):
        sgn = np.sign(rpair[i,0])*np.sign(rpair[i,1])
        i0 = int(abs(rpair[i,0]))
        i1 = int(abs(rpair[i,1]))
        if i0<i1:
            regpair[i] = [i0,sgn*i1]
            covpair[i] = [max(1,coverage[i0]),max(1,coverage[i1%mt_len])]
        else:
            regpair[i] = [i1,sgn*i0]
            covpair[i] = [max(1,coverage[i1]),max(1,coverage[i0%mt_len])]

    pset = []
    cvset = []
    cc   = []
    for i in range(len(regpair)):
        if not pset:  # empty
            pset.append(regpair[i])
            cvset.append(covpair[i])
            cc.append(1)
        else: # not empty
            dp  = pset - regpair[i]
            dp1 = dp[:,0]
            dp2 = dp[:,1]
            dpp = abs(dp1-dp2)
            #id  = np.argmin(dpp)
            id  = np.argmin(abs(dp1))
            if dpp[id]<5 and abs(dp1[id])<15:
                pset[id] = (cc[id]*pset[id]+regpair[i])/(cc[id]+1)
                cvset[id] = (cc[id]*cvset[id]+covpair[i])/(cc[id]+1)
                cc[id]  += 1
            else:
                pset.append(regpair[i])
                cvset.append(covpair[i])
                cc.append(1)
    pset = np.array(pset)
    cvset = np.array(cvset)
    cv2 = cvset[:,0]*cvset[:,1]
    c2v = cvset[:,0]+cvset[:,1]
    cc = np.array(cc)

    pairs.update({name:pset})
    covs.update({name:cvset})
    cov2s.update({name:cv2})
    c2vs.update({name:c2v})
    counts.update({name:cc})
    rats.update({name: counts[name]/cov2s[name]})

# for nuclear chromosome
npairs = {}
ncounts = {}
ncovs = {}
nc2vs = {}
ncov2s = {}
nrats  = {}
nmcovs = {}
nlcovs = {}
for name in names: 
    data = np.load('recomb/recomb_'+name+'_nucIV_70-10.npz')
    rpair = data['rpair'][()]

    data = np.load('snps/'+name+'/chrIV.npz')
    coverage = data['coverage'][()]['tot']
    nmcovs.update({name:coverage.mean()})
    nlcovs.update({name:coverage})

    regpair = np.zeros((len(rpair),2))
    covpair = np.zeros((len(rpair),2))
    for i in range(len(rpair)):
        sgn = np.sign(rpair[i,0])*np.sign(rpair[i,1])
        i0 = int(abs(rpair[i,0]))
        i1 = int(abs(rpair[i,1]))
        if i0<i1:
            regpair[i] = [i0,sgn*i1]
            covpair[i] = [max(1,coverage[i0]),max(1,coverage[i1])]
        else:
            regpair[i] = [i1,sgn*i0]
            covpair[i] = [max(1,coverage[i1]),max(1,coverage[i0])]

    pset = []
    cvset = []
    cc   = []
    for i in range(len(regpair)):
        if not pset:  # empty
            pset.append(regpair[i])
            cvset.append(covpair[i])
            cc.append(1)
        else: # not empty
            dp  = pset - regpair[i]
            dp1 = dp[:,0]
            dp2 = dp[:,1]
            dpp = abs(dp1-dp2)
            #id  = np.argmin(dpp)
            id  = np.argmin(abs(dp1))
            if dpp[id]<5 and abs(dp1[id])<15:
                pset[id] = (cc[id]*pset[id]+regpair[i])/(cc[id]+1)
                cvset[id] = (cc[id]*cvset[id]+covpair[i])/(cc[id]+1)
                cc[id]  += 1
            else:
                pset.append(regpair[i])
                cvset.append(covpair[i])
                cc.append(1)
    pset = np.array(pset)
    cvset = np.array(cvset)
    cv2 = cvset[:,0]*cvset[:,1]
    c2v = cvset[:,0]+cvset[:,1]
    cc = np.array(cc)

    npairs.update({name:pset})
    ncovs.update({name:cvset})
    nc2vs.update({name:c2v})
    ncov2s.update({name:cv2})
    ncounts.update({name:cc})
    nrats.update({name: ncounts[name]/ncov2s[name]})

pset = []
cc   = []
for name in names:
    for i,pp in enumerate(pairs[name]):
        if not pset:  # empty
            pset.append(pp)
            cc.append(counts[name][i])
        else: # not empty
            dp  = pset - pp
            dp1 = dp[:,0]
            dp2 = dp[:,1]
            dpp = abs(dp1-dp2)
            id  = np.argmin(abs(dp1))
            if dpp[id]<5 and abs(dp1[id])<15:
                pset[id] = (cc[id]*pset[id]+counts[name][i]*pp)/(cc[id]+counts[name][i])
                cc[id]  += counts[name][i]
            else:
                pset.append(pp)
                cc.append(counts[name][i])
pset = np.array(pset)
cc = np.array(cc)

(pset, indx) = np.unique(regpair,axis=0,return_inverse=True)
cc = np.bincount(indx,np.arange(pset))
rpns = pset[cc>1]

for i in range(5):
    fig = plt.figure()
    #id=counts[name]>0 #np.logical_or(np.logical_and(rats[name]>0.001,counts[name]>1),rats[name]>0.01)
    fig=plotCircle(fig,rpair[i+1],totlen=mt_len,alpha=min(0.2,100/len(rpair[i+1]))) #,wt=counts[name][id]
    fig=featureCircle(fig,features,totlen=mt_len)
    fig=GCCircle(fig,GCdict,totlen=mt_len)
    #fig=plotcoverage(fig,lcovs[name])
    fig.gca().axis('off')
    plt.savefig('Recombmap/recombs_barcode%02d_W303_cov_single.pdf' % (i+1))
    plt.close(fig)

for name in names:
    fig = plt.figure()
    #id=counts[name]>0 #np.logical_or(np.logical_and(rats[name]>0.001,counts[name]>1),rats[name]>0.01)
    fig=plotCircle(fig,rpair[name][plen[name]>1000,:],totlen=mt_len,alpha=min(0.2,100/np.sum(plen[name]>1000))) #,wt=counts[name][id]
    fig=featureCircle(fig,features,totlen=mt_len)
    fig=GCCircle(fig,GCdict,totlen=mt_len)
    #fig=plotcoverage(fig,lcovs[name])
    fig.gca().axis('off')
    plt.savefig('Recombmap/paired_'+name+'_S288c_1000.pdf')
    plt.close(fig)

for name in names[14:]:
    data = np.load('DelATP21/data/'+name+'_paired_W303_chrM_mito0.npz')
    mtcov = data['ccov'][()]['chrM']
    reccov = {}
    for i in np.arange(len(plen[name])):
        if plen[name][i]>1000:
            sgn = np.sign(rpair[name][i,0])*np.sign(rpair[name][i,1])
            if sgn not in reccov:
                reccov.update({sgn:[mtcov[int(abs(rpair[name][i,0]))],mtcov[int(abs(rpair[name][i,1]))]]})
            else:
                reccov[sgn].append(mtcov[int(abs(rpair[name][i,0]))])
                reccov[sgn].append(mtcov[int(abs(rpair[name][i,1]))])
    fig = plt.figure()
    (y,x) = np.histogram(reccov[-1],bins=20)
    plt.plot((x[:-1]+x[1:])/2/np.mean(mtcov),y/y.sum()/(x[1:]-x[:-1])*np.mean(mtcov),'C0-',label='reverse')
    (y,x) = np.histogram(reccov[1],bins=20)
    plt.plot((x[:-1]+x[1:])/2/np.mean(mtcov),y/y.sum()/(x[1:]-x[:-1])*np.mean(mtcov),'C1-',label='forward')
    (y,x) = np.histogram(mtcov,bins=20)
    plt.plot((x[:-1]+x[1:])/2/np.mean(mtcov),y/y.sum()/(x[1:]-x[:-1])*np.mean(mtcov),'k--',label='whole genome')
    #a = np.var(mtcov)/np.mean(mtcov)**2+1;
    #plt.plot([a,a],[0,max(y/y.sum()/(x[1:]-x[:-1])*np.mean(mtcov)**2)*1.5],'k--')
    plt.xlabel('depth of coverage / mt coverage')
    plt.ylabel('pdf')
    plt.legend(loc='best')
    plt.savefig('recomb_coverage_pdf_'+name+'_10.pdf')
    plt.close(fig)


(y,x) = np.histogram(sc1[noself],bins=np.arange(.5/30,1.02,1/30))
ytot = y.sum()
plt.plot((x[:-1]+x[1:])/2,y/y.sum()/pr,'-sC0')
(y,x) = np.histogram(sc1[np.intersect1d(np.where(idp>0)[0],noself)],bins=np.arange(.5/30,1.02,1/30))
plt.plot((x[:-1]+x[1:])/2,y/y.sum()/p,'-vC1')
(y,x) = np.histogram(sc2[np.where(np.isin(ks,noself))],bins=np.arange(.5/30,1.02,1./30))
plt.plot((x[:-1]+x[1:])/2,y/ytot/p1r,'-oC3')
(y,x) = np.histogram(sc2[np.where(np.logical_and(np.isin(ks,noself),np.logical_not(np.isin(ks,oid))))],bins=np.arange(.5/30,1.02,1./30))
plt.plot((x[:-1]+x[1:])/2,y/ytot/p1r,'->C2')
plt.yscale('log')
xx = np.arange(0.,1.01,0.01)
plt.plot(xx,np.exp(12.1*xx-6.),'--k')
xx = np.arange(0.6,1.01,0.01)
plt.plot(xx,np.exp(30.7*xx-19.0),'-k')
xx = np.arange(0.75,1.01,0.01)
plt.plot(xx,np.exp(30.7*xx-19.0)/15,'--k')
plt.xlabel('similarity')
plt.ylabel('likelihood to recombine')
plt.legend(['total','forward pairs','non GC','non GC and non rep/orig',r'$y=\exp(12.1x-6.0)$',r'$y=\exp(30.7x-19.0)$'],loc='best')
plt.savefig('precombsimi_1914_noGCrepori_fit_1.pdf')

ngc = np.zeros(2*len(sc1))
for i in range(len(rpair)):
    if di1[i]>30 and sc1[i]>0.8:
        sf = sf1[i]
        i0 = int(abs(rpair[i][0]))
        i1 = int(abs(rpair[i][1]))
        l00 = i0+sf-15
        l01 = i0+sf+15
        l10 = i1+sf-15
        l11 = i1+sf+15
        sq0 = seqpiece(mt_chrom,l00,l01)
        sq1 = seqpiece(mt_chrom,l10,l11)
        ngc[2*i] = ATfrac(sq0)
        ngc[2*i+1] = ATfrac(sq1)

ngc = np.zeros(2*len(sc3))
ii=0
for i in range(len(rpair)):
    if i in id3:
        if di3[ii]>30 and sc3[ii]>0.8:
            sf = sf3[ii]
            i0 = int(abs(rpair[i][0]))
            i1 = int(abs(rpair[i][1]))
            l00 = i0+sf-15
            l01 = i0+sf+15
            l10 = i1+sf-15
            l11 = i1+sf+15
            sq0 = seqpiece(mt_chrom,l00,l01)
            sq1 = seqpiece(mt_chrom,l10,l11)
            ngc[2*ii] = ATfrac(sq0)
            ngc[2*ii+1] = ATfrac(sq1)
        ii += 1


ii = []
for i in rid:
    read = potreads1[i]
    idx  = readid[read.query_name]
    end1 = mtreads[idx[0]]
    end2 = mtreads[idx[1]]
    t_start = time.time()
    i1, i2 = findrecomb(read, end1, end2)
    t_end = time.time()
    print (t_end-t_start), i1, i2
    ii.append([i1,i2])
#i1, i2 = findrecomb(read, end1, end2)
cc = 0
for i in range(len(rid)):
    [i1,i2] = ii[i]
    read = potreads1[rid[i]]
    idx = readid[read.query_name]
    end1 = mtreads[idx[0]]
    end2 = mtreads[idx[1]]
    if read.seq[0:choplen]==end2.seq or read.seq[0:choplen]==str(Seq(end2.seq).reverse_complement()):
        etmp = end1
        end1 = end2
        end2 = etmp
    if i1>0 and i2<0:
        if read.seq[0:choplen]==end1.seq:
            start = end1.reference_start
        else:
            start = end1.reference_end-stepsize*i1
        end = start+stepsize*i1
        if (start<0 and end<=0) or (start>=0 and end<=mt_len):
            endcov[start:end] += 1
        else:
            if start<0 and end>0:
                endcov[start:] += 1
                endcov[:end] += 1
            else:
                endcov[(start-mt_len):(end-mt_len)] +=1
    if i2>0 and i1<0:
        if read.seq[0:choplen]==end2.seq:
            start = end2.reference_end-stepsize*i2
        else:
            start = end1.reference_start+stepsize*(i2-1)
        end = start+stepsize
        if (start<0 and end<=0) or (start>=0 and end<=mt_len):
            endcov[start:end] += 1
        else:
            if start<0 and end>0:
                endcov[start:] += 1
                endcov[:end] += 1
            else:
                endcov[(start-mt_len):(end-mt_len)] +=1

seq = read.seq
start1 = read.reference_start
start2 = end1.reference_end
i1 = forwardloc(seq,start1,start2,rflg=True)
dp = (i1-1)*stepsize
p1 = dp
while stepsize > min_step:
    seq = seq[dp:]
    start1 += dp
    start2 -= dp
    stepsize = int(stepsize/10)
    i1 = forwardloc(seq,start1,start2,rflg=True)
    dp = (i1-1)*stepsize
    p1 += dp

res = 450  # resolution
L   = int(round(mt_len/res))
x, y  = np.meshgrid(np.arange(L)*res, np.arange(L)*res)
#kx,ky = np.meshgrid(np.arange(L)*res/mt_len, np.arange(L)*res/mt_len)

recmat = {}
for name in names:
    mat = np.zeros((L,L))
    #Gmat = np.exp(-0.5*res**2*((kx-np.arange(L)*res/mt_len/2)**2+(ky-np.arange(L)*res/mt_len/2)**2)*(2*np.pi)**2)
    #rmat = np.zeros((L,L))
    pairs = rpair[name][plen[name]>1000,:]
    sgn = np.sign(pairs[:,1])
    for i in np.arange(len(sgn)):
        #rmat += sgn[i]
        mat += sgn[i]*np.exp(-0.5*((x-pairs[i,0])**2+(y-abs(pairs[i,1]))**2)/res**2)
        mat += sgn[i]*np.exp(-0.5*((y-pairs[i,0])**2+(x-abs(pairs[i,1]))**2)/res**2)
    recmat.update({name:mat})

dc = {}
for name in names:
    rp = rpair[name]
    pl = plen[name]
    dcov = np.zeros((mt_len,))
    for v in rp[pl>1000]:
        if v[0]*v[1] > 0:
            v1 = int(v[0])
            v2 = int(v[1])
            if v1>v2 and v1>0:
                dcov[v1:] += 1
                dcov[:v2] += 1
            elif v1<v2 and v1<0:
                dcov[v1:v2] += 1
            elif v1>v2 and v1<0:
                dcov[v2:v1] += 1
            else:
                dcov[v2:] += 1
                dcov[:v1] += 1
    dc.update({name:dcov})

sc = recmat['T7_W303'][59,99]
for name in names:
    fig = plt.figure()
    plt.imshow(recmat[name]/abs(recmat[name][59,99])-recmat['T7_W303']/abs(sc))
    plt.colorbar()
    plt.xticks(np.arange(0,200,20), (np.arange(0,200,20)*450/1000).astype(int))
    plt.yticks(np.arange(0,200,20), (np.arange(0,200,20)*450/1000).astype(int))
    plt.xlabel('location (kbp)')
    plt.ylabel('location (kbp)')
    plt.savefig('recombmat_'+name+'_W303.pdf')
    plt.close(fig)

"""