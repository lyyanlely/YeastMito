"""
Created on Fri Jan 18, 2019

@author: Le Yan

Functions to compute the complexity of the sequences.

Adapted from Perl code created by Juan Caballero


  # Example 1. Run all methods for each sequence
  use SeqComplex;

  my @seq = loadSeqs($file);
  my $num = 0;
  foreach my $seq (@seq) {
	$num++;
	my %results = Complex::runAllMethods( $seq, $win, $kmer );
	foreach my $m (keys %results) {
		my $res = join "\t", @{ $results{$m} };
		print "seq_$num\t$m\t$res\n";
	}
  }


Calculate composition and complexity of a DNA sequence.

:methods exports methods  = [gc at gcs ats cpg cwf ce cz ct cl cm]
:utils   exports basic utils = [log_k pot createWords countWords randSeq]
:all     exports @methods and @utils

"""

from __future__ import division
import numpy as np
from numpy.fft import fft, ifft
import gzip
import os
import subprocess
from fnmatch import filter
from natsort import natsorted
from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib.pyplot as plt
from FWHT import FWHT
from printSeq import Seq2Bin

eps = 1e-15
winsize = 256 # close to optimal

methods = ['gc','gcs','ats','she','lsc','mkc']  #,'wfc','czip'

"""
vcov = np.var(cov_nano)
for method in methods:
  fig = plt.figure()
  for slen in 50*np.arange(3,11):
    mt_comp = SeqComplex(mt_chrom, slen, periodic=1, method=method)
    vcomp   = np.var(mt_comp)
    mt_corr = periodic_corr(cov_nano-np.mean(cov_nano),mt_comp-np.mean(mt_comp))
    fig.gca().plot(np.arange(mt_len),mt_corr/mt_len/np.sqrt(vcov*vcomp),label='window size=%d' % slen)
  fig.gca().set_xscale('log')
  fig.gca().set_xlim([0.5,50000])
  fig.gca().legend(loc='best')
  fig.gca().set_title('complexity of mt computed by %s' % method)
  fig.gca().set_xlabel('distance (bp)')
  fig.gca().set_ylabel('correlation')
  fig.savefig('correlation_complex_%s_winsizes_nano.pdf' % method)
  plt.close(fig)

FeatCovMat = np.zeros((mt_len,2*winsize+1))
for i in np.arange(mt_len):
  if i>mt_len-winsize:
    seq = np.concatenate((mt_chrom[i:],mt_chrom[:i+winsize-mt_len]))
  else:
    seq = mt_chrom[i:i+winsize]
  bseq = Seq2Bin(seq)
  feats = FWHT(bseq).T
  FeatCovMat[i] = np.concatenate((feats[0],[cov[(i+int(winsize/2))%mt_len]]))

FeatMat = np.zeros((mt_len,2*winsize))
for i in np.arange(mt_len):
  if i>mt_len-winsize:
    seq = np.concatenate((mt_chrom[i:],mt_chrom[:i+winsize-mt_len]))
  else:
    seq = mt_chrom[i:i+winsize]
  bseq = Seq2Bin(seq)
  feats = FWHT(bseq).T
  FeatMat[i] = feats[0]
"""
#
def SeqComplex(seq, winsize, wordsize=5, periodic=0, method='lsc'):  #gc composition explain most coverage variance up to 70% when winsize~300
  if method == 'gc':
    return gc(seq, winsize, periodic)
  elif method == 'gcs':
    return gcs(seq, winsize, periodic)
  if method == 'at':
    return at(seq, winsize, periodic)
  elif method == 'ats':
    return ats(seq, winsize, periodic)
  elif method == 'she':
    return she(seq, winsize, periodic)
  #elif method == 'wfc':
  #  return wfc(seq, winsize, periodic)
  elif method == 'czip':
    return czip(seq, winsize, periodic)
  elif method == 'lsc':
    return lsc(seq, winsize, wordsize, periodic)
  elif method == 'mkc':
    return mkc(seq, winsize, wordsize, periodic)
  else:
    print('Please enter a valid method!')
    return None


def gc(seq, win, periodic):
  """
   Function to calculate the GC content.
   (|G|+|C|)/win
  """
  if type(seq) == str:
    seq = np.fromstring(seq, dtype='S1')
  slen = len(seq)
  winhalf = int(win/2)

  gcid = np.logical_or(seq=='G', seq=='C')

  if periodic:
    gcf = np.convolve(np.concatenate((gcid[-winhalf:],gcid,gcid[:winhalf])),np.ones(win)/win)[win-1:slen+win-1]
  else:
    gcf = np.convolve(gcid,np.ones(win)/win,'same')
  return gcf


def at(seq, win, periodic):
  """
   Function to calculate the AT content.
   (|A|+|T|)/win
  """
  if type(seq) == str:
    seq = np.fromstring(seq, dtype='S1')
  slen = len(seq)
  winhalf = int(win/2)

  atid = np.logical_or(seq=='A', seq=='T')

  if periodic:
    atf = np.convolve(np.concatenate((atid[-winhalf:],atid,atid[:winhalf])),np.ones(win)/win)[win-1:slen+win-1]
  else:
    atf = np.convolve(atid,np.ones(win)/win,'same')
  return atf

def gcs(seq, win, periodic):
  """
   Function to calculate the GC skew.
   (|G|-|C|)/(|G|+|C|)
  """
  if type(seq) == str:
    seq = np.fromstring(seq, dtype='S1')
  slen = len(seq)
  winhalf = int(win/2)

  gid  = seq=='G'
  cid  = seq=='C'
  #gcid = np.logical_or(gid, cid)

  if periodic:
    gf = np.convolve(np.concatenate((gid[-winhalf:],gid,gid[:winhalf])),np.ones(win)/win)[win-1:slen+win-1]
    cf = np.convolve(np.concatenate((cid[-winhalf:],cid,cid[:winhalf])),np.ones(win)/win)[win-1:slen+win-1]
  else:
    gf = np.convolve(gid,np.ones(win)/win,'same')
    cf = np.convolve(cid,np.ones(win)/win,'same')
  gcf = gf+cf
  gcf[np.where(gcf==0)] = 1/win
  return (gf-cf)/gcf

def ats(seq, win, periodic):
  """
   Function to calculate the AT skew.
   (|A|-|T|)/(|A|+|T|)
  """
  if type(seq) == str:
    seq = np.fromstring(seq, dtype='S1')
  slen = len(seq)
  winhalf = int(win/2)

  aid  = seq=='A'
  tid  = seq=='T'
  #atid = np.logical_or(aid, tid)

  if periodic:
    af = np.convolve(np.concatenate((aid[-winhalf:],aid,aid[:winhalf])),np.ones(win)/win)[win-1:slen+win-1]
    tf = np.convolve(np.concatenate((tid[-winhalf:],tid,tid[:winhalf])),np.ones(win)/win)[win-1:slen+win-1]
  else:
    af = np.convolve(aid,np.ones(win)/win,'same')
    tf = np.convolve(tid,np.ones(win)/win,'same')
  atf = af+tf
  atf[np.where(atf==0)] = 1/win
  return (af-tf)/atf

def she(seq, win, periodic):
  """
   Function to calculate the Shannon Entropy.
   S = - Sum_i P_i*log_2(P_i)
  """
  if type(seq) == str:
    seq = np.fromstring(seq, dtype='S1')
  slen = len(seq)
  winhalf = int(win/2)

  aid  = seq=='A'
  tid  = seq=='T'
  gid  = seq=='G'
  cid  = seq=='C'

  if periodic:
    af = np.convolve(np.concatenate((aid[-winhalf:],aid,aid[:winhalf])),np.ones(win)/win)[win-1:slen+win-1]
    tf = np.convolve(np.concatenate((tid[-winhalf:],tid,tid[:winhalf])),np.ones(win)/win)[win-1:slen+win-1]
    gf = np.convolve(np.concatenate((gid[-winhalf:],gid,gid[:winhalf])),np.ones(win)/win)[win-1:slen+win-1]
    cf = np.convolve(np.concatenate((cid[-winhalf:],cid,cid[:winhalf])),np.ones(win)/win)[win-1:slen+win-1]
  else:
    af = np.convolve(aid,np.ones(win)/win,'same')
    tf = np.convolve(tid,np.ones(win)/win,'same')
    gf = np.convolve(gid,np.ones(win)/win,'same')
    cf = np.convolve(cid,np.ones(win)/win,'same')
  
  return -af*np.log2(af+eps)-tf*np.log2(tf+eps)-gf*np.log2(gf+eps)-cf*np.log2(cf+eps)

def czip(seq, win, periodic):
  """
   Function to calculate the Compression factor.
   S = - Sum_i P_i*log_2(P_i)
  """
  if type(seq) == str:
    seq = np.fromstring(seq, dtype='S1')
  slen = len(seq)
  winhalf = int(win/2)

  print("Not complete")

  return None

def countWords(seq, l):
  words = {}

  slen = len(seq)

  for i in np.arange(slen-l+1):
    k = seq[i:i+l]
    if words.has_key(k):
      words[k] += 1
    else:
      words.update({k:1})

  return words

def lsc(seq, win, wlen, periodic):
  """
   Function to calculate the linguistic sequence complexity.

   https://en.wikipedia.org/wiki/Linguistic_sequence_complexity
  """
  if type(seq) != str:
    seq = seq.tostring()

  slen = len(seq)
  winhalf = int(win/2)

  if periodic:
    seq = seq[-winhalf:]+seq+seq[:winhalf]
    sl = slen
  else:
    sl = slen-win+1
  ls = np.zeros(sl)

  for i in np.arange(sl):
    vl = 1
    vm = 1
    for l in np.arange(wlen):
      words = countWords(seq[i:i+win], l+1)
      vm *= min(4**(l+1),win-l+1)
      vl *= len(words.keys())
    ls[i] = vl/vm

  return ls

def mkc(seq, win, wlen, periodic):
  """
  Function to calculate the Complexity of Markov model values.

  Call: mkc( seq, win, word ) str, int, int

  Y.L. Orlov, and V.N. Potapov, "Complexity: an internet resource for analysis of DNA
  sequence complexity", 1994 NAR

  Returns: numpy.ndarray

  """
  if type(seq) != str:
    seq = seq.tostring()

  slen = len(seq)
  winhalf = int(win/2)

  m = 4**wlen

  if periodic:
    seq = seq[-winhalf:]+seq+seq[:winhalf]
    sl = slen
  else:
    sl = slen-win+1
  mk = np.zeros(sl)

  for i in np.arange(sl):
    words = countWords(seq[i:i+win], wlen)
    
    for k,v in words.items():
      r = v/(win-wlen+1)
      mk[i] -= r*np.log(r)/np.log(m)

  return mk

def periodic_corr(x, y):
    """Periodic correlation, implemented using the FFT.

    x and y must be real sequences with the same length.
    """
    return ifft(fft(x) * fft(y).conj()).real