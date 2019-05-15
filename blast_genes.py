#!/usr/bin/env python

from Bio.Blast.Applications import NcbiblastnCommandline
from StringIO import StringIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import re
import string
import os

nmis = 5
nlost= 10

bad_chars = '(){}<>'
bad_chard = '&#@$&^%'
gene_dir = 'YeastGenes/'

#genefiles = [gname.split('.')[0] for gname in os.listdir(gene_dir)]

lost_genes = {}
for read in SeqIO.parse('orf_coding_all.fasta','fasta'):
    gene_name = read.description.split(' ')[1]
    gene_name = gene_name.translate(string.maketrans("", "", ), bad_chars)
    gene_seq = read.seq
    # if gene_name not in genefiles:

    seq1 = SeqRecord(gene_seq, id=gene_name, name=read.name, description=read.description.translate(string.maketrans("", "", ), bad_chard))
    gene = gene_dir+gene_name+".fasta"
    SeqIO.write(seq1, gene, "fasta")
    # Run BLAST and parse the output as XML
    output = NcbiblastnCommandline(query=gene, subject="LacO_T0.fasta", dust="no", soft_masking="false", outfmt=5)()[0]   # -dust no -soft_masking false -outmft 6
    blast_result_record = NCBIXML.read(StringIO(output))

    # Print some information on the result
    if len(blast_result_record.alignments)==0:   # no significant alignments
        lost_genes.update({gene_name: {'LacO':'Lost'}})
    else:
        tot_len = 0
        for alignment in blast_result_record.alignments:
            for hsp in alignment.hsps:
                tot_len += hsp.align_length
                if hsp.align_length-hsp.identities>nmis:
                    if gene_name not in lost_genes:
                        lost_genes.update({gene_name: {'LacO':[(hsp.align_length-hsp.identities, hsp.align_length)]}})
                    else:
                        lost_genes[gene_name]['LacO'].append((hsp.align_length-hsp.identities, hsp.align_length))
            if len(gene_seq)-tot_len>nlost:
                if gene_name not in lost_genes:
                    lost_genes.update({gene_name: {'LacO':[(len(gene_seq)-tot_len, len(gene_seq))]}})
                else:
                    lost_genes[gene_name]['LacO'].append((len(gene_seq)-tot_len, len(gene_seq)))
                # print '****Alignment****'
                # print 'sequence:', alignment.title
                # print 'length:', alignment.length
                # print 'e value:', hsp.expect
                # print hsp.query
                # print hsp.match
                # print hsp.sbjct

    # Run BLAST and parse the output as XML
    output = NcbiblastnCommandline(query=gene, subject="W303_T7.fasta", dust="no", soft_masking="false", outfmt=5)()[0]  # outfmt=5
    blast_result_record = NCBIXML.read(StringIO(output))

    # Print some information on the result
    if len(blast_result_record.alignments)==0:   # no significant alignments
        if gene_name not in lost_genes:
            lost_genes.update({gene_name: {'W303':'Lost'}})
        else:
            lost_genes[gene_name].update({'W303':'Lost'})
    else:
        tot_len = 0
        for alignment in blast_result_record.alignments:
            for hsp in alignment.hsps:
                tot_len += hsp.align_length
                if hsp.align_length-hsp.identities>nmis:
                    if gene_name not in lost_genes:
                        lost_genes.update({gene_name: {'W303':[(hsp.align_length-hsp.identities, hsp.align_length)]}})
                    else:
                        if 'W303' not in lost_genes[gene_name]:
                            lost_genes[gene_name].update({'W303':[(hsp.align_length-hsp.identities, hsp.align_length)]})
                        else:
                            lost_genes[gene_name]['W303'].append((hsp.align_length-hsp.identities, hsp.align_length))
            if len(gene_seq)-tot_len>nlost:
                if gene_name not in lost_genes:
                    lost_genes.update({gene_name: {'W303':[(len(gene_seq)-tot_len, len(gene_seq))]}})
                else:
                    if 'W303' not in lost_genes[gene_name]:
                        lost_genes[gene_name].update({'W303':[(len(gene_seq)-tot_len, len(gene_seq))]})
                    else:
                        lost_genes[gene_name]['W303'].append((len(gene_seq)-tot_len, len(gene_seq)))

    # return lost_genes

