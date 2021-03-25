#!/usr/bin/env python3
# Author: Janky

import argparse
import os
import subprocess
import gzip
import multiprocessing as mp
import contextlib
import shutil
from datetime import datetime

@contextlib.contextmanager
def cd(cd_path):
    saved_path = os.getcwd()
    os.chdir(cd_path)
    yield
    os.chdir(saved_path)


def create_RSEM_index(args):
    cmd = 'rsem-prepare-reference'+' --gtf '+args.annotation \
		+' --star ' \
		+' -p '+str(args.nthread) \
		+' '+args.ref_genome \
		+' '+args.RSEM_index+'/RSEM'

    return cmd

def Quantify(args):
    cmd = 'rsem-calculate-expression'+' --paired-end ' \
		+' --forward-prob '+str(args.forward_prob) \
		+' --star ' \
		+' -p '+str(args.nthread) \
		+' --star-gzipped-read-file '+args.fastq1+' '+args.fastq2 \
		+' '+args.RSEM_index+'/RSEM' \
		+' '+args.prefix
		
    return cmd


parser = argparse.ArgumentParser(description='SERVE_quant: Quantify expressed ERVs')
parser.add_argument('-fq1', '--fastq1', help='Read1 in FASTQ format (required)')
parser.add_argument('-fq2', '--fastq2', help='Read1 in FASTQ format (required)')
parser.add_argument('-p', '--prefix', default='SERVE', help='Prefix for output file name (default: SERVE)')
parser.add_argument('-r', '--ref_genome', help='Reference genome in FASTA format (required)')
parser.add_argument('-a', '--annotation', help='Genome annotation in GTF format (required)')
parser.add_argument('-R', '--RSEM_index', default='./RSEM_index', help='Path to the directory where RSEM index generated (default: RSEM_index)')
parser.add_argument('-t', '--nthread', type=int, default=1, help='Number of threads to run SERVE (default: 1)')
parser.add_argument('-f', '--forward_prob', default=0.5, help='Probability of generating a read from the forward strand of a transcript. 0.5 for unstranded-specific, 0 for stranded-specific where upstream reads are all derived from the forward strand, 1 for stranded-specific where upstream reads are all derived from the reverse strand (default: 0.5)')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory (default: .)')

args = parser.parse_args()
script_dir = os.path.abspath(os.path.dirname(__file__))


if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)


print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Running SERVE_quant on {0:d} threads.'.format(args.nthread), flush=True)


with cd(args.output_dir):
	if not os.path.exists(args.RSEM_index):
		os.makedirs(args.RSEM_index)
		if not args.annotation:
			print('ERROR: Lack annotation file (--annotation)')
		if not args.ref_genome:
			print('ERROR: Lack reference genome (--ref_genome)')

		print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Start to create RSEM index.', flush=True)

		cmd = create_RSEM_index(args)
		subprocess.check_call(cmd, shell=True, executable='/bin/bash')

		print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finish (In ./RSEM_index/).', flush=True)



	print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Start to quantify RNA-seq reads.', flush=True)

	cmd = Quantify(args)
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finish.', flush=True)

	shutil.rmtree(args.prefix+'.stat')
	os.remove(args.prefix+'.transcript.bam')
	

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] ERV quantification is done. Now you can do ERV quality control with SERVE_quant_QC among multiple samples (Recommond).', flush=True)
