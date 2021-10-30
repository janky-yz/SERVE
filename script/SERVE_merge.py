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


def Taco_merge(args):
    cmd = 'taco_run'+' '+args.input_gtf_list \
		+' -p '+str(args.nthread) \
		+' --gtf-expr-attr '+'TPM' \
		+' -o '+merge_dir
    return cmd

def Cuff_merge(args):
    cmd = 'cuffmerge'+' -o '+merge_dir \
		+' -p '+str(args.nthread) \
		+' '+args.input_gtf_list
    return cmd

def StringTie_merge(args):
    cmd = 'stringtie --merge'+' '+args.input_gtf_list \
                +' -o '+assem_gtf

    return cmd

def Compare(args):
    cmd = 'cuffcompare'+' -T ' \
                +' -o '+args.prefix+'_cuff' \
                +' -r '+ERV_merge_gtf \
                +' -i '+args.input_gtf_list

    return cmd

def create_GMAP_index(args):
    cmd = 'gmap_build'+' -D '+args.GMAP_index \
                      +' -d '+args.GMAP_index_name \
                      +' '+args.ref_genome

    return cmd

def Remap(args):
    cmd = 'gmap'+' -t '+str(args.nthread) \
		+' -D '+args.GMAP_index \
		+' -d '+args.GMAP_index_name \
		+' -f gff3_gene' \
		+' '+ERV_fasta \
		+' > '+ERV_gff3

    return cmd

parser = argparse.ArgumentParser(description='SERVE_merge: merge expressed ERVs')
parser.add_argument('-i', '--input_gtf_list', help='A text file with a list of SERVE GTF files (required)')
parser.add_argument('-p', '--prefix', default='SERVE', help='Prefix for output file name (default: SERVE)')
parser.add_argument('--taco', help='Merge ERV transcripts with TACO (default: FALSE)')
parser.add_argument('--stringtie', help='Merge ERV transcripts with StringTie (default: FALSE)')
parser.add_argument('--cuffmerge', help='Merge ERV transcripts with cuffmerge (default: FALSE)')
parser.add_argument('-n', '--nsample', type=int, help='The number of samples included in the input sample list (required)')
parser.add_argument('-r', '--ref_genome', help='Reference genome in FASTA format (required)')
parser.add_argument('-G', '--GMAP_index', default='./GMAP_index', help='Path to the directory where GMAP index generated (default: GMAP_index)')
parser.add_argument('-g', '--GMAP_index_name', default='GRCh38', help='GMAP index name (default: GRCh38)')
parser.add_argument('-t', '--nthread', type=int, default=1, help='Number of threads to run SERVE_merge (default: 1)')
parser.add_argument('-l', '--length', default=200, help='Minimum ERV length (bp) (default: 200)')
parser.add_argument('--ratio', default=0.50, help='Minimum sample ratio of ERV identified (default: 0.50)')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory (default: .)')

args = parser.parse_args()
script_dir = os.path.abspath(os.path.dirname(__file__))

merge_dir= args.output_dir+'/ERV_merge'

assem_gtf = args.output_dir+'/'+args.prefix+'_assem.gtf'
ERV_merge_gtf = args.output_dir+'/'+args.prefix+'_ERV_merge.gtf'
ERV_fasta = args.output_dir+'/'+args.prefix+'_gmap.fasta'
ERV_gff3 = args.output_dir+'/'+args.prefix+'_gmap.gff3'
ERV_gtf = args.output_dir+'/'+args.prefix+'_gmap.gtf'

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)


print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Running SERVE_merge on {0:d} threads.'.format(args.nthread), flush=True)


with cd(args.output_dir):

	print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Start to merge ERV transcripts.', flush=True)

	if args.taco:
		cmd = Taco_merge(args)
		subprocess.check_call(cmd, shell=True, executable='/bin/bash')
		shutil.copyfile('./ERV_merge/assembly.gtf',assem_gtf)
		shutil.rmtree(merge_dir)

	elif args.stringtie:
		cmd = StringTie_merge(args)
		subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	elif args.cuffmerge:
		cmd = Cuff_merge(args)
		subprocess.check_call(cmd, shell=True, executable='/bin/bash')
		shutil.copyfile('./ERV_merge/merged.gtf',assem_gtf)
		shutil.rmtree(merge_dir)

	elif args.nsample>50:
		cmd = Taco_merge(args)
		subprocess.check_call(cmd, shell=True, executable='/bin/bash')
		shutil.copyfile('./ERV_merge/assembly.gtf',assem_gtf)
		shutil.rmtree(merge_dir)

	else:
		cmd = StringTie_merge(args)
		subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finish.', flush=True)


	print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Start to quality control ERVs.', flush=True)

	cmd = 'gffread'+' -T --sort-alpha'+' -o '+ERV_merge_gtf+' '+assem_gtf
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	os.remove(assem_gtf)

	cmd = Compare(args)
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	cmd = 'Rscript '+script_dir+'/SERVE_merge_QC.R 1 '+args.prefix+' '+str(args.ratio)+' '+str(args.length)
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	cmd = 'gffread'+' -g '+args.ref_genome+' -w '+ERV_fasta+' '+ERV_gtf
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	if not os.path.exists(args.GMAP_index):
		os.makedirs(args.GMAP_index)
		if not args.ref_genome:
			print('ERROR: Lack reference genome (--ref_genome)')
		print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Start to create GMAP index.', flush=True)

		cmd = create_GMAP_index(args)
		subprocess.check_call(cmd, shell=True, executable='/bin/bash')

		print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finish (In ./GMAP_index/).', flush=True)

	cmd = Remap(args)
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	cmd = 'Rscript '+script_dir+'/SERVE_merge_QC.R 2 '+args.prefix+' '+str(args.ratio)+' '+str(args.length)
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')


	print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finish.', flush=True)


print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] ERV merge is done. Now you can quantify expressed ERVs with SERVE_quant.', flush=True)
