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


def create_STAR_index(args):
    cmd = 'STAR'+' --runMode genomeGenerate' \
		+' --runThreadN '+str(args.nthread) \
		+' --genomeDir '+args.STAR_index \
		+' --genomeFastaFiles '+args.ref_genome \
		+' --sjdbGTFfile '+args.annotation \
		+' --sjdbOverhang 100' \
		+' --genomeSAindexNbases '+str(args.genomeSAindexNbases)
    return cmd

def STAR(args):
    cmd = 'STAR'+' --runThreadN '+str(args.nthread) \
		+' --genomeDir '+args.STAR_index \
		+' --outFileNamePrefix '+align_dir+'/'+args.prefix+'_' \
		+' --readFilesIn '+args.fastq1+' '+args.fastq2 \
		+' --readFilesCommand zcat' \
		+' --outFilterType BySJout' \
		+' --outFilterIntronMotifs RemoveNoncanonical' \
		+' --outSAMtype BAM SortedByCoordinate' \
        	+' --outSAMattributes NH HI AS nM NM' \
		+' --twopassMode Basic' \
		+' --outSAMstrandField intronMotif'
		
    if args.nthreadsort:
        cmd += ' --outBAMsortingThreadN '+str(args.nthreadsort)
    if args.nRAMsort:
        cmd += ' --limitBAMsortRAM '+str(args.nRAMsort)
    return cmd

def create_BAM_index(args):
    cmd = 'sambamba index'+' -t '+str(args.nthread) \
                          +' '+align_bam

    return cmd

def Extract(args):
    cmd = 'sambamba view'+' -t '+str(args.nthread) \
                         +' -L '+args.erv_bed \
                         +' -f bam' \
                         +' -o '+ERV_bam \
                         +' '+align_bam

    return cmd

def Assemble(args):
    cmd = 'Trinity'+' --genome_guided_bam '+ERV_bam \
                   +' --genome_guided_max_intron '+str(args.max_intron) \
                   +' --CPU '+str(args.nthread) \
                   +' --max_memory '+args.nRAMassem \
                   +' --output '+assem_dir+'/trinity_ERV'

    if args.stranded_type:
        cmd += ' --SS_lib_type '+args.stranded_type

    return cmd

def Estimate(args):
    cmd = 'align_and_estimate_abundance.pl'+' --transcripts '+ERV_fasta \
		+' --seqType fq'+' --left '+args.fastq1+' --right '+args.fastq2 \
		+' --est_method RSEM'+' --aln_method bowtie2' \
		+' --gene_trans_map '+ERV_gene_map \
		+' --prep_reference' \
		+' --output_dir '+assem_dir \
		+' --thread_count '+str(args.nthread)

    if args.stranded_type:
        cmd += ' --SS_lib_type '+args.stranded_type

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
		+' --max-intronlength-middle='+str(args.max_intron) \
		+' --min-identity='+str(args.min_identity) \
		+' '+ERV_fasta \
		+' > '+ERV_gff3

    return cmd


parser = argparse.ArgumentParser(description='SERVE: pipeline for detecting expressed ERVs')
parser.add_argument('-fq1', '--fastq1', help='Read1 in FASTQ format (required)')
parser.add_argument('-fq2', '--fastq2', help='Read1 in FASTQ format (required)')
parser.add_argument('-e', '--erv_bed', help='ERV position in BED format (required)')
parser.add_argument('-p', '--prefix', default='SERVE', help='Prefix for output file name (default: SERVE)')
parser.add_argument('-r', '--ref_genome', help='Reference genome in FASTA format (required)')
parser.add_argument('-a', '--annotation', help='Genome annotation in GTF format')
parser.add_argument('--genomeSAindexNbases', default=14, type=int, help='length (bases) of the SA pre-indexing string for creating STAR index. Typically between 10 and 15. For small genomes, this parameter must be scaled down to min(14, log2(GenomeLength)/2-1)')
parser.add_argument('-S', '--STAR_index', default='./STAR_index', help='Path to the directory where STAR index generated (default: STAR_index)')
parser.add_argument('-G', '--GMAP_index', default='./GMAP_index', help='Path to the directory where GMAP index generated (default: GMAP_index)')
parser.add_argument('-g', '--GMAP_index_name', default='GRCh38', help='GMAP index name (default: GRCh38)')
parser.add_argument('-t', '--nthread', type=int, default=1, help='Number of threads to run SERVE (default: 1)')
parser.add_argument('--nthreadsort', type=int, help='Number of threads for BAM sorting')
parser.add_argument('--nRAMsort', type=int, help='Maximum available RAM (bytes) for sorting BAM.')
parser.add_argument('-s', '--stranded_type', default=None, help='Strand-specific RNA-seq read orientation: RF or FR (default: None)')
parser.add_argument('-m', '--nRAMassem', default='10G', help='Maximum available RAM (Gb) for assembly (default: 10G)')
parser.add_argument('--max_intron', default=10000, type=int, help='Maximum intron length of ERVs (default: 10000)')
parser.add_argument('--min_identity', default=0.96, help='Minimum identity of ERV transcripts (default: 0.96)')
parser.add_argument('--count', default=5, help='Minimum ERV count (default: 5)')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory (default: .)')

args = parser.parse_args()
script_dir = os.path.abspath(os.path.dirname(__file__))

align_dir='./1_align'
assem_dir='./2_assem'
qc_dir='./3_qc'

align_bam = align_dir+'/'+args.prefix+'_Aligned.sortedByCoord.out.bam'
ERV_bam = align_dir+'/'+args.prefix+'_ERV.bam'
ERV_fasta = assem_dir+'/'+args.prefix+'_ERV.fasta'
ERV_gene_map = assem_dir+'/'+args.prefix+'_ERV_gene_trans_map'
ERV_gff3 = qc_dir+'/'+args.prefix+'_ERV.gff3'
ERV_gtf = qc_dir+'/'+args.prefix+'_ERV.gtf'
ERV_bed = qc_dir+'/'+args.prefix+'_ERV.bed'
ERV_exon_bed = qc_dir+'/'+args.prefix+'_ERV_exon.bed'

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)


print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Running SERVE on {0:d} threads.'.format(args.nthread), flush=True)


with cd(args.output_dir):
	if not os.path.exists(args.STAR_index):
		os.makedirs(args.STAR_index)
		if not args.annotation:
			print('ERROR: Lack annotation file (--annotation)')
		if not args.ref_genome:
			print('ERROR: Lack reference genome (--ref_genome)')

		print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Start to create STAR index.', flush=True)

		cmd = create_STAR_index(args)
		subprocess.check_call(cmd, shell=True, executable='/bin/bash')

		print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finish (In ./STAR_index/).', flush=True)


	if not os.path.exists(align_dir):	
		os.makedirs(align_dir)

	print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Start to align RNA-seq reads to reference genome.', flush=True)

	cmd = STAR(args)
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finish (In ./1_align/).', flush=True)
	

	print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Start to extract ERV reads.', flush=True)

	cmd = create_BAM_index(args)
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	cmd = Extract(args)
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finish (In ./1_align/).', flush=True)


	if not os.path.exists(assem_dir):
		os.makedirs(assem_dir)

	print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Start to de novo assemble ERVs.', flush=True)

	cmd = Assemble(args)
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	shutil.copyfile('./2_assem/trinity_ERV/Trinity-GG.fasta',ERV_fasta)

	trinity_path = os.popen('which Trinity')
	trinity_abspath = trinity_path.read()
	cmd = os.path.dirname(trinity_abspath.splitlines()[0])+'/../opt/trinity-*/util/support_scripts/get_Trinity_gene_to_trans_map.pl ./2_assem/trinity_ERV/Trinity-GG.fasta >'+ERV_gene_map
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	cmd = Estimate(args)
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')
	
	print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finish (In ./2_assem/).', flush=True)


	if not os.path.exists(args.GMAP_index):
		os.makedirs(args.GMAP_index)
		if not args.ref_genome:
			print('ERROR: Lack reference genome (--ref_genome)')
		print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Start to create GMAP index.', flush=True)

		cmd = create_GMAP_index(args)
		subprocess.check_call(cmd, shell=True, executable='/bin/bash')

		print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finish (In ./GMAP_index/).', flush=True)


	if not os.path.exists(qc_dir):
		os.makedirs(qc_dir)

	print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Start to quality control ERVs.', flush=True)

	cmd = Remap(args)
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	cmd = 'gffread'+' -i '+str(args.max_intron)+' -g '+args.ref_genome+' -N -T --sort-alpha'+' -o '+ERV_gtf+' '+ERV_gff3
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	cmd = 'gffread'+' --bed --sort-alpha'+' -o '+ERV_bed+' '+ERV_gtf
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	cmd = "grep exon "+ERV_gtf+" | awk '{print $1,$4-1,$5,$10,100,$7}' OFS='\t' | tr -d ';' >"+ERV_exon_bed
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	ERV_overlap = qc_dir+'/'+args.prefix+'_overlap.txt'
	if args.stranded_type:
		cmd = 'bedtools intersect'+' -a '+ERV_exon_bed+' -b '+args.erv_bed+' -s -wa -wb >'+ERV_overlap
	else:
		cmd = 'bedtools intersect'+' -a '+ERV_exon_bed+' -b '+args.erv_bed+' -wa -wb >'+ERV_overlap
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')
	os.remove(ERV_exon_bed)

	cmd = 'Rscript '+script_dir+'/SERVE_QC.R '+args.prefix+' '+str(args.count)
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')

	cmd = 'gffread'+' --in-bed --sort-alpha -F'+' -o '+ERV_gff3+' '+ERV_bed
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')
	cmd = 'gffread'+' -T -F'+' -o '+ERV_gtf+' '+ERV_gff3
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')
	os.remove(ERV_overlap)

	print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finish (In ./3_qc/).', flush=True)

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] ERV identification is done. Now you can do ERV merge with SERVE_merge (Recommond).', flush=True)
