#!/bin/bash

nthread=""
prefix=""

data_dir=""
wrk_dir=""
ref_dir=""

ref_genome=${ref_dir}/${GRCh}.fa
ref_annotation=${ref_dir}/${GRCh}.gtf
ref_erv=${ref_dir}/${ERV}.bed

fastq1=${data_dir}/${sample}_1.fastq.gz
fastq2=${data_dir}/${sample}_2.fastq.gz

cd ${wrk_dir}

## Step 1: Detecting expressed ERVs for each sample

SERVE.py --fastq1 ${fastq1} --fastq2 ${fastq2} -e ${ref_erv} -r ${ref_genome} -a ${ref_annotation} -t ${nthread} -p ${prefix}

## Step 2: Merging ERV transcripts of multiple samples in one condition
<<b
nsample=`ls ${wrk_dir}/3_qc/*gtf | wc -l`

ls ${wrk_dir}/3_qc/*gtf >${wrk_dir}/gtf.list
SERVE_merge.py -i ${wrk_dir}/gtf.list -n ${nsample} -r ${ref_genome} -t ${nthread}
b
## Step 3: Quantifying ERVs for each sample
<<b
ERV_annotation=${wrk_dir}/SERVE_ERV_merge.gtf ## produced by SERVE_merge.py
merge_annotation=${wrk_dir}/GRCh38_ERV_merge.gtf

cat ${ref_annotation} ${ERV_annotation} >${merge_annotation}
SERVE_quant.py --fastq1 ${fastq1} --fastq2 ${fastq2} -r ${ref_genome} -a ${merge_annotation} -t ${nthread} -p ${prefix}

ls ${wrk_dir}/*genes.results >${wrk_dir}/sample.list ## Merging quantification results
SERVE_quant_QC.py -i ${wrk_dir}/sample.list 
b
