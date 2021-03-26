# SERVE

A developed pipeline for expressed endogenous retrovirus identification.

This repository contains SERVE pipeline ([script/](https://github.com/janky-yz/SERVE/blob/main/script)) and ERV annotaion files ([Annotation/](https://github.com/janky-yz/SERVE/blob/main/Annotation))

## Software requirements

1. Python 3.x.x and R
2. STAR, sambamba, Trinity, RSEM, Bowtie2, GMAP, gffread and BEDTools (for [SERVE.py](https://github.com/janky-yz/SERVE/blob/main/script/SERVE.py))
3. TACO, StringTie, Cufflinks, GMAP, and gffread (for [SERVE_merge.py](https://github.com/janky-yz/SERVE/blob/main/script/SERVE_merge.py))
4. RSEM and STAR (for [SERVE_quant.py](https://github.com/janky-yz/SERVE/blob/main/script/SERVE_quant.py))

## Obtaining SERVE

```
git clone https://github.com/yangence/SERVE.git
chmod u+x SERVE/script/SERVE*
```

## Usage

SERVE is composed of transcript assembly ([SERVE.py](https://github.com/janky-yz/SERVE/blob/main/script/SERVE.py)), annotation merge ([SERVE_merge.py](https://github.com/janky-yz/SERVE/blob/main/script/SERVE_merge.py)) and gene quantification ([SERVE_quant.py](https://github.com/janky-yz/SERVE/blob/main/script/SERVE_quant.py) + [SERVE_quant_QC.py](https://github.com/janky-yz/SERVE/blob/main/script/SERVE_quant_QC.py)).

### Input files

```
RNA sequencing files: paired end, in FASTQ and compressed format (fastq.gz)
ERV reference bed: ERV annotation in BED format (at least 6 fields: chrom, chromStart, chromEnd, name, score, and strand). See Annotation/HERV_noncoding.bed
Reference genome: reference genome sequnce (FASTA) is required. You can create genome index before running SERVE, or provide Reference annotation (GTF) for SERVE to automatically create index.
```

!!!NOTICE: all directories and files should be provided with absolute paths.

### SERVE.py

SERVE.py is designed for ERV transcript assembly in single sample.

```bash

Example:
1. SERVE.py -fq1 ${Dir}/test_1.fastq.gz -fq2 ${Dir}/test_2.fastq.gz -e ${Dir}/ERV_ref.bed -p test -r ${Dir}/GRCh38.fa -a ${Dir}/GRCh38.gtf -t 16
2. SERVE.py -fq1 ${Dir}/test_1.fastq.gz -fq2 ${Dir}/test_2.fastq.gz -e ${Dir}/ERV_ref.bed -p test -S ${STAR_index_Dir} -G ${GMAP_index_Dir} -g GMAP -t 16
```

```bash

usage: SERVE.py [-h] [-fq1 FASTQ1] [-fq2 FASTQ2] [-e ERV_BED] [-p PREFIX]
                [-r REF_GENOME] [-a ANNOTATION]
                [--genomeSAindexNbases GENOMESAINDEXNBASES] [-S STAR_INDEX]
                [-G GMAP_INDEX] [-g GMAP_INDEX_NAME] [-t NTHREAD]
                [--nthreadsort NTHREADSORT] [--nRAMsort NRAMSORT]
                [-s STRANDED_TYPE] [-m NRAMASSEM] [--max_intron MAX_INTRON]
                [--min_identity MIN_IDENTITY] [--count COUNT] [-o OUTPUT_DIR]

SERVE: pipeline for detecting expressed ERVs

optional arguments:
  -h, --help            show this help message and exit
  -fq1 FASTQ1, --fastq1 FASTQ1
                        Read1 in FASTQ format (required)
  -fq2 FASTQ2, --fastq2 FASTQ2
                        Read1 in FASTQ format (required)
  -e ERV_BED, --erv_bed ERV_BED
                        ERV position in BED format (required)
  -p PREFIX, --prefix PREFIX
                        Prefix for output file name (default: SERVE)
  -r REF_GENOME, --ref_genome REF_GENOME
                        Reference genome in FASTA format (required)
  -a ANNOTATION, --annotation ANNOTATION
                        Genome annotation in GTF format
  --genomeSAindexNbases GENOMESAINDEXNBASES
                        length (bases) of the SA pre-indexing string for
                        creating STAR index. Typically between 10 and 15. For
                        small genomes, this parameter must be scaled down to
                        min(14, log2(GenomeLength)/2-1)
  -S STAR_INDEX, --STAR_index STAR_INDEX
                        Path to the directory where STAR index generated
                        (default: STAR_index)
  -G GMAP_INDEX, --GMAP_index GMAP_INDEX
                        Path to the directory where GMAP index generated
                        (default: GMAP_index)
  -g GMAP_INDEX_NAME, --GMAP_index_name GMAP_INDEX_NAME
                        GMAP index name (default: GRCh38)
  -t NTHREAD, --nthread NTHREAD
                        Number of threads to run SERVE (default: 1)
  --nthreadsort NTHREADSORT
                        Number of threads for BAM sorting
  --nRAMsort NRAMSORT   Maximum available RAM (bytes) for sorting BAM.
  -s STRANDED_TYPE, --stranded_type STRANDED_TYPE
                        Strand-specific RNA-seq read orientation: RF or FR
                        (default: None)
  -m NRAMASSEM, --nRAMassem NRAMASSEM
                        Maximum available RAM (Gb) for assembly (default: 10G)
  --max_intron MAX_INTRON
                        Maximum intron length of ERVs (default: 10000)
  --min_identity MIN_IDENTITY
                        Minimum identity of ERV transcripts (default: 0.96)
  --count COUNT         Minimum ERV count (default: 5)
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory (default: .)
```

In the output directory, you will see output files in three directories: "1_align/" for alignment, "2_assem/" for assembly, and "3_qc/" for quality control.

${prefix}_ERV.gtf (in "3_qc/") records ERV information in one sample. It can be used for downstream analysis (SERVE_merge.py).

### SERVE_merge.py

SERVE_merge.py is designed for annotation merge, which is suited for single sample or multiple samples in one condition.

```bash

Example:
cat *gtf >gtf.list
SERVE_merge.py -i ${Dir}/gtf.list -p test -n 10 -r ${Dir}/GRCh38.fa -t 16
```

```bash

usage: SERVE_merge.py [-h] [-i INPUT_GTF_LIST] [-p PREFIX] [--taco TACO]
                      [--stringtie STRINGTIE] [-n NSAMPLE] [-r REF_GENOME]
                      [-G GMAP_INDEX] [-g GMAP_INDEX_NAME] [-t NTHREAD]
                      [-l LENGTH] [--ratio RATIO] [-o OUTPUT_DIR]

SERVE_merge: merge expressed ERVs

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_GTF_LIST, --input_gtf_list INPUT_GTF_LIST
                        A text file with a list of SERVE GTF files (required)
  -p PREFIX, --prefix PREFIX
                        Prefix for output file name (default: SERVE)
  --taco TACO           Merge ERV transcripts with TACO (default: FALSE)
  --stringtie STRINGTIE
                        Merge ERV transcripts with StringTie (default: FALSE)
  -n NSAMPLE, --nsample NSAMPLE
                        The number of samples included in the input sample
                        list (required)
  -r REF_GENOME, --ref_genome REF_GENOME
                        Reference genome in FASTA format (required)
  -G GMAP_INDEX, --GMAP_index GMAP_INDEX
                        Path to the directory where GMAP index generated
                        (default: GMAP_index)
  -g GMAP_INDEX_NAME, --GMAP_index_name GMAP_INDEX_NAME
                        GMAP index name (default: GRCh38)
  -t NTHREAD, --nthread NTHREAD
                        Number of threads to run SERVE_merge (default: 1)
  -l LENGTH, --length LENGTH
                        Minimum ERV length (bp) (default: 200)
  --ratio RATIO         Minimum sample ratio of ERV identified (default: 0.50)
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory (default: .)
```

${prefix}_ERV_merge.gtf in the output directory records ERV gene information in specific condition. It can be used for downstream gene quantification (ERV_quant.py).

### SERVE_quant.py

SERVE_quant.py is designed for ERV quantification in single sample.

```bash

Example:
cat ${Dir}/GRCh38.gtf ${Dir}/${prefix}_ERV_merge.gtf >${Dir}/GRCh38_ERV.gtf
SERVE_quant.py -fq1 ${Dir}/test_1.fastq.gz -fq2 ${Dir}/test_2.fastq.gz -p test -r ${Dir}/GRCh38.fa -a ${Dir}/GRCh38_ERV.gtf -t 16
```

```bash

usage: SERVE_quant.py [-h] [-fq1 FASTQ1] [-fq2 FASTQ2] [-p PREFIX]
                      [-r REF_GENOME] [-a ANNOTATION] [-R RSEM_INDEX]
                      [-t NTHREAD] [-f FORWARD_PROB] [-o OUTPUT_DIR]

SERVE_quant: Quantify expressed ERVs

optional arguments:
  -h, --help            show this help message and exit
  -fq1 FASTQ1, --fastq1 FASTQ1
                        Read1 in FASTQ format (required)
  -fq2 FASTQ2, --fastq2 FASTQ2
                        Read1 in FASTQ format (required)
  -p PREFIX, --prefix PREFIX
                        Prefix for output file name (default: SERVE)
  -r REF_GENOME, --ref_genome REF_GENOME
                        Reference genome in FASTA format (required)
  -a ANNOTATION, --annotation ANNOTATION
                        Genome annotation in GTF format (required)
  -R RSEM_INDEX, --RSEM_index RSEM_INDEX
                        Path to the directory where RSEM index generated
                        (default: RSEM_index)
  -t NTHREAD, --nthread NTHREAD
                        Number of threads to run SERVE (default: 1)
  -f FORWARD_PROB, --forward_prob FORWARD_PROB
                        Probability of generating a read from the forward
                        strand of a transcript. 0.5 for unstranded-specific, 0
                        for stranded-specific where upstream reads are all
                        derived from the forward strand, 1 for stranded-
                        specific where upstream reads are all derived from the
                        reverse strand (default: 0.5)
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory (default: .)
```

${prefix}.genes.results in the output directory records gene quantification results. If you have biological replicates, you can do quality control analysis (ERV_quant_QC.py).

### SERVE_quant_QC.py

SERVE_quant_QC.py is designed for merge gene quantification results (produced by SERVE_quant.py) and quality control.

```bash

Example:
cat *genes.results >sample.list
SERVE_quant_QC.py -i ${Dir}/sample.list -p test
```

```bash

usage: SERVE_quant_QC.py [-h] [-i INPUT_SAMPLE_LIST] [-p PREFIX]
                         [--count COUNT] [--TPM TPM] [--ratio RATIO]
                         [-o OUTPUT_DIR]

SERVE_quant_QC: Quality control expressed ERVs

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_SAMPLE_LIST, --input_sample_list INPUT_SAMPLE_LIST
                        A text file with a list of SERVE quant files
                        (required)
  -p PREFIX, --prefix PREFIX
                        Prefix for output file name (default: SERVE)
  --count COUNT         Minimum ERV count (default: 5)
  --TPM TPM             Minimum ERV TPM (default: 0.1)
  --ratio RATIO         Minimum sample ratio of ERV identified (default: 0.50)
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory (default: .)
```

${prefix}_ERV_count.txt and ${prefix}_ERV_TPM.txt records gene read counts and TPM, respectively. Default: count > 5 & TPM >= 0.1 in over 50% samples.

### Example pipeline

Details in [example.sh](https://github.com/janky-yz/SERVE/blob/main/script/example.sh)

## Copyright and License Information

Copyright (C) 2021 Jianqi She (janky666@bjmu.edu.cn). See the [LICENSE](https://github.com/janky-yz/SERVE/blob/main/LICENSE) file for license rights and limitations.

