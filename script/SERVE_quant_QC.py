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


parser = argparse.ArgumentParser(description='SERVE_quant_QC: Quality control expressed ERVs')
parser.add_argument('-i', '--input_sample_list', help='A text file with a list of SERVE quant files (required)')
parser.add_argument('-p', '--prefix', default='SERVE', help='Prefix for output file name (default: SERVE)')
parser.add_argument('--count', default=5, help='Minimum ERV count (default: 5)')
parser.add_argument('--TPM', default=0.1, help='Minimum ERV TPM (default: 0.1)')
parser.add_argument('--ratio', default=0.50, help='Minimum sample ratio of ERV identified (default: 0.50)')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory (default: .)')

args = parser.parse_args()
script_dir = os.path.abspath(os.path.dirname(__file__))


if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)


print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Running SERVE_quant_QC.', flush=True)


with cd(args.output_dir):

	cmd = 'Rscript '+script_dir+'/SERVE_quant_QC.R '+args.input_sample_list+' '+args.prefix+' '+str(args.count)+' '+str(args.TPM)+' '+str(args.ratio)
	subprocess.check_call(cmd, shell=True, executable='/bin/bash')


print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] expressed ERVs quality control is done.', flush=True)
