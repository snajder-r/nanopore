import os
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

for k in config.keys():
    globals()[k] = config[k]

if "min_read_q" not in globals():
    min_read_q = 0

basecall_group = 'Basecall_1D_%s' % basecall_id
if 'unique_samples' not in config.keys():
    unique_samples = os.listdir(os.path.join(basedir, 'raw'))
if 'fastq_ending' not in config.keys():
    fastq_ending = 'fq'
if 'chroms' not in config.keys():
    with open(reference, "r") as f:
        print("Trying to infer chromosomes from reference")
        chroms = [l[1:].split(" ")[0] for l in f if l.startswith(">")]

chroms = [str(c) for c in chroms]

#basedir = Path(basedir)

'''
Detect batches for each sample. Fills the global variables:
  samples: flat array of length num_samples x num_batches
  batches: flat array of length num_samples x num_batches
  samplebatches: dict containing the list of batches per sample

samples and batches are filled such that zip(samples,batches) would result
in the full list of sample and batch tuples.
'''

class SampleBatchesFilenames:
    def __init__(self):
        if os.path.exists(os.path.join(basedir, 'raw')):
            sbf = glob_wildcards(os.path.join(basedir, 'raw', '{sample}', 'batched', '{batch}', '{filename}.fast5'))
            self.sbf_samples = sbf.sample
            self.sbf_batches = sbf.batch
            self.sbf_filenames = sbf.filename

            sbdict = {}
            for s, b in zip(sbf.sample, sbf.batch):
                sbdict[f"{s}_{b}"] = (s,b)
            self.sb_samples = [s for s,_ in sbdict.values()]
            self.sb_batches = [b for _,b in sbdict.values()]
            
        elif os.path.exists(os.path.join(basedir, 'fastq')):
            sb = glob_wildcards(os.path.join(basedir, 'fastq', '{sample}', '{batch}.%s' % fastq_ending))
            self.sb_samples = sb.sample
            self.sb_batches = sb.batch
            self.sbf_samples = None
            self.sbf_batches = None
            self.sbf_filenames = None

sbf = SampleBatchesFilenames()

def samplebatches(sample):
    return [sbf.sb_batches[i] for i in range(len(sbf.sb_batches)) if sbf.sb_samples[i] == sample]


'''
##############################################################################
# Combinators 
##############################################################################

Below we define some combinators of wildcards. These can be used in "all_"
type rules. The expand function in Snakemake assumes that the list of 
wildcard combinations is always the full cross-product. In our case, this
is not the case, as samples have different number of batches. 

Therefore we define a few custom combinators below
'''

''' 
This combinator assumes that all wildcard variables are simply zipped 
together. This is designed for the sample and batch variables.
'''


def zip_combinator(*args, **kwargs):
    # All samples and batches
    for i in range(len(args[0])):
        yield tuple(args[j][i] for j in range(len(args)))


'''
This combinator zips the first two arguments, and combines those with the 
third argument. We use this to create all sample-batch-mtype combinations.
'''


def zip2_comb3_combinator(*args, **kwargs):
    # First two wildcards
    for i in range(len(args[0])):
        # Combine with third wildcard
        for j in range(len(args[2])):
            #     wc1   -     wc2   -     wc3
            yield args[0][i], args[1][i], args[2][j]


'''
##############################################################################
# Splitting fast5 files into batches for parallel processing computation
##############################################################################

The input directory raw/samplename/guppy/ will be searched for fast5 files.
The list of fast5 files per sample is split into evenly sized batches.

The directory raw/samplename/batched will contain a subdirectory for each 
batch.

CAVEAT: The fast5 files will be **symlinked** from the batch directory to the
basecalled directory. Do NOT delete the fast5 files from the basecalled
directory!

Rule will not be performed if "batched" directory exists. Delete the "batched"
directory, if you want to redo batch splitting.

This rule will be performed locally, as it is only creating symlinks and is
not computationally expensive.
'''


def split_batches_from_file_list(all_files, outdir):
    # Create "batched" directory if it doesn't exist
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    i = 0
    b = 0
    while i < len(all_files):
        batchdir = os.path.join(outdir, '%d' % b)
        os.mkdir(batchdir)
        for _ in range(per_batch):
            if i == len(all_files):
                break
            link = os.path.join(batchdir, '%s' % os.path.basename(all_files[i]))
            os.symlink(all_files[i], link)
            i += 1
        b += 1


def split_batches_input(wildcards):
    #multi_file_path = os.path.join(basedir, 'raw', wildcards['sample'], 'multi', '{fname}.fast5')
    file_path = os.path.join(basedir, 'raw', wildcards['sample'], 'guppy', '{fname}.fast5')
    return expand(file_path, fname=glob_wildcards(file_path).fname)


rule split_batches:
    input: split_batches_input
    output: os.path.join(basedir, 'raw', '{sample}', 'batched', 'done')
    params:
          outdir=os.path.join(basedir, 'raw', '{sample}', 'batched')
    run:
        split_batches_from_file_list(input, params.outdir)
        with open(output[0], 'w') as fp:
            pass


def split_batches_from_fastq_input(wildcards):
    file_path = os.path.join(basedir, 'fastq', wildcards['sample'], 'guppy', '{fname}.%s' % fastq_ending)
    return expand(file_path, fname=glob_wildcards(file_path).fname)


rule split_batches_from_fastq:
    input: split_batches_from_fastq_input
    output: os.path.join(basedir, 'fastq', '{sample}', 'batched', 'done')
    params:
          outdir=os.path.join(basedir, 'fastq', '{sample}', 'batched')
    run:
        split_batches_from_file_list(input, params.outdir)
        with open(output[0], 'w') as fp:
            pass

rule all_split_batches:
    input: expand(rules.split_batches.output, sample=unique_samples)

rule all_split_batches_from_fastq:
    input: expand(rules.split_batches_from_fastq.output, sample=unique_samples)

include: 'rules/fastq.rules'
include: 'rules/guppy.rules'
include: 'rules/mapping.rules'
include: 'rules/nanopolish.rules'
include: 'rules/pycoqc.rules'

if "nanonome" in globals():
    include: "rules/nanonome.rules"

if "bonito" in globals():
    include: "rules/bonito.rules"

if "accessibility_model" in globals():
    include: 'rules/accessibility_prediction.rules'
else:
    print("accessibility_model not defined in config. Support for accessibility prediction disabled")

if "whatshap" in globals():
    include: "rules/whatshap.rules"

if "pycometh" in globals():
    include: "rules/pycometh.rules"

if "custom_rules" in globals().keys():
    print(custom_rules)
    for f in Path("rules.d").iterdir():
        if f.name.endswith(".rules"):
            if f.name in custom_rules:
                include: f



'''
##############################################################################
# Reports methylation histograms - only really useful for benchmarking
# datasets where you know the methylation rate (like fully methylated or
# fully unmethylated datasets)
##############################################################################
'''


def report_methylation_input(wildcards):
    return expand(rules.metcall.output, batch=samplebatches[wildcards.sample],
                  sample=wildcards.sample, mtype=wildcards.mtype)


rule report_methylation:
    input: report_methylation_input
    output: os.path.join(basedir, 'report/{sample}_{mtype}.pdf')
    params:
          jobname='report_{sample}_{mtype}',
          runtime='0:10',
          memusage='4000',
          slots='1',
          misc=''
    shell: '{python} report_met.py ' + os.path.join(basedir, 'met') + \
         ' {output} {wildcards.sample} {wildcards.mtype}'

rule all_report_methylation:
    input: expand(rules.report_methylation.output, mtype=mettypes, sample=unique_samples)

localrules: prepare_mergebams, split_batches, split_batches_from_fastq
