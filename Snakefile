import os
import sys
import h5py
import pandas as pd
import numpy as np

for k in config.keys():
    globals()[k] = config[k]

basecall_group = 'Basecall_1D_%s'%basecall_id
if not 'unique_samples' in globals().keys():
    unique_samples = os.listdir(os.path.join(basedir,'raw'))

samples = []
batches = []

'''
Detect batches for each sample. Fills the global variables:
  samples: flat array of length num_samples x num_batches
  batches: flat array of length num_samples x num_batches
  samplebatches: dict containing the list of batches per sample

samples and batches are filled such that zip(samples,batches) would result
in the full list of sample and batch tuples.
'''
samplebatches = dict()
def update_batches():
    for s in unique_samples:
        samplebatches[s] = []
        batchdir=os.path.join(basedir,'raw', s,'batched')
        if os.path.exists(batchdir):
            for b in os.listdir(batchdir):
                if os.path.isdir(os.path.join(batchdir,b)):
                    samples.append(s)
                    batches.append(b)
                    samplebatches[s].append(b)

update_batches()
print(samplebatches)


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
rule split_batches:
    output: directory(os.path.join(basedir,'raw','{sample}','batched','done'))
    params:
        outdir=os.path.join(basedir,'raw','{sample}','batched')
    run: 
        # Create "batched" directory if it doesn't exist
        if not os.path.exists(params.outdir):
            os.mkdir(params.outdir)
        i = 0
        b = 0
        sample_dir = os.path.join(basedir,'raw/%s/guppy/' % wildcards.sample)
        print(sample_dir)
        raw_batches = os.listdir(sample_dir)
        while i < len(raw_batches):
            batchdir=os.path.join(params.outdir, '%d'%b)
            os.mkdir(batchdir)
            for _ in range(per_batch):
                if i == len(raw_batches):
                     break
                src=os.path.join(basedir, os.path.join(sample_dir, 
                                 raw_batches[i]))
                dst=os.path.join(batchdir, '%s' % raw_batches[i])
                os.symlink(src,dst)
                i+=1
            b+=1
        update_batches()
        with open(output[0], 'w') as fp: 
            pass

checkpoint all_split_batches:
    input: expand(rules.split_batches.output, sample=unique_samples)


include: 'rules/fastq.rules'
include: 'rules/guppy.rules'
include: 'rules/mapping.rules'
include: 'rules/medaka.rules'
include: 'rules/megalodon.rules'
include: 'rules/nanopolish.rules'
include: 'rules/pycoqc.rules'
include: 'rules/sniffles.rules'
include: 'rules/tombo.rules'
include: 'rules/edgecase.rules'

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
    shell: '{python} report_met.py ' + os.path.join(basedir,'met') + \
           ' {output} {wildcards.sample} {wildcards.mtype}'

rule all_report_methylation:
    input: expand(rules.report_methylation.output, mtype=mettypes, sample=unique_samples)

localrules: prepare_mergebams, split_batches
