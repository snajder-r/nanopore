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
# GUPPY basecalling
##############################################################################
'''

rule guppy_basecall:
    input: os.path.join(basedir, 'raw', '{sample}', 'multi')
    output: os.path.join(basedir, 'raw', '{sample}', 'guppy')
    params:
        jobname='guppy_{sample}',
        runtime='48:00',
        memusage='16000',
        slots='1',
        misc = '-q gputest -gpu num=1:j_exclusive=yes:mode=exclusive_process:gmem=10G'
    shell: '{python} basecall_guppy_gpucluster.py {input} {output}'

rule all_guppy_basecall:
    input: expand(rules.guppy_basecall.output, sample=unique_samples)

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
    output: directory(os.path.join(basedir,'raw','{sample}','batched'))
    run: 
        # Create "batched" directory if it doesn't exist
        if not os.path.exists(output[0]):
            os.mkdir(output[0])
        i = 0
        b = 0
        sample_dir = os.path.join(basedir,'raw/%s/guppy/' % wildcards.sample)
        print(sample_dir)
        raw_batches = os.listdir(sample_dir)
        while i < len(raw_batches):
            batchdir=os.path.join(output[0], '%d'%b)
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

rule all_split_batches:
    input: expand(rules.split_batches.output, sample=unique_samples)


'''
##############################################################################
# Extracting FASTQ from basecalled FAST5 files
##############################################################################

This rule will open each basecalled fast5 for a sample and batch, and create
a single fastq file as an outputfile, containing the reads for this sample 
and batch.

The rule depends on the "basecall_id" variable in the config file. This would
typically be "Basecall_1D_000" if there was only a single basecalling 
performed. If there were multiple basecalls (e.g. basecalled with different 
versions of guppy), make sure you specify the correct basecall id in the 
config file.

If individual fast5 files could not be opened or basecalls could not be found,
the entire file will be skipped and a warning written to the log. This is
done because guppy will sometimes file to produce basecalls, and is such a
case we just skip that file and move on.
'''

rule fast5_to_fastq:
    input: rules.split_batches.output
    output: os.path.join(basedir, 'fastq', '{sample}_{batch}.fq')
    params:
        jobname='fq_{sample}{batch}',
        runtime='03:00',
        memusage='8000',
        slots='1',
        misc=''
    run:
        # Recursively finds fast5 files
        def find_fast5(indir, fl=[]):
            ll = os.listdir(indir)
            for l in ll:
                l = os.path.join(indir,l)
                if os.path.isdir(l):
                    fl = fl + find_fast5(l, fl)
                elif l.endswith('.fast5'):
                    fl.append(l)
            return fl

        indir = input[0]
        fl = find_fast5(indir)

        fastq_group='Analyses/{basecall_group}/BaseCalled_template/Fastq'.format(
                        basecall_group=basecall_group)
        with open(output[0],'w') as of:
            for fn in fl:
#                try:
                    with h5py.File(fn,'r') as f:
                        for read in f.keys():
                            print(f[read].keys())
                            fq = f[read][fastq_group][()].decode('utf8')
                            fq = fq.split('\n')
                            if read.startswith('read_'):
                                read = read[5:]
                            fq[0] = '@%s'%read

                            of.write('\n'.join(fq))
#                except:
#                    print('WARN: Could not read basecalls in %s'%fn)
#                    pass

rule all_fast5_to_fastq:
    input: expand(rules.fast5_to_fastq.output, zip_combinator, sample=samples,  batch=batches)

'''
##############################################################################
# Mapping using minimap2
##############################################################################

This rule will map fastq file for a single sample and batch to the reference. 

Output will be stored in the "mapping" directory as a bam file including a 
bai index file.

The following steps are performed:
  - Alignment
  - Filter only unique mappings
  - Sort bam file
  - Create bam index
'''
rule alignment:
    input: rules.fast5_to_fastq.output
    output: os.path.join(basedir,'mapping/{sample}_{batch}.sorted.filtered.bam')
    params:
        jobname='align_{sample}_{batch}',
        runtime='09:00',
        memusage='2000',
        slots='1',
        misc=''
    shell: '''
           module load samtools
           minimap2 -a -x map-ont {reference} {input} | samtools view -q 1 -b | samtools sort -T {output}.tmp -o {output}
           samtools index {output}
           '''

rule all_alignment:
    input: expand(rules.alignment.output, zip_combinator, sample=samples, batch=batches)

'''
##############################################################################
# Methylation calling using nanopolish
##############################################################################
'''

'''
In preparation for nanopolish methylation calling, raw fast5 entries 
for each read in the fastq files are indexed so they can be accessed more 
efficiently by nanopolish.

The output files will be created within the fastq directory 
(with .index.readdb suffix for each fq file) 
'''
rule nanopolish_index:
    input:
        f5=os.path.join(basedir, 'raw/{sample}/batched/{batch}/'),
        fq=rules.fast5_to_fastq.output
    output: os.path.join(basedir, 'fastq/{sample}_{batch}.fq.index.readdb')
    params:
        jobname='npidx_{sample}_{batch}',
        runtime='08:00',
        memusage='2000',
        slots='1',
        misc=''
    shell: 'nanopolish index -d {input.f5} {input.fq}'

rule all_nanopolish_index:
    input: expand(rules.nanopolish_index.output, zip_combinator, sample=samples, batch=batches)


'''
Performs methylation calling using nanopolish for one sample and batch, for
one type of methylation call (the "mtype" wildcard). The "mtype" wildcard 
could be one of "cpg", "gpc", or "dam".

The output will be stored in the "met" directory in tsv format.
'''
rule metcall:
    input: 
        fq=rules.fast5_to_fastq.output,
        bam=rules.alignment.output,
        fqidx=rules.nanopolish_index.output
    output: os.path.join(basedir, 'met/{sample}_{batch}_met_{mtype}.tsv')
    params:
        jobname='metcall_{sample}_{batch}',
        runtime='08:00',
        memusage='8000',
        slots='8 -R "span[hosts=1]"',
        misc=''
    shell: 'nanopolish call-methylation -t 8 -g {reference} -b {input.bam} -r {input.fq} -q {wildcards.mtype} > {output}'

rule all_metcall:
    input: expand(rules.metcall.output, zip2_comb3_combinator, sample=samples, batch=batches, mtype=mettypes)

'''
This merges the nanopolish methylation calls, loads it into a pandas dataframe
and stores it in pickled format.

Note that this job is performed per sample per chromosome. It will save one
outputfile for a sample,chromosome,methylation type combination, in order to 
break up the data and make it easier to load loater.

So while it merges batches, it also splits up chromosomes.
'''
def merge_met_perchrom_input(wildcards):
    return expand(os.path.join(basedir, 
        'met/%s_{batch}_met_%s.tsv' % (wildcards.sample, wildcards.mettype)),
        batch=samplebatches[wildcards.sample])

rule merge_met_perchrom:
    input: merge_met_perchrom_input
    output: os.path.join(basedir, 'met_merged/{sample}_chr{chrom}_met_{mettype}.pkl')
    params:
        jobname='mergemet_{sample}_{chrom}',
        runtime='12:00',
        memusage='24000',
        slots='1',
        misc=''
    run:
        metcall_dir = os.path.join(basedir, 'met') 
        sample = wildcards['sample']
        chrom = wildcards['chrom']
        mettype = wildcards['mettype']
        pickle_file = output[0]

        sample_met = None
        sample_batch_re = re.compile('%s_(.*)_met_%s\.tsv' % (sample, mettype))
        for f in os.listdir(metcall_dir):
            mitch = sample_batch_re.match(f)

            if mitch is None:
                continue

            batch = mitch.group(1)
            
            for lmet in pd.read_csv(os.path.join(metcall_dir,f), sep='\t', header=0, chunksize=1000000,
                               dtype={'chromosome':'category', 
                                      'strand':'category', 
                                      'num_calling_strands':np.uint8, 
                                      'num_motifs':np.uint8}):

                lmet = lmet.loc[lmet.chromosome==chrom].copy()
                if sample_met is None:
                    sample_met = lmet
                else:
                    sample_met = pd.concat((sample_met, lmet))

        sample_met.to_pickle(pickle_file, compression='gzip')

rule all_merge_met_perchrom:
    input: expand(rules.merge_met_perchrom.output, sample=unique_samples, chrom=chroms, mettype=mettypes)


'''
##############################################################################
# Merging of BAM files
##############################################################################

In case you need a single bam file for one sample (e.g. for SV calling), this 
rule performs the merging and re-sorting of bam files. The 
prepare_mergebams_input function returns the list of bam files, which will 
then in an intermediate step be saved to a txt. 

The following steps are performed:
 - Create filelist file
 - Merge using samtools merge
 - Use calmd to re-generate the MD tag for the merged bam file
 - Re-sort using samtools sort  (not sure if this is necessary after merge)
'''
def prepare_mergebams_input(wildcards):
    return expand((os.path.join(basedir, 'mapping/%s' % wildcards.sample) + 
        '_{batch}.sorted.filtered.bam'), 
        batch=samplebatches[wildcards.sample])

rule prepare_mergebams:
    input: prepare_mergebams_input
    output: os.path.join(basedir, 'mapping/{sample}.bam.filelist.txt')
    shell: 'echo {input} | sed \'s/ /\\n/g\' > {output}'

rule mergebams:
    input: rules.prepare_mergebams.output
    output: os.path.join(basedir, 'mapping/{sample}.sorted.bam')
    params:
        jobname='mergebam_{sample}',
        runtime='09:59',
        memusage='64000',
        slots='1',
        misc=''
    shell: '''
           module load samtools
           samtools merge -b {input} - | samtools calmd -b - {reference} | samtools sort -T {output}.tmp -o {output}
           '''

rule all_mergebams:
    input: expand(rules.mergebams.output, sample=unique_samples)

'''
##############################################################################
# Calling of SVs using sniffles
##############################################################################

Uses merged bams as input. Since sniffles considers the amount of evidence 
for an SV, it is important that it has access to all reads, hence it needs a
merged bam file. Therefore it depends on the mergebams rule and can run only
per sample and not batchwise.
'''
rule varcall:
    input:
        bam=rules.mergebams.output
    output: os.path.join(basedir, 'var/{sample}_sv.vcf')
    params:
        jobname='varcall_{sample}',
        runtime='24:00',
        memusage='16000',
        slots='16 -R "span[hosts=1]"',
        misc=''
    shell: '''
           module load gcc/5.4.0
           {sniffles} -s 4 -t 16 -m {input} -v {output}
           '''

rule all_varcall:
    input: expand(rules.varcall.output, sample=unique_samples)


'''
##############################################################################
# Methylation calling using Megalodon
##############################################################################

Use Megalodon without GPU support to perform methylation calling. This 
requires guppy methylation called fast5 files.

The config file must contain the megalodon_calibration, which points to the
calibration filename specific to the machine and modification type.

Also requires the path to the guppy basecalling server in the "guppy" 
variable.

While Megalodon supports GPU computation, I found it is typically faster to do
batched CPU computations, since I can do more parallelization than using the
limited number of GPUs available on the DKFZ cluster.
'''
rule megalodon:
    input: os.path.join(basedir, 'raw/{sample}/batched/{batch}/')
    output: os.path.join(basedir, 'megalodon/{sample}/{batch}/basecalls.modified_base_scores.hdf5')
    params: 
        jobname='megalodon_{sample}_{batch}',
        runtime='48:00',
        memusage='16000',
        slots='16 -R "span[hosts=1]"',
        misc=''
    shell: '{megalodon} --guppy-params "--cpu_threads_per_caller 16" --guppy-server-path {guppy} --overwrite --output-directory %s --mod-calibration-filename {megalodon_calibration} --outputs=mod_basecalls {input}' % os.path.join(basedir, 'megalodon/{wildcards.sample}/{wildcards.batch}/')

rule all_megalodon:
    input: expand(rules.megalodon.output, zip_combinator, sample=samples, batch=batches)

'''
#############################################################################
# MEDAKA methylation calling
#############################################################################

This pipeline performs its own alignment (albeit also just using minimap2),
which writes the methylation probabilities into the bam files
'''

'''
Alignment and extraction of methylation probabilities
'''
rule medaka_align:
    input: os.path.join(basedir, 'raw/{sample}/batched/{batch}/')
    output: bam=os.path.join(basedir, 'medaka/{sample}_{batch}.bam'),
            index=os.path.join(basedir, 'medaka/{sample}_{batch}.bam.bai')
    params:
        jobname='medaka_align_{sample}_{batch}',
        runtime='12:00',
        memusage='16000',
        slots='32 -R "span[hosts=1]"',
        misc=''
    shell: '''
           module load samtools;
           {medaka} methylation guppy2sam {input} --reference {reference} --workers 31 --recursive | samtools sort -@ 8 | samtools view -b -@ 8 > {output.bam};
           samtools index {output.bam};
           '''

rule all_medaka_align:
    input: expand(rules.medaka_align.output, zip_combinator, sample=samples, batch=batches)

rule medaka_metcall:
    input: rules.medaka_align.output.bam
    output: os.path.join(basedir, 'medaka/{sample}_{batch}_met_{chrom}.tsv')
    params:
        jobname='medaka_metcall_{sample}_{batch}_{chrom}',
        runtime='1:00',
        memusage='16000',
        slots='1',
        misc=''
    shell: '{medaka} methylation call --meth all {input} {reference} {wildcards.chrom}:0-10000000 {output}' 

rule testmedaka:
    input: expand(rules.medaka_metcall.output, zip2_comb3_combinator, sample=['all'], batch=['0'], chrom=[chroms[0]])

rule all_medaka_metcall:
    input: expand(rules.medaka_metcall.output, zip2_comb3_combinator, sample=samples, batch=batches, chrom=chroms)

'''
##############################################################################
# TOMBO Methylation calling
##############################################################################

Note that tombo requires SINGLE fast5 files. They are expected to be in the
directory raw/samplename/single.
'''

'''
Creates index of fast5 raw files for Tombo
'''
rule tombo_make_index:
    input: os.path.join(basedir, 'raw/{sample}/single/')
    output: os.path.join(basedir, 'raw/{sample}/.single.RawGenomeCorrected_000.tombo.index')
    params:
        jobname='tombo_make_index',
        runtime='24:00',
        memusage='16000',
        slots='1',
        misc=''
    shell: '{tombo} filter clear_filters --fast5-basedirs {input}'

'''
Uses Tombo resquiggle to segment the raw signal and align it to the bases.
Note that the output file ".single.resquiggled" does not really contain any
information. It's only there as a marker that resquiggling has been performed.
Tombo resquiggle modifies the input files, so the actual output is actually
written into the inputfiles.

CAVEAT: This manipulates the original fast5 files.
'''
rule tombo_resquiggle:
    input: 
        index=rules.tombo_make_index.output,
        directory=rules.tombo_make_index.input
    output: os.path.join(basedir, 'raw/{sample}/.single.resquiggled')
    params:
        jobname='tombo_resquiggle',
        runtime='48:00',
        memusage='16000',
        slots='32 -R "span[hosts=1]"',
        misc=''
    shell: '''
           {tombo} resquiggle --processes 32 --dna {input.directory} {reference};
           RESULT=$?
           if [[ $RESULT -eq 0 ]] ; then touch {output}; fi
           exit $RESULT
           '''

'''
Perform tombo methylation calling for 6mA and 5mC
'''
rule tombo_metcall:
    input: resquiggled=rules.tombo_resquiggle.output,
           directory=rules.tombo_make_index.input
    output: stats_6ma=os.path.join(basedir, 'tombo/{sample}.6mA.tombo.stats'),
            stats_5mc=os.path.join(basedir, 'tombo/{sample}.5mC.tombo.stats')
    params:
        jobname='tombo_metcall',
        runtime='48:00',
        memusage='16000',
        slots='32 -R "span[hosts=1]"',
        misc=''
    shell: '''
           {tombo} detect_modifications alternative_model --fast5-basedirs {input.directory} --statistics-file-basename {basedir}/tombo/{wildcards.sample} --alternate-bases 6mA 5mC --processes 32
           '''  
        
rule all_tombo_metcall:
    input: expand(rules.tombo_metcall.output.stats_5mc, sample=unique_samples)

'''
Perform tombo methylation calling for 6mA and 5mC
'''
rule tombo_de_novo:
    input: resquiggled=rules.tombo_resquiggle.output,
           directory=rules.tombo_make_index.input
    output: stats=os.path.join(basedir, 'tombo/{sample}_denovo.tombo.stats'),
    params:
        jobname='tombo_metcall',
        runtime='48:00',
        memusage='16000',
        slots='32 -R "span[hosts=1]"',
        misc=''
    shell: '''
           {tombo} detect_modifications de_novo --fast5-basedirs {input.directory} --statistics-file-basename {basedir}/tombo/{wildcards.sample}_denovo --processes 32
           '''  

rule all_tombo_de_novo:
    input: expand(rules.tombo_de_novo.output.stats, sample=unique_samples)

'''
#############################################################################
# Use Adrien Ledger's pycoQC to generate quality reports
#############################################################################
'''

rule fast5_to_seq_summary:
    input: os.path.join(basedir,'raw/{sample}/guppy/')
    output: os.path.join(basedir, 'raw/{sample}/guppy_summary.txt')
    params:
        jobname='fast5_to_seq_summary_{sample}',
        runtime='24:00',
        memusage='16000',
        slots='32 -R "span[hosts=1]"',
        misc=''
    shell: '{Fast5_to_seq_summary} --fast5_dir {input} --seq_summary_fn {output} --threads 32 --basecall_id {basecall_id}'

rule all_fast5_to_seq_summary:
    input: expand(rules.fast5_to_seq_summary.output, sample=unique_samples)

rule pycoqc_report:
    input: 
        summary = rules.fast5_to_seq_summary.output,
        bams = lambda wildcards: expand(rules.alignment.output, batch=samplebatches[wildcards.sample], sample=wildcards.sample)
    output:
        os.path.join(basedir, 'report', '{sample}_pycoqc.html')
    params:
        jobname='pycoqc_{sample}',
        runtime='2:00',
        memusage='16000',
        slots='1',
        misc=''
    shell: '{pycoQC} -f {input.summary} --bam_file {input.bams} --html_outfile {output}'

rule all_pycoqc_report:
    input: expand(rules.pycoqc_report.output, sample=unique_samples)

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
