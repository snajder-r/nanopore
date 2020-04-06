import os
import sys

python='/home/r933r/.conda/envs/gastrulation/bin/python'
reference='/home/r933r/data/projects/nanopore/data/krebs_2020/reference/dmel-all-chromosome-r6.29.fasta'
basedir='/home/r933r/data/projects/nanopore/data/krebs_2020/'

megalodon='/icgc/dkfzlsdf/analysis/B260/software/users/snajder/miniconda3/envs/megalodon/bin/megalodon'
guppy='/home/r933r/data/software/users/snajder/guppy/ont-guppy-cpu/bin/guppy_basecall_server'
megalodon_calibration='/icgc/dkfzlsdf/analysis/B260/software/users/snajder/miniconda3/envs/megalodon/lib/python3.7/site-packages/megalodon/model_data/dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac.cfg/megalodon_mod_calibration.npz'


per_batch=5

chroms = ['2L','2R', '3L', '3R', '4', 'X', 'Y', 'mitochondrion_genome']

#justsamples = ['invitro11_1', 'invitro12_1', 'invitro11_2', 'invitro12_2', 'footprinting']
justsamples = ['invitro11_1','invitro12_1']
samples = []
batches = []

samplebatches = dict()
def update_batches():
    for s in justsamples:
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

rule split_batches:
    output: directory(os.path.join(basedir,'raw/{sample}/batched/'))
    run: 
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
                src=os.path.join(basedir, os.path.join(sample_dir, raw_batches[i]))
                dst=os.path.join(batchdir, '%s' % raw_batches[i])
                os.symlink(src,dst)
                i+=1
            b+=1
        update_batches()

rule all_split_batches:
    input: expand(rules.split_batches.output, sample=justsamples)

rule fast5_to_fastq:
    input: os.path.join(basedir, 'raw/{sample}/batched/{batch}/')
    output: os.path.join(basedir, 'fastq/{sample}_{batch}.fq')
    params:
        jobname='fq_{sample}{batch}',
        runtime='03:00',
        memusage='8000',
        slots='1',
        misc=''
    shell: '{python} fast5_to_fastq.py {input} > {output}'


rule alignment:
    input: rules.fast5_to_fastq.output
    output: os.path.join(basedir, 'mapping/{sample}_{batch}.sorted.filtered.bam')
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

def prepare_mergebams_input(wildcards):
    return expand((os.path.join(basedir, 'mapping/%s' % wildcards.sample) + '_{batch}.sorted.filtered.bam'), batch=samplebatches[wildcards.sample])

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
        memusage='32000',
        slots='1',
        misc=''
    shell: '''
           module load samtools
           samtools merge -b {input} - | samtools calmd -b - {reference} | samtools sort -T {output}.tmp -o {output}
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
           sniffles -s 4 -t 16 -m {input} -v {output}
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

def merge_met_perchrom_input(wildcards):
    return expand(os.path.join(basedir, 'met/%s_{batch}_met.tsv' % wildcards.sample), batch=samplebatches[wildcards.sample])

rule merge_met_perchrom:
    input: merge_met_perchrom_input
    output: os.path.join(basedir, 'met_merged/{sample}_{chrom}_met.pkl')
    params:
        jobname='mergemet_{chrom}',
        runtime='4:00',
        memusage='16000',
        slots='1',
        misc=''
    shell: '{python} met_merge.py ' + os.path.join(basedir, 'met') + ' {wildcards.sample} {wildcards.chrom} {output}'


def allmetcall_combinator(*args, **kwargs):
    # All samples and batches
    for i in range(len(args[0])):
        # All mtypes
        for j in range(len(args[2])):
            #      sample   -  batch   -  mtype
            yield args[0][i],args[1][i],args[2][j]
        
def report_methylation_input(wildcards):
    return expand(rules.metcall.output, batch=samplebatches[wildcards.sample], sample=wildcards.sample, mtype=wildcards.mtype)

rule report_methylation:
    input: report_methylation_input
    output: os.path.join(basedir, 'report/{sample}_{mtype}.pdf')
    params:
        jobname='report_{sample}_{mtype}',
        runtime='0:10',
        memusage='4000',
        slots='1',
        misc=''
    shell: '{python} report_met.py ' + os.path.join(basedir,'met') + ' {output} {wildcards.sample} {wildcards.mtype}'

rule report_all_methylation:
    input: expand(rules.report_methylation.output, mtype=['cpg','gpc','dam'], sample=justsamples)


rule allmetcall:
    input: expand(os.path.join(basedir, 'met/{sample}_{batch}_met_{mtype}.tsv'), \
               allmetcall_combinator, sample=samples, batch=batches, mtype=['cpg','gpc','dam'])

def allmegalodon_combinator(*args, **kwargs):
    # All samples and batches
    for i in range(len(args[0])):
        #      sample   -  batch   -  mtype
        yield args[0][i],args[1][i]


rule allmegalodon:
    input: expand(os.path.join(basedir, 'megalodon/{sample}/{batch}/basecalls.modified_base_scores.hdf5'), \
               allmegalodon_combinator, sample=samples, batch=batches)

rule allmetmerged:
    input: expand(os.path.join(basedir, 'met_merged/{sample}_{chrom}_met.pkl'), sample=justsamples, chrom=chroms)


rule allvarcall:
    input: expand(os.path.join(basedir, 'var/{sample}_sv.vcf'), sample=justsamples)

rule test:
    input: os.path.join(basedir, 'megalodon/invitro11_1/0/basecalls.modified_base_scores.hdf5')

localrules: prepare_mergebams, split_batches
