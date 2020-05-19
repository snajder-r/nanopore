# Nanopore methylation calling pipeline

Snakemake pipeline performing the following steps:

## Nanopolish pipeline

This is the main pipeline

* Splitting of basecalled fast5 files into batches
* Fast5 to fastq conversion
* Alignment to reference using minimap2
* Methylation calling for cpg, gpc and dam using Nanopolish
* Convert methylation call files into pickled pandas dataframe split per chromosome

## Other pipelines:

These are also implemented, albeit more or less experimental:

### Tombo pipeline

Note that this requires single-fast5 files. Use the script '''multi_to_single_fast5''' from the ont_fast5_api in order to convert.

* Tombo indexing
* Tombo resquiggle
* Tombo methylation calling

### Megalodon:

This requires guppy basecall server installed.

* Alignment-free methylation calling using megalodon

Note that in this pipeline, Megalodon is configured perform batchwise and in parallel using CPUs.
GPU accellerated computation is supported, but in a setup where you have only few GPU resources
and many CPU resources, the high paralellization of a CPU cluster is actually faster.

## Data structure:

Please provide input data in the following structure. You only need to provide the *guppy* folder for the Nanopolish pipeline.
The *multi* folder is only required if you want to use the script for basecalling (not part of the main pipeline).

~~~~
basedir/
 |-raw/
    |-<sample1>/
      |-guppy/ (basecalled multi-fast5 files)
        |-*.fast5       
      |-multi/ (non-basecalled multi-fast5 files)
        |-*.fast5 
      |-single/ (basecalled single-fast5 files for tombo)
        |-0/
          |-*.fast5 
        |-1/
          |-*.fast5
        |-[...]
    |-<sample2>/
      |-[...]
    |-[...]
~~~~

The pipeline will create the following data structure, with the variables <sample> referring to the samplename, <batchnum> the batch number within that sample, <chr> for chromosome, and <mettype> the type of methylation call (cpg,gpc, or dam):

~~~~
basedir/
  |-raw/
    |-<sample>/
      |-batched/
        |-<batchnum>/
          |-*.fast5 (symlink to basedir/<sample>/guppy/*.fast5)
  |-fastq/ (fastq files and indices for samtools and nanopolish)
    |-<sample>_<batchnum>.fastq
    |-<sample>_<batchnum>.fastq.index
    |-<sample>_<batchnum>.fastq.fai
    |-<sample>_<batchnum>.fastq.gzi
    |-<sample>_<batchnum>.fastq.read_db
  |-mapping/ (bam files sorted and filtered)
    |-<sample>_<batchnum>.sorted.filtered.bam
    |-<sample>_<batchnum>.sorted.filtered.bai
  |-met/ (nanopolish methylation calls)
    |-<sample>_<batchnum>_met_<mettype>.tsv
  |-met_merge/ (batches merged back together, split by chromosome, and stored as pickled pandas frame)
    |-<sample>_<chr>_met_<mettype>.pkl
  |-logs/ 
    |-<jobname>.log
    |-<jobname>.err
~~~~

## Configuration

The pipeline is configured using a YAML file. An example is provided below.
Not all paths to programs are required if they are not part of the main pipeline.

~~~~
#### Minimum required parameters ####
# Base directory of the project
basedir: /home/r933r/snajder/nanopore/data/medulloblastoma_dna/
# Path to python that has the required packages installed
python: /home/r933r/.conda/envs/gastrulation/bin/python
# Reference (must be uncompressed and have an index generated using samtools faidx)
reference: /home/r933r/snajder/nanopore/data/medulloblastoma_dna/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa
# Number of multi-fast5 files per batch
per_batch: 15
# Names of chromosomes to be used (all chromosomes not in the list will be filtered out)
chroms: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, X, Y, MT]
# Names of samples
unique_samples: [s1,s2,s3]
# Name of the group that contains the basecalls in the fast5 file
basecall_id: Basecall_1D_000
# Which methylation types to call with nanopolish. Supported: cpg, gpc, dam
mettypes: [cpg,gpc,dam]

#### Additional configuration (optional) ####
fast5_single_to_multi: /home/r933r/.conda/envs/gastrulation/bin/multi_to_single_fast5
tombo: /home/r933r/.conda/envs/gastrulation/bin/tombo
guppy_bc_server: /home/r933r/data/software/users/snajder/guppy/ont-guppy-cpu/bin/guppy_basecall_server
sniffles: /icgc/dkfzlsdf/analysis/B260/software/users/snajder/opt/bin/sniffles
megalodon: /icgc/dkfzlsdf/analysis/B260/software/users/snajder/miniconda3/envs/megalodon/bin/megalodon
# Megalodon calibration file
megalodon_calibration: /icgc/dkfzlsdf/analysis/B260/software/users/snajder/miniconda3/envs/megalodon/lib/python3.7/site-packages/megalodon/model_data/dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac.cfg/megalodon_mod_calibration.npz
~~~~

## How to run

The pipeline is set up to run on the DKFZ cluster, but should work on any cluster system.
For ease of use, I typically use a script like "do_snakemake.sh" to run snakemake with all the parameters required to run on the cluster and save logs.

If I use the following script:

do_snakemake.sh:

~~~~
#!/bin/bash

cluster_script="bsub {params.misc} -n {params.slots} -W {params.runtime} -R \"rusage[mem={params.memusage}]\" -o {basedir}/logs/{params.jobname}.log -e {basedir}/logs/{params.jobname}.err"

snakemake --cluster "$cluster_script" --jobs 128 --latency-wait 120 $@
~~~

Then I run it using:

~~~~
do_snakemake.sh --configfile <myprojectconfig.yaml> <target>
~~~~

### Targets:

Assuming you have your basecalled fast5 files in the raw/<sample>/guppy folders, you first need to do splitting into batches. To do this, first run the target

* all_split_batches

Afterwards, you can run the full Nanopolish pipeline by running the target:

* all_merge_met_perchrom

But you can also run any of the following intermediate steps:

* all_fast5_to_fastq: to get fastq files from basecalled fast5 files
* all_alignment: runs alignment
* all_nanopolish_index: Creates index required by nanopolish
* all_metcall: Performs nanopolish methylation calling for all methylation types

Some optional targets (partly experimental - may be unstable):

* all_mergebams: Merges all bam files
* all_varcall: Performs SV calling using sniffles
* all_megalodon: Performs methylation calling using megalodon
* all_tombo_metcall: Performs methylation calling using tombo
* all_report_methylation: Creates some histograms for methylation frequency (very global view on methylation, only useful if you know the expected methylation rate)

