#!/bin/bash


cluster_script="bsub {params.misc} -n {params.slots} -W {params.runtime} -R \"rusage[mem={params.memusage}]\" -o ../logs/{params.jobname}.log -e ../logs/{params.jobname}.err"

snakemake --cluster "$cluster_script" --jobs 64 --latency-wait 120 $@

