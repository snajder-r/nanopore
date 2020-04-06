#!/bin/bash


BASEDIR=/home/r933r/data/projects/nanopore/data/krebs_202ß

for batchdir in $BASEDIR/raw/single/*; do
    batchnum=$(basename $batchdir)
    outf=$BASEDIR/raw/multi/$batchnum
    fast5dir=$batchdir/workspace/0
    logfname=$BASEDIR/log/single_to_mult_$batchnum
    bsub -o $logfname.out -e $logfname.err -n 4 -R "span[hosts=1]" -M 4000 -W 00:30 "single_to_multi_fast5 -i $fast5dir -s $outf -n 4000 -t 4"
    exit
done
    
