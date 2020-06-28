#!/bin/bash
reference=$1
read1=$2
read2=$3
samp=$4

bwa index $reference

bwa mem $reference $read1 $read2 > bwa_$samp.sam
