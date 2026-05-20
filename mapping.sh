#bin/bash

input_file=$1
ref_fa=$2
output_prefix=$3

mkdir -p mapping
minimap2 -a -x asm5 $ref_fa $input_file -o mapping/$output_prefix.sam
cd mapping
samtools sort $output_prefix.sam -O BAM --write-index -o $output_prefix.bam
rm $output_prefix.sam

