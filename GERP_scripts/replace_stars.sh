#!/bin/bash

cd /project/beetlegenes/sholtz1/GERP/analyses/last/msa_fasta


maf_dir=($(ls  *.fa ))


for fa_file in "${maf_dir[@]}"
do
chr=$(basename $fa_file)

awk '{if(/^[^>]/)gsub(/\*/,"-");print $0}' $fa_file  > "$chr".test

sed -i "1s/.*/Tribolium_castaneum/" $fa_file

done
