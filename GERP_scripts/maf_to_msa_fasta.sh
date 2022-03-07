#!/bin/bash

#SBATCH --account=beetlegenes
#SBATCH --time=120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sholtz1@uwyo.edu
#SBATCH --job-name=mad_to_fasta





module load swset/2018.05
module load swset/2018.05
module load gcc/7.3.0
module load intel/18.0.1
module load perl/5.30.1
module load kentutils/1.04.0




ref_rm=./data/sequences/ref/Tribolium_N_masked.fna

ref_rm_name=$( basename $ref_rm )

combined_maf=./analyses/last/net_axt/net_maf/combined.maf

# splitting combined maf by target sequence
mkdir -p ./analyses/last/split_maf/
outroot=./analyses/last/split_maf/
mafSplit -byTarget dummy.bed ./analyses/last/split_maf/  ./analyses/last/net_axt/net_maf/combined.maf.filtered -useFullSequenceName

# mafSplit throws an error if you don't put dummy.bed, even though byTarget
# means it's being ignored

path_to_phast=/project/beetlegenes/sholtz1/GERP/GERPsoftware/phast
msa_view=${path_to_phast}/bin/msa_view
maf_dir=($( ls -d ./analyses/last/split_maf/*.maf ))


mkdir -p ./analyses/last/msa_fasta/
mkdir -p ./data/sequences/ref/split/

faSplit byname $ref_rm ./data/sequences/ref/split/

path_to_match_masking=./gerp/scripts-alignment/matchMasking.pl

for maf_file in "${maf_dir[@]}"
do
  chr=$(basename $maf_file | sed -e 's/0\(.*\).maf/\1/')
  ref_chr=$(basename $maf_file | sed -e 's/fna\.\(.*\).maf/\1/')

  fasta=./analyses/last/msa_fasta/"$chr".fa
  ref_rm_chr=./data/sequences/ref/split/"$ref_chr".fa
  rm_fasta=./analyses/last/msa_fasta/"$ref_chr"_rm.fa
  $msa_view $maf_file -i MAF -f -G 1 --refseq $ref_rm_chr > $fasta
  sed -E 's/> />/g' $fasta > "$fasta"_tmp && mv "$fasta"_tmp $fasta
  perl $path_to_match_masking --ref $ref_rm_chr --fasta $fasta --out $rm_fasta
done
