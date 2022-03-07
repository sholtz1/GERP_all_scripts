#!/bin/bash

#SBATCH --account=beetlegenes
#SBATCH --time=120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sholtz1@uwyo.edu
#SBATCH --job-name=maf_processes



ref=./data/sequences/ref/Tribolium_ref_2020.fna
ref_name=$( basename $ref )
multiz_path=/project/beetlegenes/sholtz1/GERP/GERPsoftware/multiz-master
maf_array=($( ls -d ./analyses/last/net_axt/net_maf/*1to1.maf ))
combined_maf=./analyses/last/net_axt/net_maf/combined.maf
module load kentutils/1.04.0

mv Tribolium_ref_2020.fna_.gz.1to1.maf Tribolium_ref_2020.fna_aa.gz.1to1.maf
mv Tribolium_ref_2020.fna_.gz.1to1.maf_tmp Tribolium_ref_2020.fna_aa.gz.1to1.maf_tmp
mv Tribolium_ref_2020.fna_.gz.net.maf_w_header Tribolium_ref_2020.fna_aa.gz.net.maf_w_header


cat ${maf_array[@]:0:1} | sed -n '/##maf version=1 scoring=blastz/,$p' > \
${maf_array[@]:0:1}_tmp
cat ${maf_array[@]:1:1} | sed -n '/##maf version=1 scoring=blastz/,$p' > \
${maf_array[@]:1:1}_tmp

$multiz_path/multiz ${maf_array[@]:0:1}_tmp ${maf_array[@]:1:1}_tmp 1 > $combined_maf

for maf in ${maf_array[@]:2};
do
  echo "processing " $maf
 cat $maf | sed -n '/##maf version=1 scoring=blastz/,$p' > \
 "$maf"_tmp
$multiz_path/multiz $combined_maf "$maf"_tmp 1 > "$combined_maf"_tmp
 mv "$combined_maf"_tmp $combined_maf
done

# and filter mafs so all blocks have Zea mays and are at least 20 bp long
mafFilter -minCol=20 -needComp="$ref_name" $combined_maf > "$combined_maf".filtered
