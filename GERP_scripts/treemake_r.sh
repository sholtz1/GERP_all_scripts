#!/bin/bash

#SBATCH --account=beetlegenes
#SBATCH --time=01:50:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sholtz1@uwyo.edu
#SBATCH --job-name=make_tree_big






module load swset/2018.05
module load gcc/7.3.0
module load r/3.6.1

gff_name=Tribolium_annotation.gff 
project_path=/project/beetlegenes/sholtz1/GERP
tree="(((((((((Tribolium_madens_genomic,Tribolium_castaneum),tribolium_confusum_genomic),Asbolus_verrucosus_genomic),Tenebrio_molitor_genomic),(Hycleus_cichorii_genomic,Pyrochroa_serraticornis_genomic)),Anoplophora_glabripennis_genomic),(Hypothenemus_hampei_genomic,Dendroctonus_ponderosae_genomic)),Agrilus_planipennis_genomic),Drosophila_melanogaster_genomic);"
Rscript /project/beetlegenes/sholtz1/GERP/gerp/scripts-gerp/estimate_neutral_tree.R ${project_path}/data/sequences/ref/split/ ${project_path}/analyses/last/msa_fasta/ ${project_path}/data/annotations/${gff_name} $tree ${project_path}/analyses/tree/neutral_tree.txt








