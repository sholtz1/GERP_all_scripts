---
title: "NAM conservation analysis readme"
author: "Asher Hudson"
date: "7/20/2020"
output:
  html_document: default
  pdf_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, eval = FALSE)
```

# Conservation analysis pipeline

### Key steps

1. Align query genomes to maize B73 v5 reference genome.
2. Use GERP to identify conserved elements.
3. Statistical analysis to detect enrichment/depletion of structural variants in conserved elements.


### Software used for these analyses *(add versions)*

a) Last (lastdb, last-train, lastal, maf-convert, last-split, maf-swap, last-postmask) v 2.31.1
b) Kent-utils (axtChain, chainMergeSort, faSize, chainPreNet, chainNet, faToTwoBit, netToAxt, axtToMaf)
c) multiz (multiz.v10.6)
d) GERP++
e) bedtools 2.27.1
f) R 3.6.3
g) R packages
    i) rPHAST_1.6.9
    ii) dplyr_0.8.5
    iii) data.table_1.12.6
    iv) RColorBrewer_1.1-2
    v) tidyr_1.0.0
    vi) ggplot2_3.3.0
    vii) grid_3.5.1

## Preface: Set up file system and input files
```{bash}
project_path=
mkdir -p ${project_path}/data/sequences/ref/
mkdir -p ${project_path}/data/sequences/query/
mkdir -p ${project_path}/data/variants/
mkdir -p ${project_path}/data/chromsize/
mkdir -p ${project_path}/data/annotations/
```

Download the following files and add them to the following directories. All except query files available on Cyverse: \
`${project_path}/data/sequences/ref/` \
B73 v5 reference sequence: Zm-B73-REFERENCE-NAM-5.0.fa \
Repeat masked B73 v5 reference sequence: Zm-B73-REFERENCE-NAM-5.0_rm.fa \
`${project_path}/data/sequences/query/` \
From Phytozome: \
Acomosus: Acomosus_321_v3.softmasked.fa.gz \
Osativa: Osativa_323_v7.0.softmasked.fa.gz \
Taestivum: Taestivum_296_v2.softmasked.fa.gz \
Bdistachyon: Bdistachyon_556_v3.0.softmasked.fa.gz \
Sbicolor: Sbicolor_313_v3.0.softmasked.fa.gz \
Othomaeum: Othomaeum_386_v1.0.fa.gz \
Stuberosum: Stuberosum_448_v4.03.fa.gz \
Vvinifera: Vvinifera_457_Genoscope.12X.fa.gz \
From Ensembl: \
Athaliana: Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz \
We used four additional genomes which are currently unreleased. Chrysopogon serrulatus and Hyparrhenia diplandra (Song et al. 2020 preprint), and genomes for Coelorachis and Vossia.\
`${project_path}/data/variants/` \
NAM founders SVs in a bed file: NAM_founders_SVs.bed \
NAM founders SNPs in a bed file: B73v5.NAM-illumina-snps-only_filtered.frq.bed \
`${project_path}/data/chromsize/` \
Size of all chromosomes in B73: Zm-B73-REFERENCE-NAM-5.0.size.bed \
`${project_path}/data/annotations/` \
B73 v5 genes: zea_maysb73_core_3_87_1.gff \
Open chromatin from merging leaf and inflorescence data from Ricci et al. 2019: open_chromatin.bed \
GFF3 of repeat masked bases in B73 v5: Zm-B73-REFERENCE-NAM-5.0-repeatmask.gff3.gz \
Genetic map from Ogut et al. 2015 lifted over to AGPv5: ogut_fifthcM_map_agpv5.bed\

Run all of the following scripts from within your `${project_path}`.\

## 1) Align genomes

### Step 1: Run lastdb.sh on the reference genome

```{bash, eval=FALSE}
bash ./scripts-alignment/lastdb.sh
```

```{bash}

ref=./data/sequences/ref/Zm-B73-REFERENCE-NAM-5.0.fa
ref_name=$( basename $ref )
db_prefix=./analyses/last/lastdb/${ref_name}-MAM4

mkdir -p ./analyses/last/lastdb/
lastdb -P0 -uMAM4 -R01 $db_prefix $ref
```

### Step 2: Make alignments

Run the following for each query genome in ${project_path}/data/sequences/query/
```{bash, eval = FALSE}
bash make_alignments.sh query
```
```{bash, eval=FALSE}

ref=./data/sequences/ref/Zm-B73-REFERENCE-NAM-5.0.fa
ref_name=$( basename $ref )
db_prefix=./analyses/last/lastdb/${ref_name}-MAM4
query=$1
query_name=$(basename $query .gz)
query_name=$(basename $query_name .fa)

mkdir -p ./analyses/last/mat/
mat=./analyses/last/mat/"$ref_name"_"$query_name".mat

last-train -P0 --revsym --matsym --gapsym -E0.05 -C2 $db_prefix $query > $mat

mkdir -p ./analyses/last/maf/

maf=./analyses/last/maf/Zm-B73-REFERENCE-NAM-5.0_"$query_name".maf

lastal -m50 -E0.05 -C2 -p $mat $db_prefix $query > $maf
# m50 makes allowed multiplicity of initial hits 50, each match lengthened
# until it occurs at most this many times
# E is threshold for errors
# C2 makes lastal discard any gapless alignment in 2 or more gapped alignments
# p specifies match/mismatch score matrix
# Then reference db
# And finally the query fasta


mkdir -p ./analyses/last/axt
axt=./analyses/last/axt/"$ref_name"_"$query_name".axt

maf-convert axt $maf > $axt

mkdir -p ./analyses/last/chain
chain=./analyses/last/chain/"$ref_name"_"$query_name".chain
merged_chain=./analyses/last/chain/"$ref_name"_"$query_name".all.chain

axtChain $axt $ref $query $chain -linearGap=loose -faQ -faT

# axtChain options
# linearGap=loose for distant species
# -faQ The specified qNibDir is a fasta file with multiple sequences for query
# -faT The specified tNibDir is a fasta file with multiple sequences for target

mkdir -p ./analyses/last/chain_merged
merged_chain=./analyses/last/chain/"$ref_name"_"$query_name".all.chain

chainMergeSort $chain > $merged_chain

mkdir -p ./analyses/chromsize/
ref_size=./analyses/chromsize/"$ref_name".size
query_size=./analyses/chromsize/"$query_name".size

if [ ! -f $ref_size ]; then
  faSize $ref -detailed > $ref_size
fi

if [ ! -f $query_size ]; then
  faSize $query -detailed > $query_size
fi

# note: if your chromosomes are id'd with numbers, faSize will add everything
# before .fa in the file name to the beginning of each chromosome id in the size
# file. This causes issues in the next step. Either make all id's non-numeric
# (e.g. 'chr1' instead of '1') or go in to size file and remove prefix in front
# of chromosome ids.

mkdir -p ./analyses/last/chain_prenet/
chain_prenet=./analyses/last/chain_prenet/"$ref_name"_"$query_name".all.pre.chain
chainPreNet $merged_chain $ref_size $query_size $chain_prenet

mkdir -p ./analyses/last/target_net/
mkdir -p ./analyses/last/query_net/
target_net=./analyses/last/target_net/"$ref_name"_"$query_name".net
query_net=./analyses/last/query_net/"$query_name"_"$ref_name".net
chainNet $chain_prenet $ref_size $query_size $target_net $query_net

# making 2bit files for netToAxt
mkdir -p ./analyses/last/2bit/
query_twobit=./analyses/last/2bit/"$query_name".2bit
ref_twobit=./analyses/last/2bit/"$ref_name".2bit

if [ ! -f $query_twobit ]; then
  faToTwoBit $query $query_twobit
fi

if [ ! -f $ref_twobit ]; then
  faToTwoBit $ref $ref_twobit
fi

mkdir -p ./analyses/last/net_axt/
net_axt=./analyses/last/net_axt/"$ref_name"_"$query_name".net.axt
netToAxt $target_net $chain_prenet $ref_twobit $query_twobit $net_axt

mkdir -p ./analyses/last/net_axt/net_maf
net_maf=./analyses/last/net_axt/net_maf/"$ref_name"_"$query_name".net.maf

axtToMaf $net_axt $ref_size $query_size $net_maf

# now to make mafs one to one

head -n 29 $maf > "$net_maf"_w_header
cat $net_maf >> "$net_maf"_w_header

one_maf=./analyses/last/net_axt/net_maf/"$ref_name"_"$query_name".1to1.maf

last-split -m1 "$net_maf"_w_header |
maf-swap |
awk -v q="$query_name" -v r="$ref_name" '/^s/ {$2 = (++s % 2 ? q "." : r ".") \
$2} 1' | last-split -m1 | \
maf-swap | last-postmask > $one_maf

```

### Step 3: Combine and filter maf files

```{bash, eval=FALSE}
bash ./scripts-alignment/process_mafs.sh
```

```{bash}

maf_array=($( ls -d ./analyses/last/net_axt/net_maf/*1to1.maf ))
combined_maf=./analyses/last/net_axt/net_maf/combined.maf

cat ${maf_array[@]:0:1} | sed -n '/##maf version=1 scoring=blastz/,$p' > \
${maf_array[@]:0:1}_tmp
cat ${maf_array[@]:1:1} | sed -n '/##maf version=1 scoring=blastz/,$p' > \
${maf_array[@]:1:1}_tmp

multiz ${maf_array[@]:0:1}_tmp ${maf_array[@]:1:1}_tmp 1 > $combined_maf

for maf in ${maf_array[@]:2};
do
  echo "processing " $maf
 cat $maf | sed -n '/##maf version=1 scoring=blastz/,$p' > \
 "$maf"_tmp
 multiz $combined_maf "$maf"_tmp 1 > "$combined_maf"_tmp
 mv "$combined_maf"_tmp $combined_maf
done

# and filter mafs so all blocks have Zea mays and are at least 20 bp long
mafFilter -minCol=20 -needComp="$ref_name" $combined_maf > "$combined_maf".filtered
```

### Step 4: Convert maf files to fasta files

```{bash}
bash ./scripts-alignment/maf_to_msa_fasta.sh
```

### We need to replace all *'s in the file with N's to get GERP to read it. We also need to change the first line in the files to Tribolium_castaneum to match the reference file. 

```{bash}
bash replace_stars.sh 
```



```{bash}

ref_rm=./data/sequences/ref/Zm-B73-REFERENCE-NAM-5.0_rm.fa

ref_rm_name=$( basename $ref_rm )

combined_maf=./analyses/last/net_axt/net_maf/combined.maf

# splitting combined maf by target sequence
mkdir -p ./analyses/last/split_maf/
outroot=./analyses/last/split_maf/
mafSplit -byTarget dummy.bed $outroot "$combined_maf".filtered -useFullSequenceName

# mafSplit throws an error if you don't put dummy.bed, even though byTarget
# means it's being ignored

path_to_phast=~/bin/phast
msa_view=${path_to_phast}/bin/msa_view
maf_dir=($( ls -d ./analyses/last/split_maf/*chr*.maf ))


mkdir -p ./analyses/last/msa_fasta/
mkdir -p ./data/sequences/ref/split/

faSplit byname $ref_rm ./data/sequences/ref/split/

path_to_match_masking=./scripts-alignment/matchMasking.pl

for maf_file in "${maf_dir[@]}"
do
  chr=$(basename $maf_file | sed -e 's/0\(.*\).maf/\1/')

  fasta=./analyses/last/msa_fasta/"$chr".fa
  ref_rm_chr=./data/sequences/ref/split/"$chr".fa
  rm_fasta=./analyses/last/msa_fasta/"$chr"_rm.fa
  $msa_view $maf_file -f -G 1 --refseq $ref_rm_chr > $fasta
  sed -E 's/> />/g' $fasta > "$fasta"_tmp && mv "$fasta"_tmp $fasta
  perl $path_to_match_masking \
  --ref $ref_rm_chr \
  --fasta $fasta \
  --out $rm_fasta
done
```



## 2) Identify conserved elements

### Step 5: Estimate neutral tree

```{bash}
project_path=
tree="((Acomosus_321_v3,((Osativa_323_v7,(Taestivum_296_v2,Bdistachyon_556_v3)),((Vossia,(Zm-B73-REFERENCE-NAM-5,(Coelorachis,(Chrysopogon_serrulatus_FLYE2,(Hyparrhenia_diplandra_flye2,Sbicolor_454_v3))))),Othomaeum_386_v1))),(Solanum_tuberosum,(Vitis_vinifera,Arabidopsis_thaliana)));"
Rscript ./scripts-gerp/estimate_neutral_tree.R ${project_path}/data/sequences/ref/split/ ${project_path}/analyses/last/msa_fasta/ ${project_path}/data/annotations/${gff_name} $tree ${project_path}/analyses/tree/neutral_tree.txt

```

```{r}
library("rphast")
args <- commandArgs(trailingOnly = TRUE)
ref_root <- args[1]
msa_root <- args[2]
feat_file <- args[3]
tree <- args[4]
tree_output <- args[5]

chr_vector <- paste("chr", c(1:10), sep = "")

chr <- chr_vector[1]

ref_file <- paste(ref_root, chr, ".fa", sep = "")

msa <- paste(msa_root, chr, ".fa", sep="")

align4d <- read.msa(msa, refseq = ref_file, format = "FASTA", do.4d = TRUE, features= read.feat(feat_file))

for (chr in chr_vector[2:length(chr_vector)]){
  ref_file <- paste(ref_root, chr, ".fa", sep = "")
  msa <- paste(msa_root, chr, ".fa", sep="")
  align4d2 <- read.msa(msa, refseq = ref_file, format = "FASTA", do.4d = TRUE, features= read.feat(feat_file))
  align4d <- concat.msa(list(align4d, align4d2))
}

neutralMod <- phyloFit(align4d, tree=tree, subst.mod="REV")
sink(tree_output)
neutralMod$tree
sink()
```

### Step 6: Run GERP

Run for all 10 maize chromosomes (so chr1, chr2, etc.)
```{bash}
bash ./scripts-gerp/run_gerp.sh chr
```

run_gerp.sh
```{bash}

chr=$1
mkdir -p {project_path}/analyses/gerp/
fasta=./analyses/last/msa_fasta/"$chr"_rm.fa


ref_rm=./data/sequences/ref/Zm-B73-REFERENCE-NAM-5.0_rm.fa
ref_rm_chr=./data/sequences/ref/split/"$chr".fa
gerpcol={path/to/GERPplusplus/gerpcol}
gerpelem={path/to/GERPplusplus/gerpelem}
tree=./analyses/tree/neutral_tree.txt
$gerpcol -f $rm_fasta -t $tree -v -e Zm-B73-REFERENCE-NAM-5 -j -a


mkdir -p ./analyses/gerp/
gerp=./analyses/gerp/"$chr"_rm.rates
mv "$rm_fasta".rates $gerp

$gerpelem -f $gerp

# # make bed file
#
awk -v chr="$chr" 'BEGIN {OFS="\t"}; {print chr, NR-1, NR, $1, $2}' $gerp > "$gerp".bed
# #
# # # bed file of positive scores
awk '$5 > 0' "$gerp".bed > "$gerp".pos.bed
awk -v chr="$chr" 'BEGIN {OFS="\t"} {print chr, $2-1, $3, $4, $5, $6, $7, $8}' \
"$gerp".elems > "$gerp".elems.bed
sort -k 1,1 -k2,2n "$gerp".elems.bed > "$gerp".elems.bed_tmp
mv "$gerp".elems.bed_tmp "$gerp".elems.bed
```

###Final output is a .rates.bed file and a .rates.elems file for each chromosome that can then be combined with other files. In this case we combine them with allele frequency and fitness files from the tribolium range expansion experiemnt. 