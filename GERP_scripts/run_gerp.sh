chr=$1
mkdir -p /project/beetlegenes/sholtz1/GERP/analyses/gerp/
fasta=./analyses/last/msa_fasta/"$chr"_rm.fa.test


gerpcol=/project/beetlegenes/sholtz1/GERP/GERPsoftware/gerp++/gerpcol
gerpelem=/project/beetlegenes/sholtz1/GERP/GERPsoftware/gerp++/gerpelem
tree=./analyses/tree/neutral_tree_extra_aligns.txt

$gerpcol -f $fasta -t $tree -v -e Tribolium_castaneum -j -a


mkdir -p ./analyses/gerp/
gerp=./analyses/gerp/"$chr"_rm.rates
mv "$fasta".rates $gerp

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
