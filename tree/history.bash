#! /bin/env bash

# partis output fasta
fasta=$1
# optional output dir
outdir=$2
# optional naive name
naive=$3
#optional seed name
seed=$4

# fasta data to phylip
phylip=${outdir}/`basename ${fasta}`.phylip
module load seqmagick
seqmagick convert $fasta $phylip

echo -n "computing ML tree... "
# run phylip's dnaml, producing new outfile and outtree
dnaml <<STDIN
`pwd`/${fasta}.phylip
O
1
5
.
Y

STDIN
echo "done"

# move dnaml files to output directory
mv outfile ${phylip}.outfile.txt
mv outtree ${phylip}.outtree.txt

# parse outfile to extract all sequences
python outfile2tree.py --outfile ${phylip}.outfile.txt --naive $naive --seed $seed
