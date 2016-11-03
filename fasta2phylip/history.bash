#! /bin/env bash

# partis output fasta
fasta=$1
# optional output dir
outdir=$2

# fasta data to phylip
phylip=${outdir}/`basename ${fasta}`.phylip
seqmagick convert $fasta $phylip

echo -n "computing ML tree... "
# run phylip's dnaml, producing new outfile and outtree
dnaml <<STDIN 1> /dev/null
`pwd`/${fasta}.phylip
O
1
5
.
Y

STDIN
echo "done"

# parse outfile to extract all sequences
python outfile2tree.py > ${phylip}.outfile.extracted_sequences.fa

# move dnaml files to output directory
mv outfile ${phylip}.outfile.txt
mv outtree ${phylip}.outtree.txt
