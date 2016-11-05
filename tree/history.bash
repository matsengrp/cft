#! /bin/env bash

# command line inputs
#--------------------
# partis output fasta
fasta=$1
# output dir
outdir=$2
# naive name
naive=$3
# seed name
seed=$4
#--------------------

# clusters within *all* file can have different lengths
# let's Npad left of shorter ones (NOTE: assuming that is the correct alignment)

# file name for temp Npadded
Npad=${outdir}/`basename $fasta`.Npad.fa
# the longest sequence length
maxL=`awk '{if(substr($0,0,1)!=">") {print length}}' $fasta | sort -n | uniq | tail -1`
# left pad shorter sequences with Ns
awk -v maxL=$maxL '{if(substr($0,0,1)!=">"){diff=maxL-length; if(diff>0){for(c=0;c<diff;c++) printf "N"}print $0} else print $0}' $fasta > $Npad

# convert fasta data to phylip
# note that names will be trimmed to 10 characters
phylip=${outdir}/`basename $fasta`.phy
module load seqmagick
seqmagick convert $Npad $phylip

rm $Npad

# get the line number of the naive, for outgroup rooting in dnaml
naiveLine=`tail -n +2 $phylip | grep -n $naive | cut -f1 -d ':'`

# run phylip's dnaml, producing new outfile and outtree
dnaml <<STDIN
$phylip
O
$naiveLine
5
.
Y

STDIN

# throw away phylip's newick tree and the input file to phylip
rm outtree $phylip
# move phlip outfile
outfile=${outdir}/`basename $fasta .fa`.outfile2tree
mv outfile $outfile

# trim seed, since it was trimmed to 10 characters in converting to phylip
seed=`echo $seed | awk '{print substr(0,10)}'`
# parse outfile and make trees
python `dirname $0`/outfile2tree.py --outfile $outfile --naive $naive --seed $seed

# throw away the phylip outfile
rm $outfile
