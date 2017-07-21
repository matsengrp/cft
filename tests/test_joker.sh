#!/usr/bin/env bash

datadir=output/test-laura-mb-v9-minadcl-raxml/BF520/BF520.3-igl/BF520-l-IgL/run-viterbi-best-plus-3/cluster0/
PATH=./bin:$PATH

echo "Datadir is: $datadir"


echo $datadir


./bin/name_internal_nodes.py --splits $datadir/asr.nwk prankin.nwk
./bin/name_internal_nodes.py --splits $datadir/asr_seqs.best.dnd prankout.nwk

nw_labels prankin.nwk | sort | uniq > prankin.ids
nw_labels prankout.nwk | sort | uniq > prankout.ids


echo "prank input nodes:"
wc -l prankin.ids
echo "prank output nodes:"
wc -l prankout.ids

# These are different... problem
diff prankin.ids prankout.ids > prankinout.iddiff
#cat prankinout.iddiff

echo "Diff lines: $(wc -l prankinout.iddiff)"


echo "\nNow we have to check to see if it's just the nw_reroot"
nw_reroot $datadir/asr.nwk > rerooted.nwk

./bin/name_internal_nodes.py --splits rerooted.nwk prankin.rrtd.nwk
nw_labels prankin.rrtd.nwk | sort | uniq > prankin.rrtd.ids
echo "n nodes"
wc -l prankin.rrtd.ids

echo "Diff (should be empty):"

diff prankin.rrtd.ids prankout.ids > prankinout.rrtd.iddiff
cat prankinout.rrtd.iddiff


echo "\nOk, now we're going to test out whether we get the same ordering in the prank output"
#nw_labels prankout


# something like this... they diff to null; so 
#cat $datadir/asr_seqs.best.dnd | nw_labels -I - > out.labels
#nw_reroot $datadir/asr.nwk | nw_labels -I - > in.rrtd.labels
#diff in.rrtd.labels out.labels


# So... the splits diff to null, and the label orderings are the same, so we can concur that our rerooted
# input has labels in the right order.
# So now we must need to iterate through the internal nodes in the sequence file and pair them up with the
# labels 


