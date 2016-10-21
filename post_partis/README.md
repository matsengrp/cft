# Python modules for processing partis output

Run `bash scripts/demo.sh` for an example.

## TODO

### Selecting alternative cluster(s)

Occasionally the cluster given in the annotations file will not be the one we want.
Create a flag and handle this case (say: `--select_clustering *int*`) using partition file as input.

#### Problem

This is fine for outputting clusters, but we do not have the naive sequence directly at-hand.
We can get this from `./bin/partis --naive-vsearch --other-args`, and the `.log` file has the function call used to generate the data.

### Output cluster statistics

For queried clusters, output the cluster ID that contains an input seed and the count of clones in that particular cluster.

### Indels

Sequences in the annotations `.csv` have indels removed, but if a sequence is queried we may want the sequence with the indels reconstructed.
Do we want this as a single `.fa` file?

