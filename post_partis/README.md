# Python modules for processing partis output

Run `bash scripts/demo.sh` for an example.

## TODO

### Selecting alternative cluster(s)

Occasionally the cluster given in the annotations file will not be the one we want.
Create a flag and handle this case (say: `--select_clustering *int*`) using partition file as input.

#### Problem

This is fine for outputting clusters, but as far as I can tell there is no way to get the naive sequence in this case.
Perhaps this is not functionality we want just yet.

### Output cluster statistics

For queried clusters, output the cluster ID that contains an input seed and the count of clones in that particular cluster.

### Indels

Sequences in the annotations `.csv` have indels removed, but if a sequence is queried we may want the sequence with the indels reconstructed.
Do we want this as a single `.fa` file?

