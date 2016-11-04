
Bootstrap two-column template downloaded from https://startbootstrap.com/template-overviews/2-col-portfolio/

## Usage

```
    $ cd cft/cftweb
    $ python -m cftweb -c "sampledata"
```

The "database" of clusters is read from *any* JSON file imediately under the directory supplied to the `-c` parameter.
The format of the cluster database is a array of dictionaries, one per cluster, e.g.

```
[
    { 
	"fasta": "sample1.fa",
	"tree": "sample1.tree",
	"svg": "sample1.svg"
    }
]
```

See `sampledata/sample.json` for an example. All filenames in the index file are relative to the directory in which the database index is found.

