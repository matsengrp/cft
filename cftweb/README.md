
CFT Web application

## Installation

Run `pip install -r requirements.txt` (add `--user` if you get permissions errors).

Then install miniconda and ete using ete's flavor of instructions [here](http://etetoolkit.org/download/).

## Usage

```
cd cft/cftweb
python -m cftweb /path/to/data/files/*/metadata.json
# Typically this will be the following:
#python -m cftweb ../output/*/metadata.json
```


