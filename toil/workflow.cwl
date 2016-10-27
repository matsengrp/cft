#!/usr/bin/env cwl-runner

# single-step workflow
# run with $ cwltoil --preserve-environment LD_LIBRARY_PATH PATH --no-container  workflow.cwl workflow.yml

cwlVersion: cwl:v1.0

class: Workflow
requirements:
  - class: InlineJavascriptRequirement
inputs:
  fasta:
    type: File?
    inputBinding:
      position: 1

outputs:
  nwkout:
    type: File
    outputSource: fasttree/output

steps:
  fasttree:
    in:
      fasta: fasta
    out: [output]
    run: fasttree.cwl
    

