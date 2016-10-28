#!/usr/bin/env cwl-runner

# workflow for making trees form fasta files
# run with,
#     $ cwltoil --preserve-environment LD_LIBRARY_PATH PATH --no-container  workflow.cwl workflow.yml
#
# Expects a yaml file to configure input variables as follows,
#    fasta:
#      class: File
#      path: sample_dedup.fa
#    svgfile: sample_dedup.svg


cwlVersion: cwl:v1.0

class: Workflow
requirements:
  - class: InlineJavascriptRequirement
inputs:
  fasta:
    type: File?
  svgfile:
    type: string
    doc: "name of SVG output file."

outputs:
  svgfile:
    type: File
    outputSource: figtree/svgout

steps:
  fasttree:
    in:
      fasta: fasta
    out: [output]
    run: fasttree.cwl
    
  figtree:
    in:
      treefile: fasttree/output
      svgfile: svgfile
    out: [svgout]
    run: figtree.cwl
    

