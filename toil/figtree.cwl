#!/usr/bin/env cwl-runner
#
# usage: (all are equivalent)
#   figtree.cwl sample.yml
#   cwl_toil figtree.cwl sample.yml
#   cwl_runner figtree.cwl sample.yml
#   cwl_tool figtree.cwl sample.yml
#
# where sample.yml contains - 
#       treefile:
#         class: File
#         path: sample_dedup.nwk
#       svgfile: sample_dedup.svg
#


cwlVersion: v1.0
class: CommandLineTool

baseCommand: java
arguments: ['-client', '-Djava.awt.headless=true', '-Xms64m',
'-Xmx512m', '-jar', $(inputs.jarfile.path),
               '-graphic', 'SVG']
inputs:
  treefile:
    type: File
    # description: "tree file in newick format, e.g. output from FastTree"
    inputBinding:
      position: 2
  svgfile:
    type: string
    # description: "tree file in newick format, e.g. output from FastTree"
    inputBinding:
      position: 3
  jarfile:
    type: File
    # description: "path to figtree JAR file."
    default:
      class: File
      location: "/home/matsengrp/local/lib/figtree.jar"
  width:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "-width"
      separate: true
    default: 800
  height:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "-height"
      separate: true
    default: 600
outputs:
  svgout:
    type: File
    outputBinding:
      glob: $(inputs.svgfile)
