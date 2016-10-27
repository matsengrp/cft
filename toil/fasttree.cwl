#!/usr/bin/env cwl-runner
#
# usage: (all are equivalent)
#   fasttree.cwl sample.yml
#
# also:  (all are equivalent)
#   cwl_runner fasttree.cwl sample.yml
#   cwl_toil fasttree.cwl sample.yml
#   cwl_tool fasttree.cwl sample.yml
#
# where sample.yml contains - 
#       treefile:
#         class: File
#         path: sample_dedup.nwk
#       svgfile: sample_dedup.svg
#
cwlVersion: v1.0
class: CommandLineTool
requirements:
  - class: InlineJavascriptRequirement

baseCommand: [FastTree]
# baseCommand: [/usr/bin/ldd, /app/easybuild/software/FastTree/2.1.9-foss-2016b/bin/FastTree]
# baseCommand: [/usr/bin/env]
inputs: 
  fasta:
    type: File?
    inputBinding:
      position: 1
outputs:
  output:
    type: stdout

stdout: $(inputs.fasta.basename.replace(/\.[^/.]+$/, "") + ".nwk")


