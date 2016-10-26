cwlVersion: v1.0
class: CommandLineTool
baseCommand: [FastTree]
# baseCommand: [/usr/bin/ldd, /app/easybuild/software/FastTree/2.1.9-foss-2016b/bin/FastTree]
# baseCommand: [/usr/bin/env]
stdout: output.txt
inputs: 
  fasta:
    type: File?
    inputBinding:
      position: 1
outputs:
  output:
    type: stdout
