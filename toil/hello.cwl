cwlVersion: v1.0
class: CommandLineTool
baseCommand: echo
inputs:
  message:
    type: string
    inputBinding:
      position: 1
outputs: 
  output:
    type: stdout
#stdout: $(inputs.message).ext
stdout: results.txt

