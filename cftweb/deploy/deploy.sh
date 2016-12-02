#!/bin/bash
# Run this script to create a virtual machine on which to host a docker image.
#
# We are using this VM to host our CFT web server app.
# The server is hosted in a docker container.
#
# Commands in this script are adapted from the instructions for "How to run your own Docker host"
# https://teams.fhcrc.org/sites/citwiki/SciComp/Pages/How%20to%20run%20your%20own%20Docker%20host.aspx
#

set -eux
shopt -s nullglob

prox new --runlist runlists/docker.runlist --disk 8 --no-bootstrap  cftweb4


