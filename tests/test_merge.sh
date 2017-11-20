#!/usr/bin/env bash

../bin/merge_timepoints_and_multiplicity.py --cluster-mapping test_merge/cluster_mapping.csv \
  --partis-seqmeta test_merge/partis_seqmeta.csv \
  --upstream-seqmeta test_merge/seqmeta.csv \
  test_merge/results.csv

diff test_merge/expected_results.csv test_merge/results.csv


