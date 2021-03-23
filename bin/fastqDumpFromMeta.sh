#!/bin/bash

## metadata file is useful to identify samples and label after import

# extract first column 
for sra in $(cut -d, -f1 "SraRunTable.txt.csv" | tail -n +2); do

  echo "$sra"
  # --split-files: split paired-end fastq files into respective read files
  # --gzip: gzip compress upon fastq generation
  # --origfmt: bug catch for gzip problems
  fastq-dump --split-files --origfmt --gzip "$sra"

done
