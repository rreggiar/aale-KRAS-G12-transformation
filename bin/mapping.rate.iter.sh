#!/bin/bash

input_dir='/public/groups/kimlab/exoRNA-biomarkers-panc/data'

for salmon_dir in $input_dir/*/*/*.ucsc.rmsk.salmon.*/logs; do

  mapping_rate=$(grep 'Mapping rate' $salmon_dir/salmon_quant.log | cut -d'=' -f2)

  output_line=$(echo $(echo $salmon_dir | cut -d'/' -f8,9) $mapping_rate | \
    sed 's/ /,/g;s/\//,/g')

  echo $output_line
done

#for star_dir in $input_dir/collected.patient.data/*/star.out/pass.2; do

#  dir_name=$(echo $star_dir | cut -d'/' -f8)

#  cat $star_dir/Log.final.out | grep 'reads' | sed -e 's/|//g' | sed -e 's|^|'"$dir_name"'\t|' | sed 's/\t/,/g'
  #exit

  #mapping_rate=$(grep 'Uniquely mapped reads' $star_dir/Log.final.out | \
  #  cut -d'|' -f2 | sed 's/\ /,/g') 

  #output_line=$(echo $(echo $star_dir | cut -d'/' -f8,9) $mapping_rate | \
  #  sed 's/ /,/g;s/\//,/g') # | cut -d, -f1,2,4)

  #echo $output_line

#done
