#!/bin/bash

# CONFIG
JUICER_JAR="/home/sutripa/kajal/AFTER_FINAL_GENOME_DEVELOPMENT/A_B_compart/juicer_latest/juicer-1.6/scripts/common/juicer_tools.jar"
HIC_FILE="sorted_clean_P.capsici_isolate_WB2_Final_HiC.hic"
RESOLUTION=10000
OUTPUT_DIR="oe_matrices"
NORM="KR"

mkdir -p $OUTPUT_DIR

# List of chromosomes
CHROMS=(HiC_scaffold_{1..17})

# Extract observed contact matrices
for (( i=0; i<${#CHROMS[@]}; i++ )); do
  for (( j=i+1; j<${#CHROMS[@]}; j++ )); do
    chr1=${CHROMS[$i]}
    chr2=${CHROMS[$j]}
    outfile="${OUTPUT_DIR}/${chr1}_${chr2}_obs.txt"
    echo "Extracting OBSERVED: $chr1 vs $chr2"
    java -jar $JUICER_JAR dump observed $NORM $HIC_FILE $chr1 $chr2 BP $RESOLUTION $outfile
  done
done

# Extract expected vectors
for chr in "${CHROMS[@]}"; do
  outfile="${OUTPUT_DIR}/${chr}_exp.txt"
  echo "Extracting EXPECTED: $chr"
  java -jar $JUICER_JAR dump expected $NORM $HIC_FILE $chr BP $RESOLUTION $outfile
done

