#!/bin/bash

INPUT_R1=$1
INPUT_R2=$2
OUTPUT_PREFIX=$3
DATABASE_PATH=$4

# Check if required arguments are provided
if [ -z "$INPUT_R1" ] || [ -z "$INPUT_R2" ] || [ -z "$OUTPUT_PREFIX" ] || [ -z "$DATABASE_PATH" ]; then
    echo "Usage: $0 <input_R1.fastq.gz> <input_R2.fastq.gz> <output_prefix> <database_path>"
    exit 1
fi

# Ensure output directory exists
OUTPUT_DIR="/data/krona_kaiju_files"
mkdir -p "$OUTPUT_DIR"

# Step 1: Run Kaiju for virus-specific classification
echo "Running Kaiju for virus-specific classification..."
kaiju -t "$DATABASE_PATH/nodes.dmp" \
      -f "$DATABASE_PATH/kaiju_db_viruses.fmi" \
      -i "$INPUT_R1" \
      -j "$INPUT_R2" \
      -o "$OUTPUT_DIR/${OUTPUT_PREFIX}_viruses.kaiju"

if [ $? -ne 0 ]; then
    echo "Error: Kaiju classification failed."
    exit 1
fi

# Step 2: Convert Kaiju output to Krona-compatible format
echo "Converting Kaiju output to Krona format..."
kaiju2krona -t "$DATABASE_PATH/nodes.dmp" \
            -n "$DATABASE_PATH/names.dmp" \
            -i "$OUTPUT_DIR/${OUTPUT_PREFIX}_viruses.kaiju" \
            -o "$OUTPUT_DIR/${OUTPUT_PREFIX}_viruses.krona"

if [ $? -ne 0 ]; then
    echo "Error: Conversion to Krona format failed."
    exit 1
fi

# Step 3: Generate Krona HTML visualization
echo "Generating Krona HTML visualization..."
perl /opt/KronaTools-2.8.1/scripts/ImportText.pl \
     "$OUTPUT_DIR/${OUTPUT_PREFIX}_viruses.krona" \
     -o "$OUTPUT_DIR/${OUTPUT_PREFIX}_viruses.html"

if [ $? -ne 0 ]; then
    echo "Error: Krona HTML generation failed."
    exit 1
fi

echo "Analysis completed. Results are in $OUTPUT_DIR."
