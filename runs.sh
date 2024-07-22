#!/usr/bin/env bash

DATADIR="/data/sequencing/nextseq-investigation"

# Clean up previous run directory
test -d nextseq-inv && rm -rf nextseq-inv

# Run the Python script and capture its output
CSV_CONTENT=$(python3 ~/.local/bin/get-sequencer --directory "$DATADIR")

echo -e $CSV_CONTENT

# Read the CSV content and process each folder
echo -e "$CSV_CONTENT" | while IFS=',' read -r folder instrument color_scheme; do
    # Skip the header line
    if [ "$folder" == "Folder" ]; then
        continue
    fi

    echo "Processing ${folder}..."

    # Determine the sequencer based on the color scheme
    if [ "$color_scheme" == "4-color" ]; then
        sequencer="miseq"
    elif [ "$color_scheme" == "2-color" ]; then
        sequencer="nextseq"
    else
        echo "Unknown color scheme for ${folder}, skipping..."
        continue
    fi

    # Run the Nextflow command
    nextflow -q run main.nf \
        --readDIR "${DATADIR}/${folder}" \
        -profile docker \
        --target v4 \
        -c conf/custom.config \
        --QC_only \
        -w "nextseq-inv/work_${folder}" \
        --outDIR "nextseq-inv/results_${folder}" \
        --sequencer "$sequencer"

done
