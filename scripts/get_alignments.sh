#!/bin/bash
find data/raw_data -maxdepth 1 -type f -name "*_seqs.fasta" -print0 | xargs -0 -n1 -I{} sh -c '
    in="{}"
    sample=$(basename "$in" | cut -d _ -f 1)
    mafft --auto "$in" > "data/alignments/${sample}_alignment.fasta"
'

