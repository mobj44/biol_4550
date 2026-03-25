# Steps

1. Get Accession Numbers from Genbank. (`./scripts/get_tables.py`)
1. Create table with all accession numbers. (`./scripts/get_master_table.py`)
1. Getting Data

   - Get fasta data from Genbank. (`./scripts/get_fastas.py`)

1. Running initial alignments (`./scripts/get_alignments.sh`)
1. Trimming Alignments (`./scritps/trim_alignments.R`)
1. Create Supermatrix (`./scripts/supermat.R`)
1. Trees:
   - ML trees (see `./scripts/README.md` for more info)
   - tree files are output to `./trees`
1. Tree Visualization
   - ML tree (see `./scripts/tree_vis.R`)
