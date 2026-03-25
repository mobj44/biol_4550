# Data Folder Organization

| Folder                    | Description                                                                                                                                                           |
| :------------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **`raw_data/`**           | Contains the raw FASTA sequence files downloaded from NCBI using the script `sctipts/get. Each file corresponds to a single gene across al taxa.                      |
| **`alignments/`**         | Multiple sequence alignments (MSAs) generated from the raw FASTA files using MAFFT. Each file represents an aligned gene dataset.                                     |
| **`trimmed_alignments/`** | Cleaned and trimmed versions of the MAFFT alignments. Gaps and poorly aligned regions were filtered out using R scripts.                                              |
| **`supermatrix/`**        | Contains the concatenated supermatrix alignment and partition files built from the trimmed alignments. Input for both maximum likelihood and Bayesian tree inference. |
| **`tables/`**             | CSV tables summarizing metadata for each gene and taxon (e.g., sequence length, accession, and taxonomy), generated during the data retrieval step.                   |
