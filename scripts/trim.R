library(ape)
library(phangorn)

process_alignment <- function(gene,
                              input_dir = "data/alignments",
                              output_dir = "data/trimmed_alignments",
                              trim_threshold = 0.8) {
    # Input and output file paths
    input_file <- file.path(input_dir, paste0(gene, "_alignment.fasta"))
    output_file <- file.path(output_dir, paste0(gene, "_trimmed.fasta"))

    # Check input
    if (!file.exists(input_file)) {
        stop(paste("Input file not found:", input_file))
    }
    aln <- read.FASTA(input_file)


    # Remove duplicates
    aln_unique <- aln[!duplicated(names(aln))]

    # Convert to phyDat and trim
    aln_phy <- phyDat(aln_unique, type = "DNA")
    aln_trim <- aln_phy[, colMeans(as.character(aln_phy) == "-") < trim_threshold]

    # Convert back to DNAbin and write
    aln_dnabin <- as.DNAbin(aln_trim)
    write.FASTA(aln_dnabin, file = output_file)

    invisible(aln_dnabin)
}


batch_process_alignments <- function(input_dir = "data/alignments",
                                     output_dir = "data/trimmed_alignments",
                                     trim_threshold = 0.8) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    files <- list.files(input_dir, pattern = "_alignment\\.fasta$", full.names = TRUE)
    if (length(files) == 0) stop("No alignment files found in ", input_dir)

    for (file in files) {
        gene <- sub("_alignment\\.fasta$", "", basename(file))
        process_alignment(gene, input_dir = input_dir,
                          output_dir = output_dir,
                          trim_threshold = trim_threshold)
    }
}

batch_process_alignments(
    input_dir = "./data/alignments",
    output_dir = "./data/trimmed_alignments",
    trim_threshold = 0.8
)
