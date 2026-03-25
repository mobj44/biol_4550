from amas import AMAS
from glob import glob
from pathlib import Path
import argparse


def main():
    # set working directory
    wd = Path(__file__).resolve().parents[1]
    alignment_dir = (wd / 'data'/'trimmed_alignments').relative_to(wd)

    parser = argparse.ArgumentParser(
        prog='Make Supermatrix',
        description='Makes supermatrix from fasta files'
    )

    parser.add_argument('-c', '--cores', default=4)
    args = parser.parse_args()

    files = sorted(glob(f'{alignment_dir}/*.fasta'))

    multi_meta_aln = AMAS.MetaAlignment(
        in_files=files,
        data_type="dna",
        in_format="fasta",
        cores=args.cores
    )

    # parse each alignment into a dictionary
    parsed_algnments = multi_meta_aln.get_parsed_alignments()

    # concatenate dictionaries into supermatrix + partitions
    supermatrix, partitions = multi_meta_aln.get_concatenated(
        parsed_algnments)

    supermatrix_dir = wd / 'data'/'supermatrix'
    supermatrix_dir.mkdir(parents=True, exist_ok=True)

    nexus_str = multi_meta_aln.print_nexus(supermatrix)

    with open(supermatrix_dir / 'supermatrix.nex', 'w') as out:
        out.write(nexus_str)

    fasta_str = multi_meta_aln.print_fasta(supermatrix)

    with open(supermatrix_dir/'supermatrix.fasta', 'w') as out:
        out.write(fasta_str)

    # write partitions file for IQ-TREE (clean names + DNA prefix)
    with open(supermatrix_dir / "partitions.txt", "w") as pf:
        for gene_name, coord_range in partitions.items():

            # Handle range objects vs strings
            if hasattr(coord_range, "start") and hasattr(coord_range, "stop"):
                start = coord_range.start + 1
                end = coord_range.stop
            else:
                s = str(coord_range).strip()
                s = s.split(",")[0].split("\\")[0]
                start_str, end_str = s.split("-")
                start, end = int(start_str), int(end_str)

            # Clean locus name for IQ-TREE
            # remove path + extension
            locus = Path(gene_name).stem
            locus = locus.replace("_trimmed", "")      # remove trimming suffix
            locus = locus.replace("_alignment", "")    # optional cleanup

            pf.write(f"DNA, {locus} = {start}-{end}\n")


if __name__ == '__main__':
    main()
