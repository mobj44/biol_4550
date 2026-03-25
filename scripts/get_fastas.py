import pandas as pd
from utils import *
from Bio import Entrez as ez
from dotenv import load_dotenv
from pathlib import Path
import os
from tqdm import tqdm


def get_one_fasta(accession, scientific_name):
    '''gets and sets a fasta property for the TaxInfo class'''
    if pd.isna(accession) or accession == '':
        return None

    scientific_name = scientific_name.lower().replace(" ", "_")

    handle = ez.efetch(
        db="nuccore",
        id=accession,
        rettype="fasta",
        retmode="text"
    )
    lines = handle.readlines()

    text = [line for line in lines if line.strip(
    ) != '' and not line.startswith('>')]

    fasta = ''.join(text)
    handle.close()

    return f'>{scientific_name}\n{fasta}'


def get_fastas(acc_column, names_col, gene):
    fasta_list = []
    for acc, sci_name in tqdm(
            zip(acc_column, names_col),
            total=len(acc_column),
            desc=gene):
        if pd.isna(acc) or acc == '':
            continue
        fasta = get_one_fasta(acc, sci_name)
        if fasta:
            fasta_list.append(fasta)

    return fasta_list


def create_fasta(file_name: Path, fasta_list: list[str]):
    '''creates a multi-fasta per gene'''
    with open(file_name, 'w') as fasta:
        for seq in fasta_list:
            fasta.write(seq)


def main():
    '''This is the main function'''
    # set working directory
    wd = Path(__file__).resolve().parents[1]

    # load env variables for email and api key
    load_dotenv(wd / '.env')
    ez.email = os.getenv("ncbi_email")
    ez.api_key = os.getenv("ncbi_key")

    fasta_path = create_data_path(wd, "raw_data")
    acc_file = wd / 'data' / 'tables' / 'master_accessions.csv'

    if not os.path.exists(acc_file):
        print('ERROR: Accessions file not found.')
        return

    df = pd.read_csv(acc_file)

    for gene in df.columns[2:]:
        fastas = get_fastas(df[gene], df['scientific_name'], gene)
        create_fasta(fasta_path / f'{gene}_seqs.fasta', fastas)


if __name__ == "__main__":
    main()
