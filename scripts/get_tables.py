import time
from dataclasses import dataclass
from Bio import Entrez as ez
import json
from dotenv import load_dotenv
import os
import csv
from pathlib import Path
import argparse
from tqdm import tqdm
from utils import *


@dataclass
class Gene:
    name: str
    term: str
    file_prefix: str


@dataclass
class TaxInfo:
    tax_id: str
    scientific_name: str | None = None
    accession: str | None = None
    seq_name: str | None = None
    seq_len: str | None = None

    def get_sci_name(self):
        '''sets the scientific name of the TaxInfo class'''
        tax_sum = ez.read(ez.esummary(
            db="taxonomy",
            id=self.tax_id
        ))[0]
        self.scientific_name = tax_sum["ScientificName"]

    def nuc_processing(self, gene: str):
        '''gets and sets the accession, seq name, and seq length on the TaxInfo class'''

        nuc_search = entrez_retry(
            lambda: ez.read(
                ez.esearch(
                    db="nuccore",
                    term=f'txid{self.tax_id}[ORGN] AND {gene}',
                    retmax=1
                )
            )
        )

        nuc_id = nuc_search["IdList"][0] if len(
            nuc_search["IdList"]) >= 1 else None

        if nuc_id is None:
            return

        nuc_sum = entrez_retry(
            lambda: ez.read(
                ez.esummary(
                    db="nuccore",
                    id=nuc_id
                )
            )[0]
        )

        self.accession = nuc_sum["AccessionVersion"]
        self.seq_name = nuc_sum["Title"]
        self.seq_len = int(nuc_sum["Length"])


def entrez_retry(request_func, retries=3, delay=1.5):
    """
    Retry wrapper for Entrez calls: retry only on NCBI backend errors.
    request_func must be a lambda returning an Entrez handle or ez.read result.
    """
    for attempt in range(1, retries + 1):
        try:
            return request_func()
        except Exception as e:
            msg = str(e)

            # NCBI backend hiccups usually trigger RuntimeError with this text:
            transient = (
                "Search Backend failed" in msg or
                "HTTP Error" in msg or
                "temporarily unavailable" in msg or
                isinstance(e, RuntimeError)
            )

            if transient:
                print(f"[Entrez] transient error (attempt {attempt}): {e}")
                if attempt == retries:
                    raise  # re-raise on final failure
                time.sleep(delay)
            else:
                # unexpected error → do NOT retry
                raise


def get_tax_ids(ids: set[str]) -> set[str]:
    '''Gets list of tax ids to use for querying data'''
    tax_ids: set[str] = set()

    for i in ids:
        resp: dict[str, str] = ez.read(
            (ez.esearch(db="taxonomy", term=f'txid{i}[ORGN]', retmax=2000)))

        tax_ids.update(resp["IdList"])

    tax_ids -= ids
    return tax_ids


def create_csv(csv_name: Path, data: list[TaxInfo]):
    '''
    creates csv file per gene that contains infomration in the 
    TaxInfo class
    '''
    headers = ["tax_id", "scientific_name", "accession", "seq_name", "seq_len"]

    with open(csv_name, 'w', newline='') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=headers)

        writer.writeheader()

        for row in data:
            if row.accession is None:
                continue
            writer.writerow({
                'tax_id':  row.tax_id,
                'scientific_name': row.scientific_name,
                'accession': row.accession,
                'seq_name': row.seq_name,
                'seq_len': row.seq_len
            })


def main():
    '''This is the main function'''
    # set working directory
    wd = Path(__file__).resolve().parents[1]

    # load env variables for email and api key
    load_dotenv(wd / '.env')
    ez.email = os.getenv("ncbi_email")
    ez.api_key = os.getenv("ncbi_key")

    # parse command line arguments
    parser = argparse.ArgumentParser(
        description="A simple program to take in a json config and get data from NBCI databases for a csv files")

    parser.add_argument(
        '-f',
        '--file',
        default=wd / 'config.json',
        help='name of file for the program to injest'
    )

    args = parser.parse_args()

    tables_path = create_data_path(wd, "tables")

    # read in config file
    if not os.path.exists(wd / args.file):
        print('ERROR: Config file not found.')
        return

    with open(wd / args.file, 'r') as config_file:
        data = json.load(config_file)

    genera_ids: set[str] = set(data["genera_ids"])
    gene_json: list[dict[str, str]] = data["genes"]

    genes: list[Gene] = []

    for gene in gene_json:
        name, term, prefix = gene.values()
        genes.append(Gene(name, term, prefix))

    # gets tax ids
    txids: list[str] = list(get_tax_ids(genera_ids))

    for gene_config in genes:
        tax_data: list[TaxInfo] = []

        csv_name: Path = tables_path / f'{gene_config.file_prefix}_table.csv'

        with tqdm(total=len(txids), desc=gene_config.name) as progress:
            for tax_id in txids:
                tax = TaxInfo(tax_id)
                tax.get_sci_name()
                tax.nuc_processing(gene_config.term)
                tax_data.append(tax)

                progress.update(1)

        create_csv(csv_name, tax_data)


if __name__ == "__main__":
    main()
