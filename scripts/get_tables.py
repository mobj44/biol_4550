import argparse
import csv
import http.client
import json
import os
import socket
import time
from dataclasses import dataclass
from pathlib import Path
from urllib.error import HTTPError, URLError

from Bio import Entrez as ez
from dotenv import load_dotenv
from tqdm import tqdm

from utils import *


# Keep the script conservative even when an API key is configured.
# Without a key, NCBI's public guidance is much lower, so we stay polite.
DEFAULT_REQUEST_DELAY_SECONDS = 0.12
DEFAULT_RETRIES = 6
DEFAULT_BACKOFF_BASE_SECONDS = 1.0


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
    seq_len: int | None = None

    def get_sci_name(self):
        """Set the scientific name for this taxon."""
        tax_sum = entrez_retry(
            lambda: ez.read(
                ez.esummary(
                    db="taxonomy",
                    id=self.tax_id,
                )
            )[0]
        )
        self.scientific_name = tax_sum["ScientificName"]

    def nuc_processing(self, gene: str):
        """Set accession, sequence title, and sequence length for this taxon."""
        nuc_search = entrez_retry(
            lambda: ez.read(
                ez.esearch(
                    db="nuccore",
                    term=f"txid{self.tax_id}[ORGN] AND {gene}",
                    retmax=1,
                )
            )
        )

        nuc_id = nuc_search["IdList"][0] if nuc_search["IdList"] else None
        if nuc_id is None:
            return

        nuc_sum = entrez_retry(
            lambda: ez.read(
                ez.esummary(
                    db="nuccore",
                    id=nuc_id,
                )
            )[0]
        )

        self.accession = nuc_sum.get("AccessionVersion")
        self.seq_name = nuc_sum.get("Title")

        length = nuc_sum.get("Length")
        self.seq_len = int(length) if length is not None else None


TRANSIENT_ERROR_SUBSTRINGS = (
    "search backend failed",
    "temporarily unavailable",
    "remote end closed connection",
    "connection reset",
    "connection aborted",
    "timed out",
    "bad gateway",
    "service unavailable",
    "gateway timeout",
    "502",
    "503",
    "504",
    "429",
    "incomplete read",
)


TRANSIENT_EXCEPTIONS = (
    RuntimeError,
    http.client.RemoteDisconnected,
    http.client.IncompleteRead,
    URLError,
    HTTPError,
    TimeoutError,
    socket.timeout,
    ConnectionResetError,
    ConnectionAbortedError,
    BrokenPipeError,
    OSError,
)


def _is_transient_error(exc: Exception) -> bool:
    """Best-effort classification of retryable Entrez/network failures."""
    msg = str(exc).lower()

    if isinstance(exc, HTTPError):
        return exc.code in {429, 500, 502, 503, 504}

    if isinstance(exc, TRANSIENT_EXCEPTIONS):
        return True

    return any(fragment in msg for fragment in TRANSIENT_ERROR_SUBSTRINGS)


def entrez_retry(request_func, retries: int = DEFAULT_RETRIES, base_delay: float = DEFAULT_BACKOFF_BASE_SECONDS):
    """
    Retry wrapper for Entrez calls.

    Uses exponential backoff and handles both NCBI backend hiccups and
    transient transport-level disconnects.
    """
    for attempt in range(1, retries + 1):
        try:
            return request_func()
        except Exception as exc:
            if not _is_transient_error(exc):
                raise

            print(
                f"[Entrez] transient error (attempt {attempt}/{retries}): {exc}")

            if attempt == retries:
                raise

            sleep_for = base_delay * (2 ** (attempt - 1))
            time.sleep(sleep_for)


def polite_pause(delay_seconds: float = DEFAULT_REQUEST_DELAY_SECONDS):
    """Small delay between requests to avoid hammering Entrez."""
    time.sleep(delay_seconds)


def get_tax_ids(ids: set[str]) -> set[str]:
    """Get descendant tax IDs to query from a set of higher-level taxonomy IDs."""
    tax_ids: set[str] = set()

    for taxon_id in ids:
        resp: dict[str, list[str]] = entrez_retry(
            lambda taxon_id=taxon_id: ez.read(
                ez.esearch(
                    db="taxonomy",
                    term=f"txid{taxon_id}[ORGN]",
                    retmax=2000,
                )
            )
        )
        tax_ids.update(resp["IdList"])
        polite_pause()

    tax_ids -= ids
    return tax_ids


def create_csv(csv_name: Path, data: list[TaxInfo]):
    """Create a CSV file for a gene using the data stored in TaxInfo objects."""
    headers = ["tax_id", "scientific_name", "accession", "seq_name", "seq_len"]

    with open(csv_name, "w", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=headers)
        writer.writeheader()

        for row in data:
            if row.accession is None:
                continue

            writer.writerow(
                {
                    "tax_id": row.tax_id,
                    "scientific_name": row.scientific_name,
                    "accession": row.accession,
                    "seq_name": row.seq_name,
                    "seq_len": row.seq_len,
                }
            )


def main():
    """Main entrypoint."""
    wd = Path(__file__).resolve().parents[1]

    load_dotenv(wd / ".env")
    ez.email = os.getenv("ncbi_email")
    ez.api_key = os.getenv("ncbi_key")

    parser = argparse.ArgumentParser(
        description="Take in a JSON config and query NCBI databases to build CSV tables."
    )
    parser.add_argument(
        "-f",
        "--file",
        default=wd / "config.json",
        help="name of file for the program to ingest",
    )
    args = parser.parse_args()

    tables_path = create_data_path(wd, "tables")

    config_path = Path(args.file)
    if not config_path.is_absolute():
        config_path = wd / config_path

    if not config_path.exists():
        print("ERROR: Config file not found.")
        return

    with open(config_path, "r") as config_file:
        data = json.load(config_file)

    genera_ids: set[str] = set(data["genera_ids"])
    gene_json: list[dict[str, str]] = data["genes"]

    genes: list[Gene] = []
    for gene in gene_json:
        name, term, prefix = gene.values()
        genes.append(Gene(name, term, prefix))

    txids: list[str] = list(get_tax_ids(genera_ids))

    for gene_config in genes:
        tax_data: list[TaxInfo] = []
        csv_name: Path = tables_path / f"{gene_config.file_prefix}_table.csv"

        with tqdm(total=len(txids), desc=gene_config.name) as progress:
            for tax_id in txids:
                try:
                    tax = TaxInfo(tax_id)
                    tax.get_sci_name()
                    polite_pause()

                    tax.nuc_processing(gene_config.term)
                    polite_pause()

                    tax_data.append(tax)
                except Exception as exc:
                    print(f"[WARN] tax_id {tax_id} failed: {exc}")
                finally:
                    progress.update(1)

        create_csv(csv_name, tax_data)


if __name__ == "__main__":
    main()
