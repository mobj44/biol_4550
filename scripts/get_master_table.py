from pathlib import Path
from utils import *
import pandas as pd


def clean_sci_name(name: str) -> str:
    return str(name).lower().replace(" ", "_")


def main():
    wd = Path(__file__).resolve().parents[1]
    tables_path = create_data_path(wd, 'tables')

    tables = read_table_names(wd)
    if not tables:
        print("No table files found.")
        return

    merged = None
    sci_names_frame = []

    for file in tables:
        if not file.endswith('_table.csv'):
            continue

        new_column_name = file.split('_')[0]
        df = pd.read_csv(tables_path / file)

        accessions_to_merge = df[["tax_id", "accession"]].rename(
            columns={"accession": new_column_name})

        if merged is None:
            merged = accessions_to_merge
        else:
            merged = pd.merge(
                merged,
                accessions_to_merge,
                on="tax_id",
                how="outer"
            )

        sci_names_frame.append(df[["tax_id", "scientific_name"]])

    sci_names_frame = (
        pd.concat(sci_names_frame, ignore_index=True)
        .drop_duplicates(subset="tax_id")
    )
    # add scientific names
    final = pd.merge(
        merged,
        sci_names_frame,
        on="tax_id",
        how="left"
    )
    # reorder columns
    accession_cols = [col for col in final.columns
                      if col not in ("tax_id", "scientific_name")]
    new_order = ["tax_id", "scientific_name"] + accession_cols
    final = final[new_order]

    out_path = tables_path / 'master_accessions.csv'
    final.to_csv(out_path, index=False)


if __name__ == '__main__':
    main()
