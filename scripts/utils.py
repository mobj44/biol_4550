from pathlib import Path
import os


def create_data_path(wd: Path, folder_name: str) -> Path:
    '''set/create folder path data will go in'''
    folder_path: Path = wd / 'data' / folder_name

    try:
        folder_path.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        print(f"Error creating folder: {e}")
    return folder_path


def read_table_names(wd: Path) -> list[str] | None:
    table_dir = wd / 'data' / 'tables'

    # check if tables folder exists
    if not os.path.exists(table_dir) or not table_dir.is_dir():
        print('ERROR: Folder not found')
        return

    # get file names
    return [file for file in os.listdir(
        table_dir) if os.path.isfile(table_dir / file)]
