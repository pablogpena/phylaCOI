# -*- coding: utf-8 -*-
# +
# -*- coding: utf-8 -*-
# -

from pathlib import Path
from otu_utils import process_folder
import sys
sys.path.append('/workspace/PhylaCOI')
from src.otu_utils import process_folder

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate OTUs per folder by cleaning FASTA files and performing clustering with VSEARCH."
    )
    parser.add_argument(
        "-r", "--root",
        type=str,
        required=True,
        help="Root directory containing subfolders with an 'output' subdirectory. Created in the previus step (get_abundance.py)"
    )
    parser.add_argument(
        "-i", "--identity",
        type=float,
        default=0.97,
        help="Sequence identity threshold for VSEARCH clustering (default: 0.97)."
    )
    parser.add_argument(
        "--no-vsearch",
        action="store_true",
        help="If set, skips the VSEARCH execution (only performs FASTA cleaning)."
    )

    args = parser.parse_args()
    root_dir = Path(args.root)
    run_vsearch_flag = not args.no_vsearch

    for folder in root_dir.iterdir():
        if folder.is_dir() and (folder / "output").is_dir():
            print(f"Processing: {folder.name}")
            process_folder(folder, args.identity, run_vsearch_flag)
