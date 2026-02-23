#!/usr/bin/env python3
"""
Command-line script to run the PriorCons pipeline.

Usage example:
python intgrate_consensus.py \
  --input /path/to/alignment.aln \
  --ref REF_ID_IN_ALIGNMENT \
  --prior /path/to/priors.parquet \
  --output_dir /path/to/output_dir

This script performs:
 1) loads alignment (expects a function load_alignment in utils.py)
 2) extracts reference coordinates via extract_ref_positions (from utils.py)
 3) creates sliding windows (sliding_windows from utils.py)
 4) scores windows and decides QC (create_windows_df)
 5) builds consensus from windows (create_consensus)
 6) re-inserts mapping consensus insertions removed by initial filtering (add_mapping_insertions)
 7) writes results: windows_trace.csv, final fasta and qc.json
"""
import argparse
import logging
import sys
from pathlib import Path
import numpy as np
import pandas as pd

# Import helper functions from the package modules
from .utils_integrate_consensus import (
    create_windows_df,
    create_consensus,
    add_mapping_insertions,
    qc_process,
    logger as um_logger
)

from .utils import load_alignment, extract_ref_positions, sliding_windows


# configure top-level logger
logger = logging.getLogger("intgrate_consensus")
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setFormatter(um_logger.handlers[0].formatter)
logger.addHandler(ch)


def parse_args(argv=None):
    p = argparse.ArgumentParser(...)
    p.add_argument("--input", required=True, type=Path, help="Path to input alignment file (.aln)")
    p.add_argument("--ref", required=True, type=str, help="Reference sequence ID present in the alignment file")
    p.add_argument("--prior", required=True, type=Path, help="Path to prior parquet file")
    p.add_argument("--output_dir", required=True, type=Path, help="Output directory to write results")
    return p.parse_args(argv)



def main(argv: list[str] = None) -> int:
    argv = argv if argv is not None else sys.argv[1:]
    args = parse_args(argv)
    input_path: Path = args.input
    ref_id: str = args.ref
    prior_path: Path = args.prior
    output_dir: Path = args.output_dir
    
    try:
        # Basic validations
        if not input_path.exists():
            logger.error("Input alignment file not found: %s", input_path)
            sys.exit(2)
        if not prior_path.exists():
            logger.error("Prior file not found: %s", prior_path)
            sys.exit(2)

        output_dir.mkdir(parents=True, exist_ok=True)
        logger.info("Output directory: %s", output_dir)

        # Read prior table
        logger.info("Loading prior table from %s", prior_path)
        prior_table = pd.read_parquet(prior_path)

        # Load alignment
        logger.info("Loading alignment from %s", input_path)
        ids, seqs = load_alignment(str(input_path))
        original_seqs = {ids[index]: seqs[index] for index in range(len(ids))}

        # Obtain the ref, mapping consensus and abacas consensus sequence ids
        ref_key,mapping_key,abacas_key =  ids

        logger.info(f"Sequences loaded ids are:\n- REF: {ids[0]}\n- Mapping consensus: {ids[1]}\n- ABACAS consensus: {ids[2]}")
        
        if ref_id != ref_key:
            logger.error("Reference id '%s' not found in alignment IDs", ref_id)
            sys.exit(2)

        # Filter alignment by reference coordinates
        logger.info("Extracting reference-coordinates sequences using ref id '%s'", ref_id)
        seqs_filtered, ref_seq = extract_ref_positions(ids=ids, seqs=seqs, ref_id=ref_id)
        filtered_seqs = {ids[index]: seqs_filtered[index] for index in range(len(ids))}

        # Build sliding windows (example: win=100 ovlp=50 as in your notebook)
        logger.info("Creating sliding windows")
        windows = sliding_windows(len(filtered_seqs[ref_id]), win_size=100, overlap=50)

        # Create windows df based on priors and scores
        logger.info("Scoring windows and creating windows dataframe")
        window_df = create_windows_df(windows, filtered_seqs, prior_table=prior_table,
                                      mapp_id=mapping_key,abacas_id=abacas_key)

        window_csv = output_dir / "windows_trace.csv"
        window_df.to_csv(window_csv, index=False)
        logger.info("Windows trace written to %s", window_csv)

        # Create integrated consensus
        logger.info("Creating integrated consensus from ABACAS and Mapping sequences")
        
        cons = create_consensus(abacas_seq=filtered_seqs[abacas_key],
                                mapp_seq=filtered_seqs[mapping_key],
                                window_df=window_df)

        # Insert mapping consensus insertions found in the original aligned mapping cons vs REF
        logger.info("Detecting and re-inserting mapping consensus insertions removed by filtering")
        final_with_ins, insertions = add_mapping_insertions(mapp_aln=original_seqs[mapping_key],
                                                         ref_aln=original_seqs[ref_id],
                                                         ref_filt = filtered_seqs[ref_id],
                                                         final_seq=cons)

        # Write final fasta (remove '-' gaps if any)
        fin_arr = np.array(final_with_ins.split("-"))
        fin_text = "".join(fin_arr[fin_arr != ""])
        base_name = input_path.stem
        fasta_path = output_dir / f"integrated_consensus.fasta"
        header = f">{base_name}\n"
        with open(fasta_path, "w") as fh:
            fh.write(header)
            fh.write(fin_text + "\n")
        logger.info("Final integrated consensus FASTA written to %s", fasta_path)
        # Perform QC and write qc.json
        qc_path = output_dir / "qc.json"
        qc_process(filtered_seqs=filtered_seqs, integrated_seq=cons, insertions=insertions, write=str(qc_path),ref_name=ref_id,mapp_id=mapping_key)

        logger.info("QC written to %s", qc_path)
        logger.info("PriorCons tool finished successfully")

    except Exception as e:
        logger.exception("PriorCons tool failed: %s", e)
        sys.exit(1)


if __name__ == "__main__":
    main()
