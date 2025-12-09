#!/usr/bin/env python3
"""
Command-line script to run the PriorCons pipeline.

Usage example:
python intgrate_consensus.py \
  --input /path/to/alignment.aln \
  --ref REF_ID_IN_ALIGNMENT \
  --output_dir

"""
import argparse
import logging
import sys
from pathlib import Path

# Import helper functions from the package modules
from .utils_integrate_consensus import (

    add_gf,
    logger as um_logger
)

from .utils import load_alignment


# configure top-level logger
logger = logging.getLogger("COMPLETE_WITH_GF")
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setFormatter(um_logger.handlers[0].formatter)
logger.addHandler(ch)


def parse_args(argv=None):
    p = argparse.ArgumentParser(...)
    p.add_argument("--input", required=True, type=Path, help="Path to input alignment file (.aln)")
    p.add_argument("--ref", required=True, type=str, help="Reference sequence ID present in the alignment file")
    p.add_argument("--output_dir", required=True, type=Path, help="Output directory to write results")
    return p.parse_args(argv)



def main(argv: list[str] = None) -> int:
    argv = argv if argv is not None else sys.argv[1:]
    args = parse_args(argv)
    input_path: Path = args.input
    ref_id: str = args.ref
    output_dir: Path = args.output_dir
    
    try:
        # Basic validations
        if not input_path.exists():
            logger.error("Input alignment file not found: %s", input_path)
            sys.exit(2)


        output_dir.mkdir(parents=True, exist_ok=True)
        logger.info("Output directory: %s", output_dir)


        # Load alignment
        logger.info("Loading alignment from %s", input_path)
        ids, seqs = load_alignment(str(input_path))
        original_seqs = {ids[index]: seqs[index] for index in range(len(ids))}

        # Obtain the ref, mapping consensus and abacas consensus sequence ids
        ref_key,mapping_key,gf_key =  ids

        logger.info(f"Sequences loaded ids are:\n- REF: {ref_key}\n- Mapping consensus: {mapping_key}\n- GF mapping consensus: {gf_key}")
        
        if ref_id != ref_key:
            logger.error("Reference id '%s' not found in alignment IDs", ref_id)
            sys.exit(2)


        # Create integrated consensus
        logger.info("Adding info from GF to Mapping consensus sequences")
        
        cons_gf = add_gf(ref_seq=original_seqs[ref_key],
                      mapp_seq=original_seqs[mapping_key],
                      gf_seq=original_seqs[gf_key]
        )
                      

        base_name = input_path.stem
        fasta_path = output_dir / "PRIORCONS_GF_INTEGRATED.fasta"
        header = f">{base_name}_GF\n"
        with open(fasta_path, "w") as fh:
            fh.write(header)
            fh.write(cons_gf + "\n")
        logger.info("Final integrated consensus FASTA written to %s", fasta_path)
  
        logger.info("PriorCons tool finished successfully")

    except Exception as e:
        logger.exception("PriorCons tool failed: %s", e)
        sys.exit(1)


if __name__ == "__main__":
    main()
