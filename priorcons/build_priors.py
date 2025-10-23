import argparse
import sys
from typing import List
from .utils import (
    save_priors,
    load_alignment,
    build_window_profiles,
    build_priors
)

def main(argv: List[str] = None) -> int:
    parser = argparse.ArgumentParser(description="Build empirical priors from alignment")
    parser.add_argument("-i", "--input", required=True, help="Input alignment FASTA file")
    parser.add_argument("-r", "--ref", required=True, help="Reference sequence ID")
    parser.add_argument("-o", "--output", required=True, help="Output file (.parquet)")
    parser.add_argument("--win", type=int, default=100, help="Window size (default=100)")
    parser.add_argument("--overlap", type=int, default=50, help="Window overlap (default=50)")
    argv = argv if argv is not None else sys.argv[1:]
    args = parser.parse_args(argv)
    
    # Load alignment
    ids, seqs = load_alignment(args.input)

    # Build window profiles
    profiles, seqs_filtered, ref_seq = build_window_profiles(
        ids, seqs, args.ref, args.win, args.overlap
    )

    # Compute priors
    df = build_priors(seqs_filtered, profiles)

    # Save to Parquet
    save_priors(df, args.output)

    return 0


if __name__ == "__main__":
    main()
