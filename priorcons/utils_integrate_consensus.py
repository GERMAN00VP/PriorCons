#!/usr/bin/env python3
"""
Utility functions for the consensus integration workflow.

This module implements helper functions that:
 - compute fragmentation/missingness,
 - count substitutions,
 - summarize insertions,
 - perform QC reporting,
 - detect and insert mapping consensus insertions into a reference-coordinate consensus,
 - create a consensus from window decisions,
 - create a windows DataFrame using priors and score_window (score_window is expected
   to live in utils.py and to accept full sequences and a window_profile with start/end).
 
 NOTE: This module imports helper routines from an existing utils.py (load_alignment,
       extract_ref_positions, sliding_windows, score_window, etc). Keep that file
       available in the same environment.
"""
from typing import List, Tuple, Dict, Any
from collections import defaultdict
import logging
import json
import numpy as np
import pandas as pd

from .utils import score_window


# configure module-level logger
logger = logging.getLogger("utils_integrate_consensus")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
formatter = logging.Formatter("[%(asctime)s] %(levelname)s - %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)


def how_fragmented(seq: str) -> int:
    """
    Return number of non-empty fragments when splitting seq by '-' gaps.
    Equivalent to counting contiguous non-gap blocks.

    Args:
        seq: aligned sequence (may contain '-')

    Returns:
        int: number of fragments (non-empty segments)
    """
    if seq is None:
        logger.debug("how_fragmented received None")
        return 0
    try:
        return int((np.array(seq.split("-")) != "").sum())
    except Exception as e:
        logger.error("Error in how_fragmented: %s", e)
        raise


def count_missing(seq: str) -> int:
    """
    Count missing positions in a sequence. Missing are N, n or '-'.

    Args:
        seq: sequence string

    Returns:
        int: number of missing positions
    """
    if seq is None:
        logger.debug("count_missing received None")
        return 0
    try:
        # split on N/n/- and count empty segments -> trick: we want number of missing characters
        # simpler: count characters that are N/n or '-'
        return int(sum(1 for c in seq if c in ("N", "n", "-")))
    except Exception as e:
        logger.error("Error in count_missing: %s", e)
        raise


def count_substitutions(ref: str, seq: str) -> int:
    """
    Count substitutions of seq relative to ref (ignoring gaps and Ns in seq).

    Args:
        ref: reference sequence (ungapped, same coords as seq or aligned but used in same coordinate space)
        seq: sequence to compare (same length as ref)

    Returns:
        int: number of substitutions (positions where seq has a different base not in [-,n,N])
    """
    if ref is None or seq is None:
        raise ValueError("ref and seq must be provided")
    if len(ref) != len(seq):
        raise ValueError("ref and seq must have the same length for substitution counting")

    muts = 0
    for pos, base in enumerate(ref):
        sbase = seq[pos]
        if (base != sbase) and (sbase not in ["-", "n", "N"]):
            muts += 1
    return muts


def process_insertions(insertions: List[Tuple[int, str]]) -> Tuple[int, int]:
    """
    Return number of insertions and total length of inserted bases.

    Args:
        insertions: list of (ref_pos, insertion_string)

    Returns:
        (n_insertions, total_length)
    """
    length_insertions = 0
    n_insertions = len(insertions)
    for insertion in insertions:
        length_insertions += len(insertion[1])
    return n_insertions, length_insertions


def qc_process(filtered_seqs: Dict[str, str],
               ref_name:str,
               mapp_id:str,
               integrated_seq: str,
               insertions: List[Tuple[int, str]],
               write: Any = False) -> Dict[str, Any]:
    """
    Compute QC metrics for the integrated consensus.

    Args:
        filtered_seqs: dict of filtered sequences )
        ref_name: str, reference id.
        integrated_seq: integrated consensus sequence (reference-coordinate, no gaps)
        insertions: list of (ref_pos, insertion_string) that were added to integrated sequence
        write: False or a filename. If filename, QC dict is written as JSON to that path.

    Returns:
        dict with QC metrics (or writes it to a file if write is a path)
    """
    try:
        if ref_name not in filtered_seqs:
            raise KeyError(f"filtered_seqs must contain {ref_name}")

        if not integrated_seq:
            raise ValueError("integrated_seq is empty or None")

        ref_len = len(integrated_seq)
        coverage_inicial = round((len(integrated_seq) - count_missing(filtered_seqs[mapp_id])) / len(integrated_seq) * 100, 3)
        coverage_final = round((len(integrated_seq) - count_missing(integrated_seq)) / len(integrated_seq) * 100, 3)

        integrated_seq_substitutions = count_substitutions(ref=filtered_seqs[ref_name], seq=integrated_seq)
        mapp_substitutions = count_substitutions(filtered_seqs[ref_name], filtered_seqs[mapp_id])

        # expected substitutions scaled by coverage change
        exp_substitutions = round((coverage_final * mapp_substitutions) / coverage_inicial, 3) if coverage_inicial != 0 else float("nan")
        obs_vs_exp_subst = round(integrated_seq_substitutions - exp_substitutions, 3)

        n_insertions, length_insertions = process_insertions(insertions)

        res_dict = {
            "MAPPING_CONSENSUS_COVERAGE": coverage_inicial,
            "FINAL_COVERAGE": coverage_final,
            "MAPPING_CONSENSUS_SUBSTITUTIONS": mapp_substitutions,
            "FINAL_SUBSTITUTIONS": integrated_seq_substitutions,
            "EXPECTED_SUBSTITUTIONS": exp_substitutions,
            "OBS-EXP_SUBSTITUTIONS": obs_vs_exp_subst,
            "N_INSERTIONS": n_insertions,
            "TOTAL_INSERTIONS_LENGTH": length_insertions,
            "INSERTIONS": insertions
        }

        if write:
            with open(write, "w") as f:
                json.dump(res_dict, f, indent=2)
            logger.info("QC results written to %s", write)
        else:
            return res_dict

    except Exception as e:
        logger.exception("qc_process failed: %s", e)
        raise


def add_mapping_insertions(mapp_aln: str, ref_aln: str, final_seq: str) -> Tuple[str, List[Tuple[int, str]]]:
    """
    Detect insertions in `mapp_aln` relative to `ref_aln` (both aligned strings, with '-' gaps),
    and insert those insertions into `final_seq` (which must be in reference coordinates, i.e. the
    same length as the number of non-gap positions in ref_aln).

    Returns:
        final_with_insertions, insertions_list

    Note:
        ref_pos in insertions is 0-based index of reference base AFTER which insertion should be placed.
        -1 means insert before the first base.
    """
    try:
        # Basic checks
        if len(mapp_aln) != len(ref_aln):
            raise ValueError("mapp_aln and ref_aln must have the same length (aligned with gaps).")

        # Count non-gap positions in reference (this must match len(final_seq))
        ref_non_gap_count = sum(1 for c in ref_aln if c != '-')
        if len(final_seq) != ref_non_gap_count:
            raise ValueError(
                f"final_seq length ({len(final_seq)}) must equal number of non-gap positions in ref_aln ({ref_non_gap_count})."
            )

        # Collect insertions from mapp relative to ref
        insertions: List[Tuple[int, str]] = []
        ref_pos = -1   # index in ungapped reference of the last reference base seen (-1 if none yet)
        pos = 0
        L = len(ref_aln)

        while pos < L:
            # If this column is a reference base, advance the reference coordinate
            if ref_aln[pos] != '-':
                ref_pos += 1

            # If we find a run where ref has '-' but mapp has bases -> insertion in mapp
            if ref_aln[pos] == '-' and mapp_aln[pos] != '-':
                ins_chars: List[str] = []
                # collect contiguous insertion characters while ref_aln is '-' and mapp_aln is not '-'
                while pos < L and ref_aln[pos] == '-' and mapp_aln[pos] != '-':
                    ins_chars.append(mapp_aln[pos])
                    pos += 1
                insertion = ''.join(ins_chars)
                insertions.append((ref_pos, insertion))
                # continue the outer loop from current pos
                continue

            # otherwise advance one column
            pos += 1

        # Group insertions by reference coordinate for efficient insertion
        by_ref_pos: Dict[int, List[str]] = defaultdict(list)
        for rp, seq in insertions:
            by_ref_pos[rp].append(seq)  # keep order found (first-to-last)

        # Build final sequence with insertions
        out_parts: List[str] = []
        # Insert any sequences that belong before the first ref base (ref_pos == -1)
        if -1 in by_ref_pos:
            out_parts.extend(by_ref_pos[-1])

        # Walk through the final_seq bases and insert after each base any insertions anchored there
        for i, base in enumerate(final_seq):
            out_parts.append(base)
            if i in by_ref_pos:
                out_parts.extend(by_ref_pos[i])

        final_with_insertions = ''.join(out_parts)
        return final_with_insertions, insertions

    except Exception:
        logger.exception("add_mapping_insertions failed")
        raise


def create_consensus(abacas_seq: str, mapp_seq: str, window_df: pd.DataFrame) -> str:
    """
    Create consensus sequence based on window selections.

    The function expects two sequences of the same length (reference coordinate).
    For each window, if WINDOW_QC_PASSED is True the base from abacas_seq is used for that window,
    otherwise the base from mapp_seq is used. Windows may overlap; ABACAS takes priority where windows
    indicate ABACAS support since we fill windows in order and fallback to ABACAS when positions are
    not filled.

    Args:
        abacas_seq: ABACAS sequence (reference-coordinate)
        mapp_seq: mapping consensus sequence (reference-coordinate)
        window_df: DataFrame with columns ['start', 'end', 'WINDOW_QC_PASSED'] (or similar)

    Returns:
        consensus sequence (string)
    """
    try:
        if len(abacas_seq) != len(mapp_seq):
            raise ValueError("ABACAS and mapp sequences must have the same length")

        seq_length = len(abacas_seq)
        consensus = [''] * seq_length  # Initialize empty consensus

        # Process each window
        for _, row in window_df.iterrows():
            start = int(row['start'])
            end = int(row['end'])
            source = bool(row.get('WINDOW_QC_PASSED', False))

            # Validate window coordinates
            if start < 0 or end > seq_length or start >= end:
                logger.error("Invalid window coordinates: %s-%s", start, end)
                raise ValueError(f"Invalid window coordinates: {start}-{end}")

            # Select sequence based on source
            selected_seq = abacas_seq if source else mapp_seq

            # Fill consensus for this window
            for pos in range(start, end):
                # ABACAS has priority in overlapped regions as requested (we write abacas if source True)
                consensus[pos] = selected_seq[pos]

        # Fill any remaining empty positions with ABACAS as fallback
        for i in range(seq_length):
            if not consensus[i]:
                consensus[i] = abacas_seq[i]

        return ''.join(consensus)

    except Exception:
        logger.exception("create_consensus failed")
        raise


def create_windows_df(windows: List[Tuple[int, int]],
                      filtered_seqs: Dict[str, str],
                      prior_table: pd.DataFrame,
                      mapp_id:str,
                      abacas_id:str) -> pd.DataFrame:
    """
    For each window compute QC statistics and compare window score against prior.
    Uses score_window (imported from utils.py) which must accept full sequence (not sliced)
    and a window_profile dict containing 'start','end','distribution'.

    Args:
        windows: list of (start,end) tuples (reference coordinates)
        filtered_seqs: dict with keys mapp_id and "ordered_REF ID" (both are reference-coordinate strings)
        prior_table: DataFrame with prior rows where .profile is an iterable distribution and .nLL_p95 exists

    Returns:
        pandas.DataFrame with computed window info and QC decision in WINDOW_QC_PASSED column
    """
    window_filtering = []

    for window in windows:
        start, end = window
        # slice sequences for diagnostics
        mapp = filtered_seqs[mapp_id][start:end]
        abacas_id = str(np.array([id for id in filtered_seqs.keys()])[np.array(["ordered" in id for id in filtered_seqs.keys()])][0])
        abacas = filtered_seqs[abacas_id][start:end]
        missing_mapp= count_missing(mapp)
        missing_abacas = count_missing(abacas)
        abacas_fragments = how_fragmented(abacas)

        window_dict = {
            "start": start,
            "end": end,
            "MISSING_mapp": missing_mapp,
            "MISSING_ABACAS": missing_abacas,
            "ABACAS_MORE_INFO": missing_mapp > missing_abacas,
            "ABACA_FRAGMENTS": f"{abacas_fragments}",
            "WINDOW_PRIOR_nLL_p95": "",
            "WINDOW_SCORE_nLL": "",
            "WINDOW_QC_PASSED": False
        }

        try:
            # Apply the same filter logic (as you wrote): require limited fragmentation and abacas provides more info
            if (abacas_fragments in [1, 2] and (missing_mapp > missing_abacas)) and missing_abacas < 50:
                prior = prior_table.loc[prior_table["start"] == start]
                if prior.empty:
                    logger.warning("No prior found for window start=%s; skipping scoring", start)
                else:
                    window_profile = {
                        "start": start,
                        "end": end,
                        "distribution": prior.profile.values[0]
                    }

                    # score_window expects the full sequence and uses window_profile["start"] internally
                    window_score = score_window(seq=filtered_seqs[abacas_id],
                                                window_profile=window_profile)['nLL']

                    window_dict["WINDOW_PRIOR_nLL_p95"] = prior.nLL_p95.values[0]
                    window_dict["WINDOW_SCORE_nLL"] = window_score
                    window_dict["WINDOW_QC_PASSED"] = prior.nLL_p95.values[0] > window_score

        except Exception:
            logger.exception("Error while scoring window %s-%s", start, end)
            # keep default window_dict values (QC False)

        window_filtering.append(window_dict)

    windows_df = pd.DataFrame(window_filtering)
    return windows_df
