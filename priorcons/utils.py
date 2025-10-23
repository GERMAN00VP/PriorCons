import logging
import math
from Bio import AlignIO
import numpy as np
import pandas as pd

# ------------------------------------------------------
# Logging configuration
# ------------------------------------------------------
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# By default log to console
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("[%(levelname)s] %(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


# ------------------------------------------------------
# Alignment and window handling
# ------------------------------------------------------
def load_alignment(fasta_file):
    """Load a multiple sequence alignment in FASTA format."""
    try:
        alignment = AlignIO.read(fasta_file, "fasta")
    except FileNotFoundError:
        logger.error(f"File not found: {fasta_file}")
        raise
    except Exception as e:
        logger.error(f"Error reading alignment: {e}")
        raise

    ids = [record.id for record in alignment]
    seqs = [str(record.seq).upper() for record in alignment]

    if not ids:
        raise ValueError("Alignment file is empty.")

    seq_lengths = {len(s) for s in seqs}
    if len(seq_lengths) > 1:
        raise ValueError("Sequences in the alignment do not have the same length.")

    logger.info(f"Loaded alignment with {len(ids)} sequences, length {len(seqs[0])}")
    return ids, seqs


def extract_ref_positions(ids, seqs, ref_id):
    """Keep only positions where reference has no gap."""
    if ref_id not in ids:
        raise ValueError(f"Reference ID '{ref_id}' not found in alignment.")

    ref_index = ids.index(ref_id)
    ref_seq = seqs[ref_index]
    keep_positions = [i for i, base in enumerate(ref_seq) if base != "-"]

    filtered = ["".join(seq[i] for i in keep_positions) for seq in seqs]
    logger.info(f"Filtered reference positions: kept {len(keep_positions)} bases")
    return filtered, filtered[ref_index]


def sliding_windows(seq_len, win_size=100, overlap=10):
    """Generate start/end indices for sliding windows."""
    step = win_size - overlap
    if step <= 0:
        raise ValueError("Overlap must be smaller than window size.")
    if win_size > seq_len:
        raise ValueError("Window size larger than sequence length.")

    return [(start, start + win_size) for start in range(0, seq_len - win_size + 1, step)]


# ------------------------------------------------------
# Profiles and scoring
# ------------------------------------------------------
def compute_window_distribution(seqs, window, alpha=1):
    """
    Compute base probability distribution for each position in a window.
    """
    start, end = window
    dist = []
    for pos in range(start, end):
        counts = {b: alpha for b in "ACGT"}
        total = 4 * alpha
        for seq in seqs:
            base = seq[pos]
            if base in "ACGT":
                counts[base] += 1
                total += 1
        probs = {b: counts[b] / total for b in "ACGT"}
        dist.append(probs)
    return dist


def build_window_profiles(ids, seqs, ref_id, win_size=100, overlap=10, alpha=1):
    """
    Build probability profiles for all sliding windows.
    """
    seqs_filtered, ref_seq = extract_ref_positions(ids, seqs, ref_id)
    seq_len = len(ref_seq)
    windows = sliding_windows(seq_len, win_size, overlap)

    profiles = []
    for start, end in windows:
        dist = compute_window_distribution(seqs_filtered, (start, end), alpha)
        profiles.append({
            "start": start,
            "end": end,
            "distribution": dist
        })
    logger.info(f"Built {len(profiles)} window profiles.")
    return profiles, seqs_filtered, ref_seq


def score_window(seq, window_profile):
    """
    Compute log-likelihood of one sequence in one window.
    """
    start, end = window_profile["start"], window_profile["end"]
    dist = window_profile["distribution"]

    loglik = 0.0
    valid_positions = 0
    for i, probs in enumerate(dist):
        base = seq[start + i]
        if base in "ACGT":
            prob = probs.get(base, 0)
            if prob <= 0:
                return {"start": start, "end": end, "loglik": float("-inf"), "nLL": float("inf")}
            loglik += math.log(prob)
            valid_positions += 1
    if valid_positions > 0:
        nLL = -loglik / valid_positions
    else:
        nLL = float("nan")

    return {"start": start, "end": end, "loglik": loglik, "nLL": nLL}


def score_all_windows(seq, profiles):
    """
    Score one sequence across all windows.
    Returns list of dicts.
    """
    return [score_window(seq, prof) for prof in profiles]


# ------------------------------------------------------
# Prior table construction
# ------------------------------------------------------
def build_priors(seqs, profiles, percentiles=[95, 99]):
    """
    Compute empirical nLL distributions per window across all sequences.
    Returns a DataFrame with percentiles and window info.
    """
    records = []
    for prof in profiles:
        start, end = prof["start"], prof["end"]
        nLLs = []
        for seq in seqs:
            result = score_window(seq, prof)
            if not math.isnan(result["nLL"]):
                nLLs.append(result["nLL"])
        if not nLLs:
            continue

        rec = {"start": start, "end": end}
        for p in percentiles:
            rec[f"nLL_p{p}"] = np.percentile(nLLs, p)
        rec["profile"] = prof["distribution"]  # can be large but needed
        records.append(rec)

    df = pd.DataFrame(records)
    logger.info(f"Built priors table with {len(df)} windows.")
    return df


def save_priors(df, out_file):
    """
    Save priors DataFrame to compressed Parquet format.
    """
    try:
        df.to_parquet(out_file, compression="zstd", index=False)
        logger.info(f"Priors saved to {out_file} (compressed Parquet).")
    except Exception as e:
        logger.error(f"Error saving priors: {e}")
        raise
