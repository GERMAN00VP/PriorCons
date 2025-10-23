# PriorCons 

This repository provides tools to:

1. **Generate Integrated Consensus (`integrate_consensus.py`)**  
   Produces a high-quality viral consensus by strategically using **ABACAS** sequences to fill missing regions in the **mapping** consensus. It employs a sliding-window approach that verifies the evolutionary plausibility of ABACAS content against empirical priors before incorporation.

2. **Build Evolutionary Priors (`build_priors.py`)**  
   Constructs empirical prior distributions from large multiple-sequence alignments. These priors model expected genetic variation across genomic windows and provide likelihood thresholds for quality control during consensus integration.

3. **Access Supporting Utilities (`utils` scripts)**  
   Provides modular helper functions for alignment processing, window scoring, and consensus construction used by both main workflows.

---
### Installation

```bash
pip install priorcons
```
---

### CLI usage 
```bash
# Create priors
priorcons build-priors --input sequences.fasta --ref REF_ID --output priors.parquet

# Run consensus integration
priorcons integrate-consensus --input alignment.aln --ref REF_ID --prior priors.parquet --output_dir results

```

---

## 🚀 Main Script: integrate_consensus.py

This is the entrypoint of the tool. It creates a *integrated consensus sequence* by combining mapping consensus and ABACAS output, both aligned to a reference sequence, but only after performing quality control (QC) at the window level.  

### 🔑 Inputs

- `--input` → path to an alignment file (`.aln`) containing at least:  
  - **1º Reference sequence**  
  - **2º Mapping consensus sequence**  
  - **3º ABACAS consensus sequence**  

 **The sequences in the alignment file must be provided in the specified order, as they will be identified by their position.**


- `--ref` → ID of the reference sequence in the alignment.  

- `--prior` → path to a priors table (`.parquet`) generated with `build_priors.py`.  

- `--output_dir` → directory to save the results.  


---

### 🧪 Workflow

1. **Start with mapping consensus** as the baseline
2. **Identify missing/unreliable regions** in mapping consensus
3. **For each window**:
   - If mapping has coverage → keep mapping sequence
   - If mapping has missing data → evaluate ABACAS for that window:
     * Check fragmentation and quality
     * Verify evolutionary plausibility using priors [(nLL score)](README.md#-methodology-build_priorspy)
     * If ABACAS passes QC → use ABACAS to fill missing regions
4. **Construct final consensus** combining mapping baseline with validated ABACAS fills
5. **Restore mapping-specific insertions**
6. **QC reporting**: compute coverage, substitutions, and insertion metrics comparing the final integrated consensus to MAPPING.  

---

### 📦 Outputs

The script produces **three files** inside `--output_dir`:

1. **Integrated consensus FASTA**  
   - File: `<basename>-INTEGRATED.fasta`  
   - Contains the final consensus sequence after merging and reinserting insertions.  

2. **Window QC trace (CSV)**  
   - File: `windows_trace.csv`  
   - One row per window, recording:  
     - `start`, `end` → genomic coordinates.  
     - `MISSING_MAPPING`, `MISSING_ABACAS` → counts of missing bases.  
     - `ABACAS_MORE_INFO` → whether ABACAS has fewer missing bases than MAPPING.  
     - `ABACAS_FRAGMENTS` → fragmentation level of ABACAS in this window (keep: 0 < n fragments < 3 ).  
     - `WINDOW_PRIOR_nLL_p95` → threshold from priors.  
     - `WINDOW_SCORE_nLL` → score of ABACAS in this window.  
     - `WINDOW_QC_PASSED` → True/False decision.  

3. **Consensus QC summary (JSON)**  
   - File: `qc.json`  
   - Provides overall metrics comparing the MAPPING consensus and the integrated consensus:  
     - `MAPPING_COVERAGE` → % of genome covered in MAPPING.  
     - `FINAL_COVERAGE` → % of genome covered in integrated consensus.  
     - `MAPPING_SUBSTITUTIONS` → substitutions vs. reference in MAPPING.  
     - `FINAL_SUBSTITUTIONS` → substitutions vs. reference in integrated consensus.  
     - `EXPECTED_SUBSTITUTIONS` → expected number of substitutions, extrapolated from mapping.  
     - `OBS-EXP_SUBSTITUTIONS` → difference between observed and expected substitutions.  
     - `N_INSERTIONS` → number of insertions added back.  
     - `TOTAL_INSERTIONS_LENGTH` → total inserted length.  
     - `INSERTIONS` → list of insertions with their coordinates.  

---

### ▶️ Example run

```bash
python integrate_consensus.py \
  --input /path/to/250694-RSVWGS.aln \
  --ref RSV_BD \
  --prior /path/to/RSVBD_win100_ovlp50_priors.parquet \
  --output_dir results
```
----

This will generate:  

- `results/250694-RSVWGS-INTEGRATED.fasta`  
- `results/windows_trace.csv`  
- `results/qc.json`  

---

## 🛠 Script: build_priors.py

This script creates **empirical priors** (overlapped windows) from a large multiple sequence alignment.  
These priors are later used by `integrate_consensus.py` to evaluate windows.

### 🔑 Inputs

- `-i / --input` → aligned FASTA file with multiple sequences.  
- `-r / --ref` → ID of the reference sequence.  
- `-o / --output` → output file (`.parquet`).  
- `--win` → window size (default: 100).  
- `--overlap` → overlap size (default: 10).  

### ▶️ Example run

```bash
python build_priors.py \
  -i alignment.fasta \
  -r ReferenceID \
  -o priors.parquet \
  --win 100 \
  --overlap 10
```
----

### 📦 Output

A `.parquet` file with one row per window, containing:  

- `start`, `end` → window coordinates.  
- `nLL_p95`, `nLL_p99` → empirical thresholds.  
- `profile` → base probability distributions for each position in the window.  

## 🧮 Methodology (build_priors.py)

### 1. Probability distributions per position

For each window of size `W` bases (e.g., `W = 100`), and for each position `j` within that window, we compute the probability of observing each nucleotide:

![P_j(b)](https://latex.codecogs.com/png.latex?P_j(b)=\frac{c_j(b)+\alpha}{\sum_{x\in\{A,C,G,T\}}(c_j(x)+\alpha)})

Where:  
- ![c_j(b)](https://latex.codecogs.com/png.latex?c_j(b)) = number of sequences with base ![b](https://latex.codecogs.com/png.latex?b) at position ![j](https://latex.codecogs.com/png.latex?j).  
- ![\alpha](https://latex.codecogs.com/png.latex?\alpha) = pseudocount (Laplace smoothing, default ![\alpha=1](https://latex.codecogs.com/png.latex?\alpha=1)) to avoid zero probabilities.  
- Bases `N` are ignored in the counts.  

This gives a **per-position categorical distribution**.

---

### 2. Log-likelihood of a sequence in a window

Given a query sequence ![Q](https://latex.codecogs.com/png.latex?Q), we compute its probability under the window profile.  
For each valid (non-`N`) position ![j](https://latex.codecogs.com/png.latex?j) with observed base ![q_j](https://latex.codecogs.com/png.latex?q_j):

![logL](https://latex.codecogs.com/png.latex?\log%20L(Q\mid%20\text{window})=\sum_{j=1}^{W}\log%20P_j(q_j))

The **normalized negative log-likelihood (nLL)** is:

![nLL](https://latex.codecogs.com/png.latex?\text{nLL}(Q)=-\frac{1}{N_{\text{valid}}}\sum_{j=1}^{W}\log%20P_j(q_j))

Where:  
- ![N_valid](https://latex.codecogs.com/png.latex?N_{\text{valid}}) = number of positions in the window where ![Q](https://latex.codecogs.com/png.latex?Q) has a non-`N` base.  

Smaller nLL values indicate sequences more likely under the empirical profile.


### 3. Empirical priors

To characterize "normal variation" for each window:

1. Score **all sequences** from the alignment against the window profile.  
2. Collect the distribution of nLL values.  
3. Extract percentiles (e.g., 95th and 99th) to serve as thresholds.  

Thus, for each window we store:  
- The **distribution (profile)**.  
- Empirical thresholds: `nLL_p95` and `nLL_p99`.  

A new sequence can later be compared:  
- If `nLL < nLL_p95` → typical.  
- If `nLL > nLL_p99` → unusually variable, possibly unreliable region.

---

---

## � Supporting utils

Several utility scripts provide reusable functions for both processes:  

- **utils.py** → basic alignment and scoring functions:  
  - `load_alignment`, `extract_ref_positions`, `sliding_windows`, `score_window`.  

- **utils_integrate_consensus.py** → additional helpers for consensus integration:  
  - missingness and fragmentation counts,  
  - insertion handling,  
  - QC calculations,  
  - consensus merging,  
  - window evaluation wrapper.  

These modular functions keep the pipeline clean and reusable.


