# PriorCons

**Prior‑guided consensus integration for viral genomes**

---

## 🧭 Introduction

PriorCons improves viral consensus sequences by safely recovering missing information while preserving reliability. 

The software integrates:
* A **high‑confidence consensus sequence** (FASTA) generated using a stringent pipeline. This sequence is trusted but may contain masked regions (Ns).
* The **reference genome** used during assembly.
* A **candidate consensus sequence** that is less conservative but potentially more informative (for example, produced with relaxed filtering or alternative assembly).

The objective is to fill gaps in the high‑confidence consensus using information from the candidate sequence — but only when supported by evolutionary evidence — so that coverage increases without introducing sequencing artefacts.

To achieve this, PriorCons uses **evolutionary priors** derived from large collections of genomes for the same virus or subtype aligned to the reference. These priors model expected variation and provide statistical thresholds that guide integration decisions.

---

## 📦 Installation

PriorCons can be installed via **Conda** (recommended for bioinformatics) or **PyPI**:

### Using Conda
```bash
conda install -c bioconda priorcons
```
[View on Bioconda](https://anaconda.org/bioconda/priorcons)

### Using Pip
```bash
pip install priorcons
```
[View on PyPI](https://pypi.org/project/priorcons/)

---

## ⚡ Quickstart + CLI Examples

Follow these steps to generate an integrated consensus using PriorCons.

### 1. Prepare the Priors Database
You need a collection of viral sequences (e.g., from GISAID or NCBI) relevant to your sample.
* **Alignment is critical:** Use MAFFT in reference-anchored mode (e.g. `--add --keeplength`) to keep coordinates consistent when building priors.
* **Include the Reference:** Ensure your reference sequence is included in this FASTA file.

### 2. Build the Priors
Run the build-priors command to create the empirical distribution of variation.

```bash
priorcons build-priors --input database_aligned.fasta --output virus_priors.json
```

### 3. Run integrate-consensus
Once you have the priors, align your three sequences (Trusted, Candidate, and Reference) and run the integration.

**Alignment Recommendation:** Since you are only aligning 3 sequences, use a high-sensitivity strategy. We recommend **MAFFT** with the following parameters:

```bash
mafft --localpair --maxiterate 1000 input.fasta > aligned_input.fasta
```
**Running the integration:**

```bash
priorcons integrate-consensus \
    --aligned-fasta aligned_input.fasta \
    --priors virus_priors.json \
    --output integrated_consensus.fasta
```
---

## 🔬 Workflow Overview

*PriorCons uses a window-based approach to statistically validate and fill gaps in viral assemblies.*

1.  **Slide** overlapping windows across the genome.
2.  **Detect** windows with missing regions (Ns) in the trusted consensus.
3.  **Evaluate** the corresponding candidate window using the priors.
4.  **Accept** candidate window only if the score is evolutionarily plausible (below the statistical threshold).
5.  **Produce** an integrated consensus with increased completeness and maintained accuracy.

---

## 🧮 Methodology

### 1. Probability distributions per position
For each window of size $W$ bases, and each position $j$:

$$P_j(b)=\frac{c_j(b)+\alpha}{\sum_{x\in\{A,C,G,T\}}(c_j(x)+\alpha)}$$

Where:
* $c_j(b)$ is the count of base $b$.
* $\alpha$ is a pseudocount.
* Bases N are ignored.

### 2. Log‑likelihood of a sequence
Given a sequence $Q$:

$$\log L(Q \mid \text{window}) = \sum_j \log P_j(q_j)$$

Normalized negative log‑likelihood:

$$\text{nLL}(Q) = -\frac{1}{N_{\text{valid}}} \sum_j \log P_j(q_j)$$

Lower values indicate sequences consistent with expected variation.

### 3. Empirical thresholds
All sequences are scored to obtain an nLL distribution. The 95th percentile is used as a cutoff: windows exceeding this threshold are considered atypical and rejected during integration.

---

## 📊 Outputs

* **Integrated consensus FASTA:** The final integrated sequence.
* **Window‑level QC trace:** A file containing scores for each window.
* **Summary QC metrics:** Summary metrics regarding coverage and changes performed.

---

## 📚 Citing

This software was developed by Germán Vallejo Palma at the Instituto de Salud Carlos III
(ISCIII) — National Centre of Microbiology, Respiratory Viruses and Influenza Unit.

If you use this software in a publication, report, or product, please cite the
appropriate authors and include the above attribution.
