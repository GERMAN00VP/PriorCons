# utils_qc.py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats
import seaborn as sns
import json
from pathlib import Path
from glob import glob
import os
import re
import importlib.resources

def get_builtin_gff_path() -> Path:
    """
    Returns the path to the packaged GFF file using importlib.resources.
    """

    resource_name = "rsv.gff"

    try:
        # Compatible with Bioconda, Containers, and PyPI
        return importlib.resources.files('priorcons.data') / resource_name
    except Exception as e:
        raise RuntimeError(f"Could not locate built-in GFF file ({resource_name}). Ensure package data is installed correctly. Error: {e}")

# --- Core Logic Functions (Hotspot Analysis) ---

def select_top_windows_nosolap(df, window_size=1000):
    """
    Selects non-overlapping windows of size window_size with the highest Recovery_score.
    """

    if df.empty:
        return pd.DataFrame(columns=["start", "end", "Recovery_score"])

    df = df.sort_values("start").reset_index(drop=True)
    selected_windows = []

    # Genomic coverage range covered by hotspot windows
    coverage_start = df["start"].min()
    coverage_end   = df["end"].max()
    genome_length  = coverage_end - coverage_start
    
    # Array of flags: True = free, False = occupied
    free_positions = np.ones(genome_length, dtype=bool)

    # Convert start/end to indices inside free_positions vector
    df_indices = df.copy()
    df_indices["start_idx"] = df_indices["start"] - coverage_start
    df_indices["end_idx"]   = df_indices["end"] - coverage_start

    while True:
        max_sum = -np.inf
        max_start_idx = None
        
        for start_idx in range(0, genome_length - window_size + 1):
            if not free_positions[start_idx:start_idx+window_size].all():
                continue

            overlapping = df_indices[
                (df_indices["end_idx"] > start_idx) & (df_indices["start_idx"] < start_idx + window_size)
            ]
            if overlapping.empty:
                continue

            score_sum = 0
            for _, row in overlapping.iterrows():
                overlap_start = max(start_idx, row["start_idx"])
                overlap_end   = min(start_idx + window_size, row["end_idx"])
                proportion = (overlap_end - overlap_start) / (row["end_idx"] - row["start_idx"])
                score_sum += row["Recovery_score"] * proportion

            if score_sum > max_sum:
                max_sum = score_sum
                max_start_idx = start_idx

        if max_start_idx is None:
            break

        free_positions[max_start_idx:max_start_idx+window_size] = False

        selected_windows.append({
            "start": max_start_idx + coverage_start,
            "end":   max_start_idx + coverage_start + window_size,
            "Recovery_score": max_sum
        })

    result_df = pd.DataFrame(selected_windows)
    if not result_df.empty:
        result_df = result_df.sort_values("Recovery_score", ascending=False).reset_index(drop=True)
    return result_df

def extract_gene_info(gff_path: Path) -> pd.DataFrame:
    """
    Reads the GFF file and extracts gene/CDS information.
    """
    try:
        gff = pd.read_csv(gff_path, sep="\t", comment="#",
                          names=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"],
                          on_bad_lines='skip')
    except Exception as e:
        print(f"Error reading GFF file: {e}")
        return pd.DataFrame()

    genes = gff[gff["type"].isin(["CDS", "gene"])].copy()

    def get_name(attr_str):
        attributes = {att.split("=")[0]: att.split("=")[1] for att in attr_str.split(";") if "=" in att}
        return attributes.get("gene", attributes.get("Name", attributes.get("ID", "")))

    genes["gene_name"] = genes["attributes"].apply(get_name)
    
    genes = genes.loc[genes["gene_name"] != ""].drop_duplicates(subset="gene_name")
    return genes.reset_index(drop=True)

def load_and_process_data(priorcons_path: Path, gff_path: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Loads window trace files, computes hotspots, and loads the GFF.
    Returns (hotspots_df, top_windows_df, genes_df).
    """
    windows_files = glob(str(priorcons_path / "**" / "windows_trace.csv"), recursive=True)
    if not windows_files:
        raise FileNotFoundError(f"No 'windows_trace.csv' files found in {priorcons_path}")
        
    windows_df = pd.concat([pd.read_csv(wind) for wind in windows_files])

    windows_df["Recovery"] = windows_df.MISSING_mapp - windows_df.MISSING_ABACAS

    valid_windows = windows_df.dropna()
    count_df = valid_windows.groupby(["start","end"])["Recovery"].count().reset_index(name="count")
    sum_df   = valid_windows.groupby(["start","end"])["Recovery"].sum().reset_index(name="sum")

    hotspots = count_df.merge(sum_df, on=["start","end"])

    hotspots["Recovery_score"] = hotspots["count"] * hotspots["sum"]
    hotspots = hotspots.sort_values("start")

    top_windows = select_top_windows_nosolap(hotspots, window_size=1000)

    genes = extract_gene_info(gff_path)

    return hotspots, top_windows, genes

def plot_hotspots(hotspots: pd.DataFrame, top_windows: pd.DataFrame, genes: pd.DataFrame, output_dir: Path):
    """
    Generates the recovery hotspot plot with gene annotations.
    """
    if hotspots.empty:
        print("No hotspot data to plot.")
        return

    plt.figure(figsize=(24, 6))

    plt.bar(
        hotspots["start"],
        hotspots["Recovery_score"],
        width=hotspots["end"] - hotspots["start"],
        color="steelblue",
        align="edge",
        label="Recovery score"
    )

    max_score = hotspots["Recovery_score"].max()
    colors = plt.get_cmap("tab10")
    if not genes.empty:
        for j, gene in enumerate(genes.itertuples()):
            plt.axvspan(gene.start, gene.end, color=colors(j % 10), alpha=0.3)
            if max_score > 0:
                plt.text((gene.start + gene.end) / 2, -0.05 * max_score, gene.gene_name,
                         ha="center", va="top", rotation=90, color=colors(j % 10))

    y_top = max_score * 1.05
    for i, row in top_windows.head(5).iterrows():
        plt.hlines(y=y_top, xmin=row["start"], xmax=row["end"], colors='red', linewidth=3,
                   label=f"Top {i+1} Window" if i == 0 else None)
        plt.vlines(x=row["start"], ymin=y_top - 0.02 * y_top, ymax=y_top + 0.02 * y_top, colors='red', linewidth=3)
        plt.vlines(x=row["end"], ymin=y_top - 0.02 * y_top, ymax=y_top + 0.02 * y_top, colors='red', linewidth=3)
        plt.text((row["start"] + row["end"])/2, y_top + 0.02 * y_top, str(i+1),
                 ha="center", va="bottom", color="red", fontsize=10, fontweight='bold')

    plt.xlabel("Genomic position (bp)")
    plt.ylabel("Recovery Score")
    plt.title("Sequence Recovery Hotspots with Gene Annotations")
    
    if max_score > 0:
        plt.ylim(bottom=-0.1 * max_score, top=y_top * 1.05)

    handles, labels = plt.gca().get_legend_handles_labels()
    
    if not genes.empty:
        unique_genes = genes["gene_name"].unique()
        gene_patches = [plt.Rectangle((0, 0), 1, 1, fc=colors(j % 10), alpha=0.3, ec="none")
                        for j in range(len(unique_genes))]
        handles += gene_patches
        labels += list(unique_genes)

    plt.legend(handles, labels, bbox_to_anchor=(1.01, 1), loc="upper left", ncol=1)
    plt.tight_layout(rect=[0, 0.05, 1, 1])

    output_path = output_dir / "recovery_hotspots_plot.png"
    plt.savefig(output_path)
    plt.close()
    print(f"Plot saved to {output_path}")

# --- New Performance Analysis Functions ---

def load_qc_data(priorcons_path: Path) -> pd.DataFrame:
    """
    Loads all qc.json files, consolidates them, and computes coverage metrics.
    """
    qc_files = glob(str(priorcons_path / "**" / "qc.json"), recursive=True)
    if not qc_files:
        raise FileNotFoundError(f"No 'qc.json' files found in {priorcons_path}")
        
    qc_out = pd.DataFrame([json.load(open(q, "r")) for q in qc_files])
    qc_out["File"] = qc_files
    qc_out = qc_out.set_index("File")

    # Metrics
    qc_out['COVERAGE_INCREASE_ABS'] = qc_out['FINAL_COVERAGE'] - qc_out['MAPPING_CONSENSUS_COVERAGE']
    qc_out['COVERAGE_INCREASE_PERC'] = (qc_out['COVERAGE_INCREASE_ABS'] /
                                        qc_out['MAPPING_CONSENSUS_COVERAGE'].replace(0, np.nan)) * 100
    qc_out['COVERAGE_INCREASE_PERC'] = qc_out['COVERAGE_INCREASE_PERC'].fillna(0)
    
    return qc_out

def plot_performance(qc_out: pd.DataFrame, output_dir: Path):
    """
    Generates the 6-panel performance analysis plot.
    """
    if qc_out.empty:
        print("No QC data to plot performance.")
        return
    
    # Statistical tests
    t_stat, p_value = stats.ttest_rel(qc_out['MAPPING_CONSENSUS_COVERAGE'],
                                      qc_out['FINAL_COVERAGE'], nan_policy='omit')
    w_stat, w_pvalue = stats.wilcoxon(qc_out['MAPPING_CONSENSUS_COVERAGE'],
                                     qc_out['FINAL_COVERAGE'], nan_policy='omit')

    fig = plt.figure(figsize=(20, 14))
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.35)

    # 1. SCATTER PLOT
    ax1 = fig.add_subplot(gs[0, 0])
    scatter = ax1.scatter(qc_out['MAPPING_CONSENSUS_COVERAGE'], 
                         qc_out['COVERAGE_INCREASE_PERC'], 
                         c=qc_out['MAPPING_CONSENSUS_COVERAGE'],
                         cmap='viridis', s=100, alpha=0.7, edgecolors='black')

    z = np.polyfit(qc_out['MAPPING_CONSENSUS_COVERAGE'],
                   qc_out['COVERAGE_INCREASE_PERC'], 1)
    p = np.poly1d(z)
    sorted_x = sorted(qc_out['MAPPING_CONSENSUS_COVERAGE'])
    ax1.plot(sorted_x, p(sorted_x), "r--", alpha=0.8, linewidth=2,
            label=f'Regression: y={z[0]:.3f}x+{z[1]:.2f}')

    ax1.set_xlabel('Initial Coverage (X)', fontweight='bold')
    ax1.set_ylabel('Relative Increase (%)', fontweight='bold')
    ax1.set_title('Efficiency vs Initial Coverage', fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3, linestyle='--')
    plt.colorbar(scatter, ax=ax1, label='Initial Coverage')

    # 2. PAIRED CONNECTION DIAGRAM
    ax2 = fig.add_subplot(gs[0, 1])
    for idx in qc_out.index:
        x_vals = [0, 1]
        y_vals = [qc_out.loc[idx, 'MAPPING_CONSENSUS_COVERAGE'], 
                  qc_out.loc[idx, 'FINAL_COVERAGE']]
        ax2.plot(x_vals, y_vals, 'grey', alpha=0.2, linewidth=0.5)

    box_data = [qc_out['MAPPING_CONSENSUS_COVERAGE'].dropna(),
                qc_out['FINAL_COVERAGE'].dropna()]
    ax2.boxplot(box_data, positions=[0, 1], widths=0.2, patch_artist=True,
                boxprops=dict(facecolor='lightblue', alpha=0.7, linewidth=1.5),
                medianprops=dict(color='red', linewidth=2.5))

    ax2.scatter([0, 1], 
               [qc_out['MAPPING_CONSENSUS_COVERAGE'].mean(),
                qc_out['FINAL_COVERAGE'].mean()],
               s=100, color='darkblue', marker='D', edgecolor='white', linewidth=1.5,
               label='Mean', zorder=5)

    ax2.set_xticks([0, 1])
    ax2.set_xticklabels(['Initial', 'Final'], fontweight='bold')
    ax2.set_ylabel('Coverage', fontweight='bold')
    ax2.set_title(f'Coverage Comparison (Paired)\nPaired t-test p-value: {p_value:.4f}\n'
                  f'Wilcoxon p-value: {w_pvalue:.4f}', fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3, linestyle='--')

    # 3. HISTOGRAM
    ax3 = fig.add_subplot(gs[0, 2])
    n, bins, patches = ax3.hist(qc_out['COVERAGE_INCREASE_ABS'], bins=15,
                               edgecolor='black', alpha=0.7, color='skyblue', linewidth=1.5)
    
    for i in range(len(patches)):
        patches[i].set_facecolor(plt.cm.viridis(i / len(patches)))
        
    mean_val = qc_out['COVERAGE_INCREASE_ABS'].mean()
    median_val = qc_out['COVERAGE_INCREASE_ABS'].median()
    ax3.axvline(mean_val, color='red', linestyle='--', linewidth=2.5, 
               label=f'Mean: {mean_val:.2f}')
    ax3.axvline(median_val, color='green', linestyle=':', linewidth=2.5,
               label=f'Median: {median_val:.2f}')
    ax3.set_xlabel('Absolute Coverage Increase', fontweight='bold')
    ax3.set_ylabel('Frequency', fontweight='bold')
    ax3.set_title('Distribution of Absolute Increase', fontweight='bold')
    ax3.legend()
    ax3.grid(True, alpha=0.3, linestyle='--')

    # 4. CORRELATION HEATMAP
    ax4 = fig.add_subplot(gs[1, 0])
    corr_cols = ['MAPPING_CONSENSUS_COVERAGE', 'FINAL_COVERAGE', 
                 'COVERAGE_INCREASE_PERC', 'OBS-EXP_SUBSTITUTIONS',
                 'MAPPING_CONSENSUS_SUBSTITUTIONS', 'FINAL_SUBSTITUTIONS']
    corr_matrix = qc_out[corr_cols].corr()

    short_names = ['Initial Cov', 'Final Cov', '% Increase', 'Obs-Exp Subs',
                   'Initial Subs', 'Final Subs']

    sns.heatmap(corr_matrix, ax=ax4, cmap='RdBu_r', vmin=-1, vmax=1, square=True,
                annot=True, fmt='.2f', annot_kws={'size': 10, 'weight': 'bold'},
                cbar_kws={'shrink': 0.8, 'label': 'Correlation'},
                linewidths=0.5, linecolor='white',
                xticklabels=short_names, yticklabels=short_names)

    ax4.set_xticklabels(short_names, rotation=45, ha='right', fontweight='bold')
    ax4.set_yticklabels(short_names, rotation=0, fontweight='bold')
    ax4.set_title('Correlation Matrix', fontweight='bold', pad=15)

    # 5. BAR PLOT (Top 10)
    ax5 = fig.add_subplot(gs[1, 1])
    top_10 = qc_out.nlargest(10, 'COVERAGE_INCREASE_PERC')
    x_pos = np.arange(len(top_10))
    colors = plt.cm.YlGnBu(np.linspace(0.4, 0.9, len(top_10)))
    bars = ax5.barh(x_pos, top_10['COVERAGE_INCREASE_PERC'], 
                   color=colors, edgecolor='black', linewidth=1)

    sample_names = []
    for idx in top_10.index:
        parts = Path(idx).parent.name
        sample_name = parts if len(parts) <= 15 else parts[:12] + '...'
        sample_names.append(sample_name)

    ax5.set_yticks([])

    for i, (bar, name) in enumerate(zip(bars, sample_names)):
        text_x = 0.5
        text_y = bar.get_y() + bar.get_height()/2
        
        ax5.text(text_x, text_y, name, va='center', ha='left', fontsize=9.5, 
                fontweight='bold', color='black',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7, edgecolor='none'))

    for i, bar in enumerate(bars):
        width = bar.get_width()
        ax5.text(width + 0.5, bar.get_y() + bar.get_height()/2, 
                f'+{width:.1f}%', va='center', ha='left', fontsize=10, 
                fontweight='bold', color='darkgreen')

    ax5.set_xlabel('Relative Increase (%)', fontweight='bold')
    ax5.set_title('Top 10 Samples with Greatest Improvement', fontweight='bold', pad=15)
    ax5.grid(True, alpha=0.3, axis='x', linestyle='--')
    ax5.set_xlim([0, top_10['COVERAGE_INCREASE_PERC'].max() * 1.2])
    ax5.set_facecolor('#f8f9fa')

    # 6. STATISTICAL SUMMARY
    ax6 = fig.add_subplot(gs[1, 2])
    ax6.axis('off')

    samples_improved = len(qc_out[qc_out["COVERAGE_INCREASE_ABS"] > 0])
    success_rate = samples_improved/len(qc_out)*100
    corr_initial_increase = qc_out['MAPPING_CONSENSUS_COVERAGE'].corr(qc_out['COVERAGE_INCREASE_PERC'])
    substitution_change = qc_out['FINAL_SUBSTITUTIONS'].mean() - qc_out['MAPPING_CONSENSUS_SUBSTITUTIONS'].mean()

    summary_text = f"""
STATISTICAL SUMMARY
=============================================

COVERAGE METRICS:
• Initial (mean): {qc_out['MAPPING_CONSENSUS_COVERAGE'].mean():.2f} ± {qc_out['MAPPING_CONSENSUS_COVERAGE'].std():.2f}
• Final (mean): {qc_out['FINAL_COVERAGE'].mean():.2f} ± {qc_out['FINAL_COVERAGE'].std():.2f}
• Absolute increase: {qc_out['COVERAGE_INCREASE_ABS'].mean():.2f} ± {qc_out['COVERAGE_INCREASE_ABS'].std():.2f}
• Relative increase: {qc_out['COVERAGE_INCREASE_PERC'].mean():.2f}% ± {qc_out['COVERAGE_INCREASE_PERC'].std():.2f}%

STATISTICAL TESTS:
• Paired t-test: p-value = {p_value:.6f}
• Wilcoxon signed-rank: p-value = {w_pvalue:.6f}
• Initial-Increase correlation: r = {corr_initial_increase:.3f}

TOOL PERFORMANCE:
• Improved samples: {samples_improved}/{len(qc_out)} ({success_rate:.1f}%)
• Substitution change: {substitution_change:.2f} (Negative = improvement)
"""

    ax6.text(0.5, 0.98, summary_text, transform=ax6.transAxes,
            fontsize=10.5, verticalalignment='top', horizontalalignment='center',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9,
                     edgecolor='gold', linewidth=2),
            fontfamily='monospace')

    plt.suptitle('COMPREHENSIVE ANALYSIS: PriorCons Performance', 
                fontsize=18, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    output_path = output_dir / "performance_analysis_plot.png"
    plt.savefig(output_path)
    plt.close()
    print(f"Performance plot saved to {output_path}")
