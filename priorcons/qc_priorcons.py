#!/usr/bin/env python3
# qc_priorcons.py
"""
CLI script to run post-processing and Quality Control (QC) analysis for PriorCons results.
Generates both Hotspot (window_trace.csv aggregation) and Performance (qc.json aggregation) plots.

This script is executed via the main CLI entry point:
priorcons qc --input_dir <path> --gff_file <path> --output_dir <path>
"""
import argparse
import sys
import logging
from pathlib import Path
from typing import List

# Import helper functions from utils_qc.py
try:
    from .utils_qc import (
        load_and_process_data, plot_hotspots, 
        load_qc_data, plot_performance,get_builtin_gff_path
    )
except ImportError:
    # Fallback para pruebas directas si no está como paquete
    from utils_qc import (
        load_and_process_data, plot_hotspots, 
        load_qc_data, plot_performance,get_builtin_gff_path
    )

# Configuración básica del logger
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger("qc_priorcons")

# --- Funciones Core de Ejecución ---

def run_hotspots_analysis_core(input_dir: Path, gff_file: Path, output_dir: Path):
    """Ejecuta la lógica de agregación y ploteo de Hotspots de recuperación."""
    logger.info("Starting Hotspots analysis...")
    
    # 1. Cargar y calcular hotspots
    logger.info("Loading and processing window trace files...")
    hotspots_df, top_windows_df, genes_df = load_and_process_data(input_dir, gff_file)
    
    if hotspots_df.empty:
        logger.warning("No valid data found to calculate hotspots.")
        return 0
    
    # Guardar datos intermedios
    hotspots_df.to_csv(output_dir / "hotspots_aggregated.csv", index=False)
    top_windows_df.to_csv(output_dir / "top_windows_nosolap.csv", index=False)
    
    # 2. Generar plot
    logger.info("Generating hotspots plot...")
    plot_hotspots(hotspots_df, top_windows_df, genes_df, output_dir)
    logger.info("Hotspots analysis finished.")
    return 0

def run_performance_analysis_core(input_dir: Path, output_dir: Path):
    """Ejecuta la lógica de agregación y ploteo del rendimiento de la herramienta."""
    logger.info("Starting Performance analysis...")

    # 1. Cargar y procesar QC JSON files
    logger.info("Loading and processing QC JSON files...")
    qc_out = load_qc_data(input_dir)
    
    if qc_out.empty:
        logger.warning("No QC data found.")
        return 0
        
    # Guardar datos intermedios
    qc_out.to_csv(output_dir / "qc_aggregated_metrics.csv")
    
    # 2. Generar plot
    logger.info("Generating performance analysis plot...")
    plot_performance(qc_out, output_dir)
    logger.info("Performance analysis finished.")
    return 0


def main(argv: List[str] = None) -> int:
    """Función principal para el subcomando 'qc', que ejecuta ambos análisis."""
    argv = argv if argv is not None else sys.argv[1:]
    
    p = argparse.ArgumentParser(
        description="Run PriorCons Quality Control (QC) and Post-processing analysis, generating both Hotspot and Performance plots.",
        usage="priorcons qc --input_dir <path> --gff_file <path> --output_dir <path>"
    )
    p.add_argument(
        "--input_dir", required=True, type=Path,
        help="Root directory containing PriorCons results. Requires both window_trace.csv and qc.json files."
    )
    p.add_argument(
            "--gff_file", type=Path, default=None, # Hazlo opcional, con un default None
            help="Path to the GFF/GFF3 file (e.g., rsv.gff). If not provided, a built-in GFF will be used."
        )
    p.add_argument(
        "--output_dir", required=True, type=Path,
        help="Output directory to write the generated plots and data."
    )
    
    args = p.parse_args(argv)
    
    input_dir: Path = args.input_dir
    gff_path_input: Path = args.gff_file
    output_dir: Path = args.output_dir
        
    try:
        # Validaciones
        if not input_dir.is_dir():
            logger.error("Input directory not found: %s", input_dir)
            sys.exit(2)

        if gff_path_input is None:
            # Usa la función para obtener la ruta del GFF empaquetado
            gff_path_to_use = get_builtin_gff_path()
            logger.info("Using built-in GFF file for annotation.")
        else:
            gff_path_to_use = gff_path_input

        if not gff_path_to_use.exists():
            logger.error("GFF file not found: %s", gff_path_to_use)
            sys.exit(2)

        output_dir.mkdir(parents=True, exist_ok=True)
        logger.info("Output directory: %s", output_dir)

        # 1. Ejecutar Análisis de Hotspots
        run_hotspots_analysis_core(input_dir, gff_path_to_use, output_dir)

        # 2. Ejecutar Análisis de Rendimiento
        run_performance_analysis_core(input_dir, output_dir)
        
        logger.info("PriorCons QC finished successfully (Hotspots and Performance plots generated).")
        return 0

    except Exception as e:
        logger.exception("PriorCons QC failed during execution: %s", e)
        return 1

if __name__ == "__main__":
    sys.exit(main())