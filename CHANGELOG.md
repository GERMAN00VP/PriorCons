# Changelog

All notable changes to this project will be documented in this file.

---

## [v0.1.0] - 2025-10-23
### First version upload

## [v0.1.1] - 2025-12-03 

### New Features and Improvements

    Integrated Quality Control (QC) Analysis: Added a new comprehensive QC module (priorcons qc) for post-processing and evaluation of results. This includes:

        Aggregation and plotting of Performance Metrics (qc.json).

        Identification and plotting of Sequence Recovery Hotspots (window_trace.csv).

    Built-in Data: Implemented the infrastructure to include static data files (like GFF annotations) within the package using priorcons/data/.

    CLI Enhancement: The main CLI (priorcons/cli.py) was updated to include the new qc sub-command, and the package entry points were reorganized in priorcons/__init__.py.

###  Bug Fixes and Technical Changes

    Fixed Dependencies/Configuration: Resolved minor configuration issues and updated dependencies in pyproject.toml.

    Removed Legacy Files: The redundant setup.cfg file was removed in favor of the standardized pyproject.toml.

    Cleanup of Sample Data: All large, specific sample files, prior parquet files, test data, and reference FASTA sequences were removed from the repository (rsv_files/ directory cleanup).

    Internal File Cleanup: Removed Python cache files (__pycache__/) and unused utility files (utils.cpython-312.pyc, etc.).

###  Refactoring

    Package Structure: Introduced .gitignore and consolidated package configuration for easier installation and distribution (PyPI/Bioconda/Containers).