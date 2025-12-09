#!/usr/bin/env python3
# cli.py
"""
priorCons CLI
"""
import sys

from . import __version__
from . import build_priors as bp
from . import integrate_consensus as ic
from . import qc_priorcons as qc 
from . import complete_with_gf as gf


def main(argv=None):
    argv = argv if argv is not None else sys.argv[1:]

    # Si el usuario no da subcomando o pide --help general
    if len(argv) == 0 or argv[0] in ("-h", "--help"):
        print(f"""
priorCons {__version__}

Usage:
  priorcons <subcommand> [options]

Available subcommands:
  build-priors         Build priors parquet file
  integrate-consensus  Run consensus integration workflow
  qc                   Run Quality Control and post-processing 
  complete-gf          Run complete with GF workflow

Use 'priorcons <subcommand> -h' for details on each one.
""")
        sys.exit(0)

    if argv[0] in ("--version", "-v", "-V"):
        print(f"priorcons {__version__}")
        sys.exit(0)

    # Delegar completamente al módulo correspondiente
    subcmd = argv[0]
    subargs = argv[1:]

    if subcmd == "build-priors":
        sys.exit(bp.main(subargs))
    elif subcmd == "integrate-consensus":
        sys.exit(ic.main(subargs))
    elif subcmd == "qc":
        sys.exit(qc.main(subargs)) 
    elif subcmd == "complete-gf":
        sys.exit(gf.main(subargs)) 
    else:
        print(f"Unknown command: {subcmd}")
        print("Use 'priorcons --help' for available commands.")
        sys.exit(1)


if __name__ == "__main__":
    main()
