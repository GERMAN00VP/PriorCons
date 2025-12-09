"""
priorcons package
Author: Germán Vallejo Palma
Developed at: Instituto de Salud Carlos III - National Centre of Microbiology
"""

__version__ = "0.1.2"
__author__ = "Germán Vallejo Palma"
__email__ = "german.vallejo@isciii.es"

# Expose key functions/modules for convenience
from .build_priors import main as build_priors_main
from .integrate_consensus import main as integrate_consensus_main
from .qc_priorcons import main as qc_priorcons_main
from .complete_with_gf import main as complete_with_gf_main

