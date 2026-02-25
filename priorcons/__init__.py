"""
priorcons package
Author: Germán Vallejo Palma
Developed at: Instituto de Salud Carlos III - National Centre of Microbiology
"""

__version__ = "0.1.4"
__author__ = "Germán Vallejo Palma"
__email__ = "german.vallejo@isciii.es"

# Expose key functions/modules for convenience
from .build_priors import main as build_priors_main
from .integrate_consensus import main as integrate_consensus_main