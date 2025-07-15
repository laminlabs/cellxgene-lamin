"""Manage cellxgene metadata using LaminDB.

Import the package::

   import cellxgene_lamin

This is the complete API reference:

.. autosummary::
   :toctree: .

   datasets
"""

__version__ = "0.3.2"  # denote a pre-release for 0.1.0 with 0.1rc1

from lamindb.curators import CellxGeneAnnDataCatManager as Curator
from lamindb.curators import CellxGeneFields

from .dev import datasets
