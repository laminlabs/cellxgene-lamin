"""Manage cellxgene metadata using LaminDB.

Import the package::

   import cellxgene_lamin

This is the complete API reference:

.. autosummary::
   :toctree: .

   Curate
   CellxGeneFields
   datasets
"""

__version__ = "0.3.1"  # denote a pre-release for 0.1.0 with 0.1rc1

from .curate import Curate
from .dev import datasets
from .fields import CellxGeneFields
