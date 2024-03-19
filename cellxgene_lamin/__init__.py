"""Manage cellxgene metadata using LaminDB.

Import the package::

   import cellxgene_lamin

This is the complete API reference:

.. autosummary::
   :toctree: .

   Validator
   CellxGeneFields
   datasets
"""

__version__ = "0.1.0"  # denote a pre-release for 0.1.0 with 0.1rc1

from . import datasets
from ._cxg_rest import get_collections_from_cxg, get_datasets_from_cxg
from ._fields import CellxGeneFields
from ._validator import Validator
