"""Manage cellxgene metadata using LaminDB.

Import the package::

   import cellxgene_lamin

"""

__version__ = "0.0.1"  # denote a pre-release for 0.1.0 with 0.1rc1

from ._cxg_rest import get_collections_from_cxg, get_datasets_from_cxg
