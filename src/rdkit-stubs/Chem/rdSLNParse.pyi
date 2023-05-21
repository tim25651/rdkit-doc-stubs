"""
rdkit.Chem.rdSLNParse module¶
Module containing classes and functions for working with Sybyl line notation (SLN).
"""
from typing import Any

from rdkit.Chem.rdchem import Mol

def MolFromQuerySLN(
    self, SLN: str, mergeHs: bool = True, debugParser: bool = False
) -> Mol:
    """
    Construct a query molecule from an SLN string.

    ARGUMENTS:

    SLN: the SLN string
    mergeHs: (optional) toggles the merging of explicit Hs in the query into the attached
    heavy atoms. Defaults to False.

    RETURNS:

    a Mol object suitable for using in substructure queries, None on failure.

    C++ signature :RDKit::ROMol* MolFromQuerySLN(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=True [,bool=False]])
    """
    ...

def MolFromSLN(self, SLN: str, sanitize: bool = True, debugParser: bool = False) -> Mol:
    """
    Construct a molecule from an SLN string.

    ARGUMENTS:

    SLN: the SLN string
    sanitize: (optional) toggles sanitization of the molecule.
    Defaults to True.

    RETURNS:

    a Mol object, None on failure.

    NOTE: the SLN should not contain query information or properties. To build aquery from SLN, use MolFromQuerySLN.

    C++ signature :RDKit::ROMol* MolFromSLN(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=True [,bool=False]])
    """
    ...
