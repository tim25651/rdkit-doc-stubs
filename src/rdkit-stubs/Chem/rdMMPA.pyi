"""
rdkit.Chem.rdMMPA module¶
Module containing a C++ implementation of code for doing MMPA
"""
from typing import Any, overload

from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdMolDescriptors import AtomPairsParameters

@overload
def FragmentMol(
    self,
    mol: Mol,
    maxCuts: int = 3,
    maxCutBonds: int = 20,
    pattern: str = "[#6+0;!$(*=, #[!#6])]!@! = !#[*]",
    resultsAsMols: bool = True,
) -> tuple:
    """
    Does the fragmentation necessary for an MMPA analysis

    C++ signature :boost::python::tuple FragmentMol(RDKit::ROMol,unsigned int,unsigned int,unsigned int [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’[#6+0;!$(=,#[!#6])]!@!=!#[]’ [,bool=True]])

    FragmentMol( (Mol)mol, (AtomPairsParameters)bondsToCut [, (int)minCuts=1 [, (int)maxCuts=3 [, (bool)resultsAsMols=True]]]) -> tuple :Does the fragmentation necessary for an MMPA analysis

    C++ signature :boost::python::tuple FragmentMol(RDKit::ROMol,boost::python::api::object [,unsigned int=1 [,unsigned int=3 [,bool=True]]])
    """
    ...

@overload
def FragmentMol(
    self,
    mol: Mol,
    minCuts: int,
    maxCuts: int,
    maxCutBonds: int,
    pattern: str = "[#6+0;!$(=,#[!#6])]!@!=!#[]",
    resultsAsMols: bool = True,
) -> tuple:
    """
    Does the fragmentation necessary for an MMPA analysis

    C++ signature :boost::python::tuple FragmentMol(RDKit::ROMol,boost::python::api::object [,unsigned int=1 [,unsigned int=3 [,bool=True]]])
    """
    ...

@overload
def FragmentMol(
    self,
    mol: Mol,
    bondsToCut: AtomPairsParameters,
    minCuts: int = 1,
    maxCuts: int = 3,
    resultsAsMols: bool = True,
) -> tuple: ...
