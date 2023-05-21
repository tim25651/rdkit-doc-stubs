"""
rdkit.Chem.rdReducedGraphs module¶
Module containing functions to generate and work with reduced graphs
"""
from typing import Any

from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdMolDescriptors import AtomPairsParameters

def GenerateErGFingerprintForReducedGraph(
    self,
    mol: Mol,
    atomTypes: AtomPairsParameters = 0,
    fuzzIncrement: float = 0.3,
    minPath: int = 1,
    maxPath: int = 15,
) -> object:
    """
    Returns the ErG fingerprint vector for a reduced graph

    C++ signature :_object* GenerateErGFingerprintForReducedGraph(RDKit::ROMol [,boost::python::api::object=0 [,double=0.3 [,int=1 [,int=15]]]])
    """
    ...

def GenerateMolExtendedReducedGraph(
    self, mol: Mol, atomTypes: AtomPairsParameters = 0
) -> Mol:
    """
    Returns the reduced graph for a molecule

    C++ signature :RDKit::ROMol* GenerateMolExtendedReducedGraph(RDKit::ROMol [,boost::python::api::object=0])
    """
    ...

def GetErGFingerprint(
    self,
    mol: Mol,
    atomTypes: AtomPairsParameters = 0,
    fuzzIncrement: float = 0.3,
    minPath: int = 1,
    maxPath: int = 15,
) -> object:
    """
    Returns the ErG fingerprint vector for a molecule

    C++ signature :_object* GetErGFingerprint(RDKit::ROMol [,boost::python::api::object=0 [,double=0.3 [,int=1 [,int=15]]]])
    """
    ...
