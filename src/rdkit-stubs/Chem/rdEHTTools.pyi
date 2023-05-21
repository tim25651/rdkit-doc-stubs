"""
rdkit.Chem.rdEHTTools module¶
Module containing interface to the YAeHMOP extended Hueckel library.
Please note that this interface should still be considered experimental and may
change from one release to the next.
"""
from typing import Any

import Boost.Python
from rdkit.Chem.rdchem import Mol

class EHTResults(Boost.Python.instance):
    """
    Raises an exception
    This class cannot be instantiated from Python"""

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def GetAtomicCharges(self, arg1: EHTResults) -> object:
        """
        returns the calculated atomic charges

        C++ signature :_object* GetAtomicCharges(RDKit::EHTTools::EHTResults {lvalue})
        """
        ...
    def GetHamiltonian(self, arg1: EHTResults) -> object:
        """
        returns the symmetric Hamiltonian matrix

        C++ signature :_object* GetHamiltonian(RDKit::EHTTools::EHTResults {lvalue})"""
        ...
    def GetOrbitalEnergies(self, arg1: EHTResults) -> object:
        """
        returns the energies of the molecular orbitals as a vector

        C++ signature :_object* GetOrbitalEnergies(RDKit::EHTTools::EHTResults {lvalue})
        """
        ...
    def GetOverlapMatrix(self, arg1: EHTResults) -> object:
        """
        returns the symmetric overlap matrix

        C++ signature :_object* GetOverlapMatrix(RDKit::EHTTools::EHTResults {lvalue})
        """
        ...
    def GetReducedChargeMatrix(self, arg1: EHTResults) -> object:
        """
        returns the reduced charge matrix

        C++ signature :_object* GetReducedChargeMatrix(RDKit::EHTTools::EHTResults {lvalue})
        """
        ...
    def GetReducedOverlapPopulationMatrix(self, arg1: EHTResults) -> object:
        """
        returns the reduced overlap population matrix

        C++ signature :_object* GetReducedOverlapPopulationMatrix(RDKit::EHTTools::EHTResults {lvalue})
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @property
    def fermiEnergy(self) -> Any: ...
    @property
    def numElectrons(self) -> Any: ...
    @property
    def numOrbitals(self) -> Any: ...
    @property
    def totalEnergy(self) -> Any: ...

def RunMol(
    self, mol: Mol, confId: int = -1, keepOverlapAndHamiltonianMatrices: bool = False
) -> tuple:
    """
    Runs an extended Hueckel calculation for a molecule.
    The molecule should have at least one conformation

    ARGUMENTS:
    mol: molecule to use
    confId: (optional) conformation to use
    keepOverlapAndHamiltonianMatrices: (optional) triggers storing the overlap
    and hamiltonian matrices in the EHTResults object

    RETURNS: a 2-tuple:
    a boolean indicating whether or not the calculation succeeded
    an EHTResults object with the results

    C++ signature :boost::python::tuple RunMol(RDKit::ROMol [,int=-1 [,bool=False]])"""
    ...
