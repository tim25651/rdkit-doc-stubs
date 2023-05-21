"""
Module containing functions to handle force fields
"""
from typing import Any

from rdkit.Chem.rdchem import Mol
from rdkit.ForceField.rdForceField import ForceField, MMFFMolProperties

def GetUFFAngleBendParams(self, mol: Mol, idx1: int, idx2: int, idx3: int) -> object:
    """
    Retrieves UFF angle bend parameters for atoms with indexes idx1, idx2, idx3 as a (ka, theta0) tuple, or None if no parameters could be found

    C++ signature :_object* GetUFFAngleBendParams(RDKit::ROMol,unsigned int,unsigned int,unsigned int)
    """
    ...

def GetUFFBondStretchParams(self, mol: Mol, idx1: int, idx2: int) -> object:
    """
    Retrieves UFF bond stretch parameters for atoms with indexes idx1, idx2 as a (kb, r0) tuple, or None if no parameters could be found

    C++ signature :_object* GetUFFBondStretchParams(RDKit::ROMol,unsigned int,unsigned int)
    """
    ...

def GetUFFInversionParams(
    self, mol: Mol, idx1: int, idx2: int, idx3: int, idx4: int
) -> object:
    """
    Retrieves UFF inversion parameters for atoms with indexes idx1, idx2, idx3, idx4 as a K float value, or None if no parameters could be found

    C++ signature :_object* GetUFFInversionParams(RDKit::ROMol,unsigned int,unsigned int,unsigned int,unsigned int)
    """
    ...

def GetUFFTorsionParams(
    self, mol: Mol, idx1: int, idx2: int, idx3: int, idx4: int
) -> object:
    """
    Retrieves UFF torsion parameters for atoms with indexes idx1, idx2, idx3, idx4 as a V float value, or None if no parameters could be found

    C++ signature :_object* GetUFFTorsionParams(RDKit::ROMol,unsigned int,unsigned int,unsigned int,unsigned int)
    """
    ...

def GetUFFVdWParams(self, mol: Mol, idx1: int, idx2: int) -> object:
    """
    Retrieves UFF van der Waals parameters for atoms with indexes idx1, idx2 as a (x_ij, D_ij) tuple, or None if no parameters could be found

    C++ signature :_object* GetUFFVdWParams(RDKit::ROMol,unsigned int,unsigned int)"""
    ...

def MMFFGetMoleculeForceField(
    self,
    mol: Mol,
    pyMMFFMolProperties: MMFFMolProperties,
    nonBondedThresh: float = 100.0,
    confId: int = -1,
    ignoreInterfragInteractions: bool = True,
) -> ForceField:
    """
    returns a MMFF force field for a molecule

    ARGUMENTS:

    mol : the molecule of interest

    pyMMFFMolPropertiesPyMMFFMolProperties object as returnedby MMFFGetMoleculeProperties()

    nonBondedThreshused to exclude long-range non-bondedinteractions (defaults to 100.0)

    confId : indicates which conformer to optimize

    ignoreInterfragInteractionsif true, nonbonded terms betweenfragments will not be added to the forcefield

    C++ signature :ForceFields::PyForceField* MMFFGetMoleculeForceField(RDKit::ROMol {lvalue},ForceFields::PyMMFFMolProperties* [,double=100.0 [,int=-1 [,bool=True]]])
    """
    ...

def MMFFGetMoleculeProperties(
    self, mol: Mol, mmffVariant: str = "MMFF94", mmffVerbosity: int = 0
) -> MMFFMolProperties:
    """
    returns a PyMMFFMolProperties object for amolecule, which is required by MMFFGetMoleculeForceField()
    and can be used to get/set MMFF properties
    ARGUMENTS:

    mol : the molecule of interest

    mmffVariant“MMFF94” or “MMFF94s”(defaults to “MMFF94”)

    mmffVerbosity : 0: none; 1: low; 2: high (defaults to 0).

    C++ signature :ForceFields::PyMMFFMolProperties* MMFFGetMoleculeProperties(RDKit::ROMol {lvalue} [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’MMFF94’ [,unsigned int=0]])
    """
    ...

def MMFFHasAllMoleculeParams(self, mol: Mol) -> bool:
    """
    checks if MMFF parameters are available for all of a molecule’s atoms

    ARGUMENTS:

    mol : the molecule of interest

    C++ signature :bool MMFFHasAllMoleculeParams(RDKit::ROMol)"""
    ...

def MMFFOptimizeMolecule(
    self,
    mol: Mol,
    mmffVariant: str = "MMFF94",
    maxIters: int = 200,
    nonBondedThresh: float = 100.0,
    confId: int = -1,
    ignoreInterfragInteractions: bool = True,
) -> int:
    """
    uses MMFF to optimize a molecule’s structure

    ARGUMENTS:

    mol : the molecule of interest
    mmffVariant : “MMFF94” or “MMFF94s”
    maxIters : the maximum number of iterations (defaults to 200)

    nonBondedThreshused to exclude long-range non-bondedinteractions (defaults to 100.0)

    confId : indicates which conformer to optimize

    ignoreInterfragInteractionsif true, nonbonded terms betweenfragments will not be added to the forcefield

    RETURNS: 0 if the optimization converged, -1 if the forcefield couldnot be set up, 1 if more iterations are required.

    C++ signature :int MMFFOptimizeMolecule(RDKit::ROMol {lvalue} [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’MMFF94’ [,int=200 [,double=100.0 [,int=-1 [,bool=True]]]]])
    """
    ...

def MMFFOptimizeMoleculeConfs(
    self: Mol,
    numThreads: int = 1,
    maxIters: int = 200,
    mmffVariant: str = "MMFF94",
    nonBondedThresh: float = 100.0,
    ignoreInterfragInteractions: bool = True,
) -> object:
    """
    uses MMFF to optimize all of a molecule’s conformations

    ARGUMENTS:

    mol : the molecule of interest

    numThreadsthe number of threads to use, only has an effect if the RDKitwas built with thread support (defaults to 1)
    If set to zero, the max supported by the system will be used.

    maxIters : the maximum number of iterations (defaults to 200)
    mmffVariant : “MMFF94” or “MMFF94s”

    nonBondedThreshused to exclude long-range non-bondedinteractions (defaults to 100.0)

    ignoreInterfragInteractionsif true, nonbonded terms betweenfragments will not be added to the forcefield.

    RETURNS: a list of (not_converged, energy) 2-tuples. If not_converged is 0 the optimization converged for that conformer.

    C++ signature :boost::python::api::object MMFFOptimizeMoleculeConfs(RDKit::ROMol {lvalue} [,int=1 [,int=200 [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’MMFF94’ [,double=100.0 [,bool=True]]]]])
    """
    ...

def MMFFSanitizeMolecule(self, mol: Mol) -> int:
    """
    sanitizes a molecule according to MMFF requirements.

    mol : the molecule of interest.

    C++ signature :unsigned int MMFFSanitizeMolecule(RDKit::ROMol {lvalue})"""
    ...

def OptimizeMolecule(self, ff: ForceField, maxIters: int = 200) -> int:
    """
    uses the supplied force field to optimize a molecule’s structure

    ARGUMENTS:

    ff : the force field
    maxIters : the maximum number of iterations (defaults to 200)

    RETURNS: 0 if the optimization converged, 1 if more iterations are required.

    C++ signature :int OptimizeMolecule(ForceFields::PyForceField {lvalue} [,int=200])
    """
    ...

def OptimizeMoleculeConfs(
    self, mol: Mol, ff: ForceField, numThreads: int = 1, maxIters: int = 200
) -> object:
    """
    uses the supplied force field to optimize all of a molecule’s conformations

    ARGUMENTS:

    mol : the molecule of interest
    ff : the force field

    numThreadsthe number of threads to use, only has an effect if the RDKitwas built with thread support (defaults to 1)
    If set to zero, the max supported by the system will be used.

    maxIters : the maximum number of iterations (defaults to 200)

    RETURNS: a list of (not_converged, energy) 2-tuples. If not_converged is 0 the optimization converged for that conformer.

    C++ signature :boost::python::api::object OptimizeMoleculeConfs(RDKit::ROMol {lvalue},ForceFields::PyForceField {lvalue} [,int=1 [,int=200]])
    """
    ...

def UFFGetMoleculeForceField(
    self,
    mol: Mol,
    vdwThresh: float = 10.0,
    confId: int = -1,
    ignoreInterfragInteractions: bool = True,
) -> ForceField:
    """
    returns a UFF force field for a molecule

    ARGUMENTS:

    mol : the molecule of interest

    vdwThreshused to exclude long-range van der Waals interactions(defaults to 10.0)

    confId : indicates which conformer to optimize

    ignoreInterfragInteractionsif true, nonbonded terms betweenfragments will not be added to the forcefield.

    C++ signature :ForceFields::PyForceField* UFFGetMoleculeForceField(RDKit::ROMol {lvalue} [,double=10.0 [,int=-1 [,bool=True]]])
    """
    ...

def UFFHasAllMoleculeParams(self, mol: Mol) -> bool:
    """
    checks if UFF parameters are available for all of a molecule’s atoms

    ARGUMENTS:

    mol : the molecule of interest.

    C++ signature :bool UFFHasAllMoleculeParams(RDKit::ROMol)"""
    ...

def UFFOptimizeMolecule(
    self: Mol,
    maxIters: int = 200,
    vdwThresh: float = 10.0,
    confId: int = -1,
    ignoreInterfragInteractions: bool = True,
) -> int:
    """
    uses UFF to optimize a molecule’s structure

    ARGUMENTS:

    mol : the molecule of interest
    maxIters : the maximum number of iterations (defaults to 200)

    vdwThreshused to exclude long-range van der Waals interactions(defaults to 10.0)

    confId : indicates which conformer to optimize

    ignoreInterfragInteractionsif true, nonbonded terms betweenfragments will not be added to the forcefield.

    RETURNS: 0 if the optimization converged, 1 if more iterations are required.

    C++ signature :int UFFOptimizeMolecule(RDKit::ROMol {lvalue} [,int=200 [,double=10.0 [,int=-1 [,bool=True]]]])
    """
    ...

def UFFOptimizeMoleculeConfs(
    self: Mol,
    numThreads: int = 1,
    maxIters: int = 200,
    vdwThresh: float = 10.0,
    ignoreInterfragInteractions: bool = True,
) -> object:
    """
    uses UFF to optimize all of a molecule’s conformations

    ARGUMENTS:

    mol : the molecule of interest

    numThreadsthe number of threads to use, only has an effect if the RDKitwas built with thread support (defaults to 1)
    If set to zero, the max supported by the system will be used.

    maxIters : the maximum number of iterations (defaults to 200)

    vdwThreshused to exclude long-range van der Waals interactions(defaults to 10.0)

    ignoreInterfragInteractionsif true, nonbonded terms betweenfragments will not be added to the forcefield.

    RETURNS: a list of (not_converged, energy) 2-tuples. If not_converged is 0 the optimization converged for that conformer.

    C++ signature :boost::python::api::object UFFOptimizeMoleculeConfs(RDKit::ROMol {lvalue} [,int=1 [,int=200 [,double=10.0 [,bool=True]]]])
    """
    ...
