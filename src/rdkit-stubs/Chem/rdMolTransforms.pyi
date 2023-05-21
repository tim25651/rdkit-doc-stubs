"""
rdkit.Chem.rdMolTransforms module¶
Module containing functions to perform 3D operations like rotate and translate conformations
"""
from typing import Any

from rdkit.Chem.rdchem import Conformer, Mol
from rdkit.Chem.rdMolDescriptors import AtomPairsParameters
from rdkit.Geometry.rdGeometry import Point3D
from rdkit.rdBase import _vectd

def CanonicalizeConformer(
    self,
    conf: Conformer,
    center: Point3D = None,
    normalizeCovar: bool = False,
    ignoreHs: bool = True,
) -> None:
    """
    Canonicalize the orientation of a conformer so that its principal axes
    around the specified center point coincide with the x, y, z axes

    ARGUMENTS:
    conf : conformer of interest

    centeroptionally center point about which the principal axes are computed if not specified the centroid of the conformer will be used

    normalizeCovar : Optionally normalize the covariance matrix by the number of atoms

    C++ signature :void CanonicalizeConformer(RDKit::Conformer {lvalue} [,RDGeom::Point3D const*=None [,bool=False [,bool=True]]])
    """
    ...

def CanonicalizeMol(
    self, mol: Mol, normalizeCovar: bool = False, ignoreHs: bool = True
) -> None:
    """
    Loop over the conformers in a molecule and canonicalize their orientation

    C++ signature :void CanonicalizeMol(RDKit::ROMol {lvalue} [,bool=False [,bool=True]])
    """
    ...

def ComputeCanonicalTransform(
    self,
    conf: Conformer,
    center: Point3D = None,
    normalizeCovar: bool = False,
    ignoreHs: bool = True,
) -> object:
    """
    Compute the transformation required aligna conformer so that
    the principal axes align up with the x,y, z axes
    The conformer itself is left unchanged

    ARGUMENTS:
    conf : the conformer of interest
    center : optional center point to compute the principal axes around (defaults to the centroid)
    normalizeCovar : optionally normalize the covariance matrix by the number of atoms

    C++ signature :_object* ComputeCanonicalTransform(RDKit::Conformer [,RDGeom::Point3D const*=None [,bool=False [,bool=True]]])
    """
    ...

def ComputeCentroid(
    self, conf: Conformer, ignoreHs: bool = True, weights: _vectd = None
) -> Point3D:
    """
    Compute the centroid of the conformation - hydrogens are ignored and no attentionis paid to the difference in sizes of the heavy atoms; however,
    an optional vector of weights can be passed.

    C++ signature :RDGeom::Point3D ComputeCentroid(RDKit::Conformer [,bool=True [,std::vector<double, std::allocator<double> > const*=None]])
    """
    ...

def ComputePrincipalAxesAndMoments(
    self, conf: Conformer, ignoreHs: bool = True, weights: AtomPairsParameters = None
) -> object:
    """
    Compute principal axes and moments of inertia for a conformer
    These values are calculated from the inertia tensor:
    Iij = - sum_{s=1..N}(w_s * r_{si} * r_{sj}) i != j
    Iii = sum_{s=1..N} sum_{j!=i} (w_s * r_{sj} * r_{sj})
    where the coordinates are relative to the center of mass.

    ARGUMENTS:
    conf : the conformer of interest
    ignoreHs : if True, ignore hydrogen atoms
    weights : if present, used to weight the atomic coordinates

    Returns a (principal axes, principal moments) tuple

    C++ signature :_object* ComputePrincipalAxesAndMoments(RDKit::Conformer [,bool=True [,boost::python::api::object=None]])
    """
    ...

def ComputePrincipalAxesAndMomentsFromGyrationMatrix(
    self, conf: Conformer, ignoreHs: bool = True, weights: AtomPairsParameters = None
) -> object:
    """
    Compute principal axes and moments from the gyration matrix of a conformer
    These values are calculated from the gyration matrix/tensor:
    Iij = sum_{s=1..N}(w_s * r_{si} * r_{sj}) i != j
    Iii = sum_{s=1..N} sum_{t!=s}(w_s * r_{si} * r_{ti})
    where the coordinates are relative to the center of mass.

    ARGUMENTS:
    conf : the conformer of interest
    ignoreHs : if True, ignore hydrogen atoms
    weights : if present, used to weight the atomic coordinates

    Returns a (principal axes, principal moments) tuple

    C++ signature :_object* ComputePrincipalAxesAndMomentsFromGyrationMatrix(RDKit::Conformer [,bool=True [,boost::python::api::object=None]])
    """
    ...

def GetAngleDeg(
    self, conf: Conformer, iAtomId: int, jAtomId: int, kAtomId: int
) -> float:
    """
    Returns the angle in degrees between atoms i, j, k

    C++ signature :double GetAngleDeg(RDKit::Conformer,unsigned int,unsigned int,unsigned int)
    """
    ...

def GetAngleRad(
    self, conf: Conformer, iAtomId: int, jAtomId: int, kAtomId: int
) -> float:
    """
    Returns the angle in radians between atoms i, j, k

    C++ signature :double GetAngleRad(RDKit::Conformer,unsigned int,unsigned int,unsigned int)
    """
    ...

def GetBondLength(self, conf: Conformer, iAtomId: int, jAtomId: int) -> float:
    """
    Returns the bond length in angstrom between atoms i, j

    C++ signature :double GetBondLength(RDKit::Conformer,unsigned int,unsigned int)"""
    ...

def GetDihedralDeg(
    self, conf: Conformer, iAtomId: int, jAtomId: int, kAtomId: int, lAtomId: int
) -> float:
    """
    Returns the dihedral angle in degrees between atoms i, j, k, l

    C++ signature :double GetDihedralDeg(RDKit::Conformer,unsigned int,unsigned int,unsigned int,unsigned int)
    """
    ...

def GetDihedralRad(
    self, conf: Conformer, iAtomId: int, jAtomId: int, kAtomId: int, lAtomId: int
) -> float:
    """
    Returns the dihedral angle in radians between atoms i, j, k, l

    C++ signature :double GetDihedralRad(RDKit::Conformer,unsigned int,unsigned int,unsigned int,unsigned int)
    """
    ...

def SetAngleDeg(
    self, conf: Conformer, iAtomId: int, jAtomId: int, kAtomId: int, value: float
) -> None:
    """
    Sets the angle in degrees between atoms i, j, k; all atoms bonded to atom k are moved

    C++ signature :void SetAngleDeg(RDKit::Conformer {lvalue},unsigned int,unsigned int,unsigned int,double)
    """
    ...

def SetAngleRad(
    self, conf: Conformer, iAtomId: int, jAtomId: int, kAtomId: int, value: float
) -> None:
    """
    Sets the angle in radians between atoms i, j, k; all atoms bonded to atom k are moved

    C++ signature :void SetAngleRad(RDKit::Conformer {lvalue},unsigned int,unsigned int,unsigned int,double)
    """
    ...

def SetBondLength(
    self, conf: Conformer, iAtomId: int, jAtomId: int, value: float
) -> None:
    """
    Sets the bond length in angstrom between atoms i, j; all atoms bonded to atom j are moved

    C++ signature :void SetBondLength(RDKit::Conformer {lvalue},unsigned int,unsigned int,double)
    """
    ...

def SetDihedralDeg(
    self, arg1: Conformer, arg2: int, arg3: int, arg4: int, arg5: int, arg6: float
) -> None:
    """
    Sets the dihedral angle in degrees between atoms i, j, k, l; all atoms bonded to atom l are moved

    C++ signature :void SetDihedralDeg(RDKit::Conformer {lvalue},unsigned int,unsigned int,unsigned int,unsigned int,double)
    """
    ...

def SetDihedralRad(
    self,
    conf: Conformer,
    iAtomId: int,
    jAtomId: int,
    kAtomId: int,
    lAtomId: int,
    value: float,
) -> None:
    """
    Sets the dihedral angle in radians between atoms i, j, k, l; all atoms bonded to atom l are moved

    C++ signature :void SetDihedralRad(RDKit::Conformer {lvalue},unsigned int,unsigned int,unsigned int,unsigned int,double)
    """
    ...

def TransformConformer(self, arg1: Conformer, arg2: AtomPairsParameters) -> None:
    """
    Transform the coordinates of a conformer

    C++ signature :void TransformConformer(RDKit::Conformer {lvalue},boost::python::api::object)
    """
    ...
