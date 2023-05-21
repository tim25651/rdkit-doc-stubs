"""
rdkit.Chem.rdchem module¶
Module containing the core chemistry functionality of the RDKit
"""
from typing import Any, ClassVar, overload

import Boost.Python
from rdkit.Chem.rdchem import _ROAtomSeq, _ROBondSeq, _ROConformerSeq, _ROQAtomSeq
from rdkit.Chem.rdMolDescriptors import AtomPairsParameters
from rdkit.DataStructs.cDataStructs import ExplicitBitVect
from rdkit.Geometry.rdGeometry import Point3D
from rdkit.rdBase import (
    _vecti,
    _vectj,
    _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE,
)

ALLOW_CHARGE_SEPARATION: ResonanceFlags
ALLOW_INCOMPLETE_OCTETS: ResonanceFlags
AllProps: PropertyPickleOptions
AtomProps: PropertyPickleOptions
BondProps: PropertyPickleOptions
CHI_ALLENE: ChiralType
CHI_OCTAHEDRAL: ChiralType
CHI_OTHER: ChiralType
CHI_SQUAREPLANAR: ChiralType
CHI_TETRAHEDRAL: ChiralType
CHI_TETRAHEDRAL_CCW: ChiralType
CHI_TETRAHEDRAL_CW: ChiralType
CHI_TRIGONALBIPYRAMIDAL: ChiralType
CHI_UNSPECIFIED: ChiralType
COMPOSITE_AND: CompositeQueryType
COMPOSITE_OR: CompositeQueryType
COMPOSITE_XOR: CompositeQueryType
ComputedProps: PropertyPickleOptions
CoordsAsDouble: PropertyPickleOptions
KEKULE_ALL: ResonanceFlags
MolProps: PropertyPickleOptions
NoConformers: PropertyPickleOptions
NoProps: PropertyPickleOptions
PrivateProps: PropertyPickleOptions
QueryAtomData: PropertyPickleOptions
STEREO_ABSOLUTE: StereoGroupType
STEREO_AND: StereoGroupType
STEREO_OR: StereoGroupType
UNCONSTRAINED_ANIONS: ResonanceFlags
UNCONSTRAINED_CATIONS: ResonanceFlags

class Atom(Boost.Python.instance):
    """
    The class to store Atoms.
    Note that, though it is possible to create one, having an Atom on its own
    (i.e not associated with a molecule) is not particularly useful.

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)

    __init__( (object)arg1, (Atom)arg2) -> None :

    C++ signature :void __init__(_object*,RDKit::Atom)

    __init__( (object)arg1, (int)arg2) -> None :Constructor, takes either an int (atomic number) or a string (atomic symbol).

    C++ signature :void __init__(_object*,unsigned int)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def ClearProp(self, arg1: Atom, arg2: str) -> None:
        """
        Removes a particular property from an Atom (does nothing if not already set).

        ARGUMENTS:
        key: the name of the property to be removed.

        C++ signature :void ClearProp(RDKit::Atom const*,char const*)"""
        ...
    def DescribeQuery(self, arg1: Atom) -> str:
        """
        returns a text description of the query. Primarily intended for debugging purposes.

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > DescribeQuery(RDKit::Atom const*)
        """
        ...
    def GetAtomMapNum(self, arg1: Atom) -> int:
        """
        Gets the atoms map number, returns 0 if not set

        C++ signature :int GetAtomMapNum(RDKit::Atom {lvalue})"""
        ...
    def GetAtomicNum(self, arg1: Atom) -> int:
        """
        Returns the atomic number.

        C++ signature :int GetAtomicNum(RDKit::Atom {lvalue})"""
        ...
    def GetBonds(self, arg1: Atom) -> tuple:
        """
        Returns a read-only sequence of the atom’s bonds

        C++ signature :boost::python::tuple GetBonds(RDKit::Atom*)"""
        ...
    def GetBoolProp(self, arg1: Atom, arg2: str) -> bool:
        """
        Returns the value of the property.

        ARGUMENTS:
        key: the name of the property to return (a bool).

        RETURNS: a bool

        NOTE:
        If the property has not been set, a KeyError exception will be raised.

        C++ signature :bool GetBoolProp(RDKit::Atom*,char const*)"""
        ...
    def GetChiralTag(self, arg1: Atom) -> ChiralType:
        """
        C++ signature :RDKit::Atom::ChiralType GetChiralTag(RDKit::Atom {lvalue})"""
        ...
    def GetDegree(self, arg1: Atom) -> int:
        """
        Returns the degree of the atom in the molecule.

        The degree of an atom is defined to be its number of
        directly-bonded neighbors.
        The degree is independent of bond orders, but is dependent

        on whether or not Hs are explicit in the graph.

        C++ signature :unsigned int GetDegree(RDKit::Atom {lvalue})"""
        ...
    def GetDoubleProp(self, arg1: Atom, arg2: str) -> float:
        """
        Returns the value of the property.

        ARGUMENTS:
        key: the name of the property to return (a double).

        RETURNS: a double

        NOTE:
        If the property has not been set, a KeyError exception will be raised.

        C++ signature :double GetDoubleProp(RDKit::Atom*,char const*)"""
        ...
    def GetExplicitBitVectProp(self, arg1: Atom, arg2: str) -> ExplicitBitVect:
        """
        Returns the value of the property.

        ARGUMENTS:
        key: the name of the property to return (a ExplicitBitVect).

        RETURNS: an ExplicitBitVect

        NOTE:
        If the property has not been set, a KeyError exception will be raised.

        C++ signature :ExplicitBitVect GetExplicitBitVectProp(RDKit::Atom*,char const*)
        """
        ...
    def GetExplicitValence(self, arg1: Atom) -> int:
        """
        Returns the explicit valence of the atom.

        C++ signature :int GetExplicitValence(RDKit::Atom {lvalue})"""
        ...
    def GetFormalCharge(self, arg1: Atom) -> int:
        """
        C++ signature :int GetFormalCharge(RDKit::Atom {lvalue})"""
        ...
    def GetHybridization(self, arg1: Atom) -> HybridizationType:
        """
        Returns the atom’s hybridization.

        C++ signature :RDKit::Atom::HybridizationType GetHybridization(RDKit::Atom {lvalue})
        """
        ...
    def GetIdx(self, arg1: Atom) -> int:
        """
        Returns the atom’s index (ordering in the molecule)

        C++ signature :unsigned int GetIdx(RDKit::Atom {lvalue})"""
        ...
    def GetImplicitValence(self, arg1: Atom) -> int:
        """
        Returns the number of implicit Hs on the atom.

        C++ signature :int GetImplicitValence(RDKit::Atom {lvalue})"""
        ...
    def GetIntProp(self, arg1: Atom, arg2: str) -> int:
        """
        Returns the value of the property.

        ARGUMENTS:
        key: the name of the property to return (an int).

        RETURNS: an int

        NOTE:
        If the property has not been set, a KeyError exception will be raised.

        C++ signature :int GetIntProp(RDKit::Atom*,char const*)"""
        ...
    def GetIsAromatic(self, arg1: Atom) -> bool:
        """
        C++ signature :bool GetIsAromatic(RDKit::Atom {lvalue})"""
        ...
    def GetIsotope(self, arg1: Atom) -> int:
        """
        C++ signature :unsigned int GetIsotope(RDKit::Atom {lvalue})"""
        ...
    def GetMass(self, arg1: Atom) -> float:
        """
        C++ signature :double GetMass(RDKit::Atom {lvalue})"""
        ...
    def GetMonomerInfo(self, arg1: Atom) -> AtomMonomerInfo:
        """
        Returns the atom’s MonomerInfo object, if there is one.

        C++ signature :RDKit::AtomMonomerInfo* GetMonomerInfo(RDKit::Atom*)"""
        ...
    def GetNeighbors(self, arg1: Atom) -> tuple:
        """
        Returns a read-only sequence of the atom’s neighbors

        C++ signature :boost::python::tuple GetNeighbors(RDKit::Atom*)"""
        ...
    def GetNoImplicit(self, arg1: Atom) -> bool:
        """
        Returns whether or not the atom is allowed to have implicit Hs.

        C++ signature :bool GetNoImplicit(RDKit::Atom {lvalue})"""
        ...
    def GetNumExplicitHs(self, arg1: Atom) -> int:
        """
        C++ signature :unsigned int GetNumExplicitHs(RDKit::Atom {lvalue})"""
        ...
    def GetNumImplicitHs(self, arg1: Atom) -> int:
        """
        Returns the total number of implicit Hs on the atom.

        C++ signature :unsigned int GetNumImplicitHs(RDKit::Atom {lvalue})"""
        ...
    def GetNumRadicalElectrons(self, arg1: Atom) -> int:
        """
        C++ signature :unsigned int GetNumRadicalElectrons(RDKit::Atom {lvalue})"""
        ...
    def GetOwningMol(self, arg1: Atom) -> Mol:
        """
        Returns the Mol that owns this atom.

        C++ signature :RDKit::ROMol {lvalue} GetOwningMol(RDKit::Atom {lvalue})"""
        ...
    def GetPDBResidueInfo(self, arg1: Atom) -> AtomPDBResidueInfo:
        """
        Returns the atom’s MonomerInfo object, if there is one.

        C++ signature :RDKit::AtomPDBResidueInfo* GetPDBResidueInfo(RDKit::Atom*)"""
        ...
    def GetProp(self, arg1: Atom, arg2: str) -> str:
        """
        Returns the value of the property.

        ARGUMENTS:
        key: the name of the property to return (a string).

        RETURNS: a string

        NOTE:
        If the property has not been set, a KeyError exception will be raised.

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetProp(RDKit::Atom*,char const*)
        """
        ...
    def GetPropNames(
        self: Atom, includePrivate: bool = False, includeComputed: bool = False
    ) -> _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE:
        """
        Returns a list of the properties set on the Atom.

        C++ signature :std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > GetPropNames(RDKit::Atom {lvalue} [,bool=False [,bool=False]])
        """
        ...
    def GetPropsAsDict(
        self: Atom, includePrivate: bool = True, includeComputed: bool = True
    ) -> dict:
        """
        Returns a dictionary of the properties set on the Atom.n.b. some properties cannot be converted to python types.

        C++ signature :boost::python::dict GetPropsAsDict(RDKit::Atom [,bool=True [,bool=True]])
        """
        ...
    def GetQueryType(self, arg1: Atom) -> str:
        """
        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetQueryType(RDKit::Atom {lvalue})
        """
        ...
    def GetSmarts(
        self: Atom,
        doKekule: bool = False,
        allHsExplicit: bool = False,
        isomericSmiles: bool = True,
    ) -> str:
        """
        returns the SMARTS (or SMILES) string for an Atom

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetSmarts(RDKit::Atom const* [,bool=False [,bool=False [,bool=True]]])
        """
        ...
    def GetSymbol(self, arg1: Atom) -> str:
        """
        Returns the atomic symbol (a string)

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetSymbol(RDKit::Atom {lvalue})
        """
        ...
    def GetTotalDegree(self, arg1: Atom) -> int:
        """
        Returns the degree of the atom in the molecule including Hs.

        The degree of an atom is defined to be its number of
        directly-bonded neighbors.
        The degree is independent of bond orders.

        C++ signature :unsigned int GetTotalDegree(RDKit::Atom {lvalue})"""
        ...
    def GetTotalNumHs(self: Atom, includeNeighbors: bool = False) -> int:
        """
        Returns the total number of Hs (explicit and implicit) on the atom.

        ARGUMENTS:

        includeNeighbors: (optional) toggles inclusion of neighboring H atoms in the sum.
        Defaults to 0.

        C++ signature :unsigned int GetTotalNumHs(RDKit::Atom {lvalue} [,bool=False])"""
        ...
    def GetTotalValence(self, arg1: Atom) -> int:
        """
        Returns the total valence (explicit + implicit) of the atom.

        C++ signature :unsigned int GetTotalValence(RDKit::Atom {lvalue})"""
        ...
    def GetUnsignedProp(self, arg1: Atom, arg2: str) -> int:
        """
        Returns the value of the property.

        ARGUMENTS:
        key: the name of the property to return (an unsigned integer).

        RETURNS: an integer (Python has no unsigned type)

        NOTE:
        If the property has not been set, a KeyError exception will be raised.

        C++ signature :unsigned int GetUnsignedProp(RDKit::Atom*,char const*)"""
        ...
    def HasOwningMol(self, arg1: Atom) -> bool:
        """
        Returns whether or not this instance belongs to a molecule.

        C++ signature :bool HasOwningMol(RDKit::Atom {lvalue})"""
        ...
    def HasProp(self, arg1: Atom, arg2: str) -> int:
        """
        Queries a Atom to see if a particular property has been assigned.

        ARGUMENTS:
        key: the name of the property to check for (a string).

        C++ signature :int HasProp(RDKit::Atom const*,char const*)"""
        ...
    def HasQuery(self, arg1: Atom) -> bool:
        """
        Returns whether or not the atom has an associated query

        C++ signature :bool HasQuery(RDKit::Atom {lvalue})"""
        ...
    def InvertChirality(self, arg1: Atom) -> bool:
        """
        C++ signature :bool InvertChirality(RDKit::Atom {lvalue})"""
        ...
    def IsInRing(self, arg1: Atom) -> bool:
        """
        Returns whether or not the atom is in a ring

        C++ signature :bool IsInRing(RDKit::Atom const*)"""
        ...
    def IsInRingSize(self, arg1: Atom, arg2: int) -> bool:
        """
        Returns whether or not the atom is in a ring of a particular size.

        ARGUMENTS:
        size: the ring size to look for

        C++ signature :bool IsInRingSize(RDKit::Atom const*,int)"""
        ...
    def Match(self, arg1: Atom, arg2: Atom) -> bool:
        """
        Returns whether or not this atom matches another Atom.

        Each Atom (or query Atom) has a query function which is
        used for this type of matching.

        ARGUMENTS:
        other: the other Atom to which to compare

        C++ signature :bool Match(RDKit::Atom {lvalue},RDKit::Atom const*)"""
        ...
    def NeedsUpdatePropertyCache(self: Atom) -> bool:
        """
        Returns true or false depending on whether implicit and explicit valence of the molecule have already been calculated.

        C++ signature :bool NeedsUpdatePropertyCache(RDKit::Atom {lvalue})"""
        ...
    def SetAtomMapNum(self: Atom, mapno: int, strict: bool = False) -> None:
        """
        Sets the atoms map number, a value of 0 clears the atom map

        C++ signature :void SetAtomMapNum(RDKit::Atom {lvalue},int [,bool=False])"""
        ...
    def SetAtomicNum(self, arg1: Atom, arg2: int) -> None:
        """
        Sets the atomic number, takes an integer value as an argument

        C++ signature :void SetAtomicNum(RDKit::Atom {lvalue},int)"""
        ...
    def SetBoolProp(self: Atom, key: str, val: bool) -> None:
        """
        Sets an atomic property

        ARGUMENTS:
        key: the name of the property to be set (a bool).
        value: the property value (a bool).

        C++ signature :void SetBoolProp(RDKit::Atom const*,char const*,bool)"""
        ...
    def SetChiralTag(self, arg1: Atom, arg2: ChiralType) -> None:
        """
        C++ signature :void SetChiralTag(RDKit::Atom {lvalue},RDKit::Atom::ChiralType)
        """
        ...
    def SetDoubleProp(self: Atom, key: str, val: float) -> None:
        """
        Sets an atomic property

        ARGUMENTS:
        key: the name of the property to be set (a double).
        value: the property value (a double).

        C++ signature :void SetDoubleProp(RDKit::Atom const*,char const*,double)"""
        ...
    def SetExplicitBitVectProp(self: Atom, key: str, val: ExplicitBitVect) -> None:
        """
        Sets an atomic property

        ARGUMENTS:
        key: the name of the property to be set (an ExplicitBitVect).
        value: the property value (an ExplicitBitVect).

        C++ signature :void SetExplicitBitVectProp(RDKit::Atom const*,char const*,ExplicitBitVect)
        """
        ...
    def SetFormalCharge(self, arg1: Atom, arg2: int) -> None:
        """
        C++ signature :void SetFormalCharge(RDKit::Atom {lvalue},int)"""
        ...
    def SetHybridization(self, arg1: Atom, arg2: HybridizationType) -> None:
        """
        Sets the hybridization of the atom.The argument should be a HybridizationType

        C++ signature :void SetHybridization(RDKit::Atom {lvalue},RDKit::Atom::HybridizationType)
        """
        ...
    def SetIntProp(self: Atom, key: str, val: int) -> None:
        """
        Sets an atomic property

        ARGUMENTS:
        key: the name of the property to be set (a int).
        value: the property value (a int).

        C++ signature :void SetIntProp(RDKit::Atom const*,char const*,int)"""
        ...
    def SetIsAromatic(self, arg1: Atom, arg2: bool) -> None:
        """
        C++ signature :void SetIsAromatic(RDKit::Atom {lvalue},bool)"""
        ...
    def SetIsotope(self, arg1: Atom, arg2: int) -> None:
        """
        C++ signature :void SetIsotope(RDKit::Atom {lvalue},unsigned int)"""
        ...
    def SetMonomerInfo(self, arg1: Atom, arg2: AtomMonomerInfo) -> None:
        """
        Sets the atom’s MonomerInfo object.

        C++ signature :void SetMonomerInfo(RDKit::Atom*,RDKit::AtomMonomerInfo const*)
        """
        ...
    def SetNoImplicit(self, arg1: Atom, arg2: bool) -> None:
        """
        Sets a marker on the atom that disallows implicit Hs.This holds even if the atom would otherwise have implicit Hs added.

        C++ signature :void SetNoImplicit(RDKit::Atom {lvalue},bool)"""
        ...
    def SetNumExplicitHs(self, arg1: Atom, arg2: int) -> None:
        """
        C++ signature :void SetNumExplicitHs(RDKit::Atom {lvalue},unsigned int)"""
        ...
    def SetNumRadicalElectrons(self, arg1: Atom, arg2: int) -> None:
        """
        C++ signature :void SetNumRadicalElectrons(RDKit::Atom {lvalue},unsigned int)"""
        ...
    def SetPDBResidueInfo(self, arg1: Atom, arg2: AtomMonomerInfo) -> None:
        """
        Sets the atom’s MonomerInfo object.

        C++ signature :void SetPDBResidueInfo(RDKit::Atom*,RDKit::AtomMonomerInfo const*)
        """
        ...
    def SetProp(self: Atom, key: str, val: str) -> None:
        """
        Sets an atomic property

        ARGUMENTS:
        key: the name of the property to be set (a string).
        value: the property value (a string).

        C++ signature :void SetProp(RDKit::Atom const*,char const*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    def SetUnsignedProp(self: Atom, key: str, val: int) -> None:
        """
        Sets an atomic property

        ARGUMENTS:
        key: the name of the property to be set (an unsigned integer).
        value: the property value (a int >= 0).

        C++ signature :void SetUnsignedProp(RDKit::Atom const*,char const*,unsigned int)
        """
        ...
    def UpdatePropertyCache(self: Atom, strict: bool = True) -> None:
        """
        Regenerates computed properties like implicit valence and ring information.

        C++ signature :void UpdatePropertyCache(RDKit::Atom {lvalue} [,bool=True])"""
        ...
    @classmethod
    def __copy__(cls, RDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class AtomKekulizeException(AtomSanitizeException):
    ...
    ...

class AtomMonomerInfo(Boost.Python.instance):
    """
    The class to store monomer information attached to Atoms

    C++ signature :void __init__(_object*)

    __init__( (object)arg1, (AtomMonomerType)type [, (str)name=’’]) -> None :

    C++ signature :void __init__(_object*,RDKit::AtomMonomerInfo::AtomMonomerType [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’’])
    """

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def GetMonomerType(self, arg1: AtomMonomerInfo) -> AtomMonomerType:
        """
        C++ signature :RDKit::AtomMonomerInfo::AtomMonomerType GetMonomerType(RDKit::AtomMonomerInfo {lvalue})
        """
        ...
    def GetName(self, arg1: AtomMonomerInfo) -> str:
        """
        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetName(RDKit::AtomMonomerInfo {lvalue})
        """
        ...
    def SetMonomerType(self, arg1: AtomMonomerInfo, arg2: AtomMonomerType) -> None:
        """
        C++ signature :void SetMonomerType(RDKit::AtomMonomerInfo {lvalue},RDKit::AtomMonomerInfo::AtomMonomerType)
        """
        ...
    def SetName(self, arg1: AtomMonomerInfo, arg2: str) -> None:
        """
        C++ signature :void SetName(RDKit::AtomMonomerInfo {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class AtomMonomerType(Boost.Python.enum):
    OTHER: AtomMonomerType = ...
    PDBRESIDUE: AtomMonomerType = ...
    UNKNOWN: AtomMonomerType = ...
    names: dict[str, AtomMonomerType] = ...
    values: dict[int, AtomMonomerType] = ...
    __slots__: ClassVar[tuple] = ...

class AtomPDBResidueInfo(AtomMonomerInfo):
    """
    The class to store PDB residue information attached to Atoms

    C++ signature :void __init__(_object*)

    __init__( (object)arg1, (str)atomName [, (int)serialNumber=1 [, (str)altLoc=’’ [, (str)residueName=’’ [, (int)residueNumber=0 [, (str)chainId=’’ [, (str)insertionCode=’’ [, (float)occupancy=1.0 [, (float)tempFactor=0.0 [, (bool)isHeteroAtom=False [, (int)secondaryStructure=0 [, (int)segmentNumber=0]]]]]]]]]]]) -> None :

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,int=1 [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’’ [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’’ [,int=0 [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’’ [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’’ [,double=1.0 [,double=0.0 [,bool=False [,unsigned int=0 [,unsigned int=0]]]]]]]]]]])
    """

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def GetAltLoc(self, arg1: AtomPDBResidueInfo) -> str:
        """
        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetAltLoc(RDKit::AtomPDBResidueInfo {lvalue})
        """
        ...
    def GetChainId(self, arg1: AtomPDBResidueInfo) -> str:
        """
        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetChainId(RDKit::AtomPDBResidueInfo {lvalue})
        """
        ...
    def GetInsertionCode(self, arg1: AtomPDBResidueInfo) -> str:
        """
        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetInsertionCode(RDKit::AtomPDBResidueInfo {lvalue})
        """
        ...
    def GetIsHeteroAtom(self, arg1: AtomPDBResidueInfo) -> bool:
        """
        C++ signature :bool GetIsHeteroAtom(RDKit::AtomPDBResidueInfo {lvalue})"""
        ...
    def GetOccupancy(self, arg1: AtomPDBResidueInfo) -> float:
        """
        C++ signature :double GetOccupancy(RDKit::AtomPDBResidueInfo {lvalue})"""
        ...
    def GetResidueName(self, arg1: AtomPDBResidueInfo) -> str:
        """
        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetResidueName(RDKit::AtomPDBResidueInfo {lvalue})
        """
        ...
    def GetResidueNumber(self, arg1: AtomPDBResidueInfo) -> int:
        """
        C++ signature :int GetResidueNumber(RDKit::AtomPDBResidueInfo {lvalue})"""
        ...
    def GetSecondaryStructure(self, arg1: AtomPDBResidueInfo) -> int:
        """
        C++ signature :unsigned int GetSecondaryStructure(RDKit::AtomPDBResidueInfo {lvalue})
        """
        ...
    def GetSegmentNumber(self, arg1: AtomPDBResidueInfo) -> int:
        """
        C++ signature :unsigned int GetSegmentNumber(RDKit::AtomPDBResidueInfo {lvalue})
        """
        ...
    def GetSerialNumber(self, arg1: AtomPDBResidueInfo) -> int:
        """
        C++ signature :int GetSerialNumber(RDKit::AtomPDBResidueInfo {lvalue})"""
        ...
    def GetTempFactor(self, arg1: AtomPDBResidueInfo) -> float:
        """
        C++ signature :double GetTempFactor(RDKit::AtomPDBResidueInfo {lvalue})"""
        ...
    def SetAltLoc(self, arg1: AtomPDBResidueInfo, arg2: str) -> None:
        """
        C++ signature :void SetAltLoc(RDKit::AtomPDBResidueInfo {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    def SetChainId(self, arg1: AtomPDBResidueInfo, arg2: str) -> None:
        """
        C++ signature :void SetChainId(RDKit::AtomPDBResidueInfo {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    def SetInsertionCode(self, arg1: AtomPDBResidueInfo, arg2: str) -> None:
        """
        C++ signature :void SetInsertionCode(RDKit::AtomPDBResidueInfo {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    def SetIsHeteroAtom(self, arg1: AtomPDBResidueInfo, arg2: bool) -> None:
        """
        C++ signature :void SetIsHeteroAtom(RDKit::AtomPDBResidueInfo {lvalue},bool)"""
        ...
    def SetOccupancy(self, arg1: AtomPDBResidueInfo, arg2: float) -> None:
        """
        C++ signature :void SetOccupancy(RDKit::AtomPDBResidueInfo {lvalue},double)"""
        ...
    def SetResidueName(self, arg1: AtomPDBResidueInfo, arg2: str) -> None:
        """
        C++ signature :void SetResidueName(RDKit::AtomPDBResidueInfo {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    def SetResidueNumber(self, arg1: AtomPDBResidueInfo, arg2: int) -> None:
        """
        C++ signature :void SetResidueNumber(RDKit::AtomPDBResidueInfo {lvalue},int)"""
        ...
    def SetSecondaryStructure(self, arg1: AtomPDBResidueInfo, arg2: int) -> None:
        """
        C++ signature :void SetSecondaryStructure(RDKit::AtomPDBResidueInfo {lvalue},unsigned int)
        """
        ...
    def SetSegmentNumber(self, arg1: AtomPDBResidueInfo, arg2: int) -> None:
        """
        C++ signature :void SetSegmentNumber(RDKit::AtomPDBResidueInfo {lvalue},unsigned int)
        """
        ...
    def SetSerialNumber(self, arg1: AtomPDBResidueInfo, arg2: int) -> None:
        """
        C++ signature :void SetSerialNumber(RDKit::AtomPDBResidueInfo {lvalue},int)"""
        ...
    def SetTempFactor(self, arg1: AtomPDBResidueInfo, arg2: float) -> None:
        """
        C++ signature :void SetTempFactor(RDKit::AtomPDBResidueInfo {lvalue},double)"""
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class AtomSanitizeException(MolSanitizeException):
    ...
    ...

class AtomValenceException(AtomSanitizeException):
    ...
    ...

class Bond(Boost.Python.instance):
    """
    The class to store Bonds.
    Note: unlike Atoms, is it currently impossible to construct Bonds from
    Python.
    Raises an exception
    This class cannot be instantiated from Python"""

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def ClearProp(self, arg1: Bond, arg2: str) -> None:
        """
        Removes a particular property from an Bond (does nothing if not already set).

        ARGUMENTS:
        key: the name of the property to be removed.

        C++ signature :void ClearProp(RDKit::Bond const*,char const*)"""
        ...
    def DescribeQuery(self, arg1: Bond) -> str:
        """
        returns a text description of the query. Primarily intended for debugging purposes.

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > DescribeQuery(RDKit::Bond const*)
        """
        ...
    def GetBeginAtom(self, arg1: Bond) -> Atom:
        """
        Returns the bond’s first atom.

        C++ signature :RDKit::Atom* GetBeginAtom(RDKit::Bond {lvalue})"""
        ...
    def GetBeginAtomIdx(self, arg1: Bond) -> int:
        """
        Returns the index of the bond’s first atom.

        C++ signature :unsigned int GetBeginAtomIdx(RDKit::Bond {lvalue})"""
        ...
    def GetBondDir(self, arg1: Bond) -> BondDir:
        """
        Returns the type of the bond as a BondDir

        C++ signature :RDKit::Bond::BondDir GetBondDir(RDKit::Bond {lvalue})"""
        ...
    def GetBondType(self, arg1: Bond) -> BondType:
        """
        Returns the type of the bond as a BondType

        C++ signature :RDKit::Bond::BondType GetBondType(RDKit::Bond {lvalue})"""
        ...
    def GetBondTypeAsDouble(self, arg1: Bond) -> float:
        """
        Returns the type of the bond as a double (i.e. 1.0 for SINGLE, 1.5 for AROMATIC, 2.0 for DOUBLE)

        C++ signature :double GetBondTypeAsDouble(RDKit::Bond {lvalue})"""
        ...
    def GetBoolProp(self, arg1: Bond, arg2: str) -> bool:
        """
        Returns the value of the property.

        ARGUMENTS:
        key: the name of the property to return (a boolean).

        RETURNS: a boolean

        NOTE:
        If the property has not been set, a KeyError exception will be raised.

        C++ signature :bool GetBoolProp(RDKit::Bond*,char const*)"""
        ...
    def GetDoubleProp(self, arg1: Bond, arg2: str) -> float:
        """
        Returns the value of the property.

        ARGUMENTS:
        key: the name of the property to return (a double).

        RETURNS: a double

        NOTE:
        If the property has not been set, a KeyError exception will be raised.

        C++ signature :double GetDoubleProp(RDKit::Bond*,char const*)"""
        ...
    def GetEndAtom(self, arg1: Bond) -> Atom:
        """
        Returns the bond’s second atom.

        C++ signature :RDKit::Atom* GetEndAtom(RDKit::Bond {lvalue})"""
        ...
    def GetEndAtomIdx(self, arg1: Bond) -> int:
        """
        Returns the index of the bond’s first atom.

        C++ signature :unsigned int GetEndAtomIdx(RDKit::Bond {lvalue})"""
        ...
    def GetIdx(self, arg1: Bond) -> int:
        """
        Returns the bond’s index (ordering in the molecule)

        C++ signature :unsigned int GetIdx(RDKit::Bond {lvalue})"""
        ...
    def GetIntProp(self, arg1: Bond, arg2: str) -> int:
        """
        Returns the value of the property.

        ARGUMENTS:
        key: the name of the property to return (an int).

        RETURNS: an int

        NOTE:
        If the property has not been set, a KeyError exception will be raised.

        C++ signature :int GetIntProp(RDKit::Bond*,char const*)"""
        ...
    def GetIsAromatic(self, arg1: Bond) -> bool:
        """
        C++ signature :bool GetIsAromatic(RDKit::Bond {lvalue})"""
        ...
    def GetIsConjugated(self, arg1: Bond) -> bool:
        """
        Returns whether or not the bond is considered to be conjugated.

        C++ signature :bool GetIsConjugated(RDKit::Bond {lvalue})"""
        ...
    def GetOtherAtom(self, arg1: Bond, arg2: Atom) -> Atom:
        """
        Given one of the bond’s atoms, returns the other one.

        C++ signature :RDKit::Atom* GetOtherAtom(RDKit::Bond {lvalue},RDKit::Atom const*)
        """
        ...
    def GetOtherAtomIdx(self, arg1: Bond, arg2: int) -> int:
        """
        Given the index of one of the bond’s atoms, returns the
        index of the other.

        C++ signature :unsigned int GetOtherAtomIdx(RDKit::Bond {lvalue},unsigned int)
        """
        ...
    def GetOwningMol(self, arg1: Bond) -> Mol:
        """
        Returns the Mol that owns this bond.

        C++ signature :RDKit::ROMol {lvalue} GetOwningMol(RDKit::Bond {lvalue})"""
        ...
    def GetProp(self, arg1: Bond, arg2: str) -> str:
        """
        Returns the value of the property.

        ARGUMENTS:
        key: the name of the property to return (a string).

        RETURNS: a string

        NOTE:
        If the property has not been set, a KeyError exception will be raised.

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetProp(RDKit::Bond*,char const*)
        """
        ...
    def GetPropNames(
        self: Bond, includePrivate: bool = False, includeComputed: bool = False
    ) -> _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE:
        """
        Returns a list of the properties set on the Bond.

        C++ signature :std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > GetPropNames(RDKit::Bond {lvalue} [,bool=False [,bool=False]])
        """
        ...
    def GetPropsAsDict(
        self: Bond, includePrivate: bool = True, includeComputed: bool = True
    ) -> dict:
        """
        Returns a dictionary of the properties set on the Bond.n.b. some properties cannot be converted to python types.

        C++ signature :boost::python::dict GetPropsAsDict(RDKit::Bond [,bool=True [,bool=True]])
        """
        ...
    def GetSmarts(self, bond: Bond, allBondsExplicit: bool = False) -> str:
        """
        returns the SMARTS (or SMILES) string for a Bond

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetSmarts(RDKit::Bond const* [,bool=False])
        """
        ...
    def GetStereo(self, arg1: Bond) -> BondStereo:
        """
        Returns the stereo configuration of the bond as a BondStereo

        C++ signature :RDKit::Bond::BondStereo GetStereo(RDKit::Bond {lvalue})"""
        ...
    def GetStereoAtoms(self, arg1: Bond) -> _vecti:
        """
        Returns the indices of the atoms setting this bond’s stereochemistry.

        C++ signature :std::vector<int, std::allocator<int> > GetStereoAtoms(RDKit::Bond const*)
        """
        ...
    def GetUnsignedProp(self, arg1: Bond, arg2: str) -> int:
        """
        Returns the value of the property.

        ARGUMENTS:
        key: the name of the property to return (an unsigned integer).

        RETURNS: an int (Python has no unsigned type)

        NOTE:
        If the property has not been set, a KeyError exception will be raised.

        C++ signature :unsigned int GetUnsignedProp(RDKit::Bond*,char const*)"""
        ...
    def GetValenceContrib(self, arg1: Bond, arg2: Atom) -> float:
        """
        Returns the contribution of the bond to the valence of an Atom.

        ARGUMENTS:

        atom: the Atom to consider.

        C++ signature :double GetValenceContrib(RDKit::Bond {lvalue},RDKit::Atom const*)
        """
        ...
    def HasOwningMol(self, arg1: Bond) -> bool:
        """
        Returns whether or not this instance belongs to a molecule.

        C++ signature :bool HasOwningMol(RDKit::Bond {lvalue})"""
        ...
    def HasProp(self, arg1: Bond, arg2: str) -> int:
        """
        Queries a Bond to see if a particular property has been assigned.

        ARGUMENTS:
        key: the name of the property to check for (a string).

        C++ signature :int HasProp(RDKit::Bond const*,char const*)"""
        ...
    def HasQuery(self, arg1: Bond) -> bool:
        """
        Returns whether or not the bond has an associated query

        C++ signature :bool HasQuery(RDKit::Bond {lvalue})"""
        ...
    def IsInRing(self, arg1: Bond) -> bool:
        """
        Returns whether or not the bond is in a ring of any size.

        C++ signature :bool IsInRing(RDKit::Bond const*)"""
        ...
    def IsInRingSize(self, arg1: Bond, arg2: int) -> bool:
        """
        Returns whether or not the bond is in a ring of a particular size.

        ARGUMENTS:
        size: the ring size to look for

        C++ signature :bool IsInRingSize(RDKit::Bond const*,int)"""
        ...
    def Match(self, arg1: Bond, arg2: Bond) -> bool:
        """
        Returns whether or not this bond matches another Bond.

        Each Bond (or query Bond) has a query function which is
        used for this type of matching.

        ARGUMENTS:
        other: the other Bond to which to compare

        C++ signature :bool Match(RDKit::Bond {lvalue},RDKit::Bond const*)"""
        ...
    def SetBondDir(self, arg1: Bond, arg2: BondDir) -> None:
        """
        Set the type of the bond as a BondDir

        C++ signature :void SetBondDir(RDKit::Bond {lvalue},RDKit::Bond::BondDir)"""
        ...
    def SetBondType(self, arg1: Bond, arg2: BondType) -> None:
        """
        Set the type of the bond as a BondType

        C++ signature :void SetBondType(RDKit::Bond {lvalue},RDKit::Bond::BondType)"""
        ...
    def SetBoolProp(self: Bond, key: str, val: bool) -> None:
        """
        Sets a bond property

        ARGUMENTS:
        key: the name of the property to be set (a string).
        value: the property value (a boolean).

        C++ signature :void SetBoolProp(RDKit::Bond const*,char const*,bool)"""
        ...
    def SetDoubleProp(self: Bond, key: str, val: float) -> None:
        """
        Sets a bond property

        ARGUMENTS:
        key: the name of the property to be set (a string).
        value: the property value (a double).

        C++ signature :void SetDoubleProp(RDKit::Bond const*,char const*,double)"""
        ...
    def SetIntProp(self: Bond, key: str, val: int) -> None:
        """
        Sets a bond property

        ARGUMENTS:
        key: the name of the property to be set (a string).
        value: the property value (an int).

        C++ signature :void SetIntProp(RDKit::Bond const*,char const*,int)"""
        ...
    def SetIsAromatic(self, arg1: Bond, arg2: bool) -> None:
        """
        C++ signature :void SetIsAromatic(RDKit::Bond {lvalue},bool)"""
        ...
    def SetIsConjugated(self, arg1: Bond, arg2: bool) -> None:
        """
        C++ signature :void SetIsConjugated(RDKit::Bond {lvalue},bool)"""
        ...
    def SetProp(self: Bond, key: str, val: str) -> None:
        """
        Sets a bond property

        ARGUMENTS:
        key: the name of the property to be set (a string).
        value: the property value (a string).

        C++ signature :void SetProp(RDKit::Bond const*,char const*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    def SetStereo(self, arg1: Bond, arg2: BondStereo) -> None:
        """
        Set the stereo configuration of the bond as a BondStereo

        C++ signature :void SetStereo(RDKit::Bond {lvalue},RDKit::Bond::BondStereo)"""
        ...
    def SetStereoAtoms(self, arg1: Bond, arg2: int, arg3: int) -> None:
        """
        Set the indices of the atoms setting this bond’s stereochemistry.

        C++ signature :void SetStereoAtoms(RDKit::Bond {lvalue},unsigned int,unsigned int)
        """
        ...
    def SetUnsignedProp(self: Bond, key: str, val: int) -> None:
        """
        Sets a bond property

        ARGUMENTS:
        key: the name of the property to be set (a string).
        value: the property value (an int >= 0).

        C++ signature :void SetUnsignedProp(RDKit::Bond const*,char const*,unsigned int)
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class BondDir(Boost.Python.enum):
    BEGINDASH: BondDir = ...
    BEGINWEDGE: BondDir = ...
    EITHERDOUBLE: BondDir = ...
    ENDDOWNRIGHT: BondDir = ...
    ENDUPRIGHT: BondDir = ...
    NONE: BondDir = ...
    UNKNOWN: BondDir = ...
    names: dict[str, BondDir] = ...
    values: dict[int, BondDir] = ...
    __slots__: ClassVar[tuple] = ...

class BondStereo(Boost.Python.enum):
    STEREOANY: BondStereo = ...
    STEREOCIS: BondStereo = ...
    STEREOE: BondStereo = ...
    STEREONONE: BondStereo = ...
    STEREOTRANS: BondStereo = ...
    STEREOZ: BondStereo = ...
    names: dict[str, BondStereo] = ...
    values: dict[int, BondStereo] = ...
    __slots__: ClassVar[tuple] = ...

class BondType(Boost.Python.enum):
    AROMATIC: BondType = ...
    DATIVE: BondType = ...
    DATIVEL: BondType = ...
    DATIVEONE: BondType = ...
    DATIVER: BondType = ...
    DOUBLE: BondType = ...
    FIVEANDAHALF: BondType = ...
    FOURANDAHALF: BondType = ...
    HEXTUPLE: BondType = ...
    HYDROGEN: BondType = ...
    IONIC: BondType = ...
    ONEANDAHALF: BondType = ...
    OTHER: BondType = ...
    QUADRUPLE: BondType = ...
    QUINTUPLE: BondType = ...
    SINGLE: BondType = ...
    THREEANDAHALF: BondType = ...
    THREECENTER: BondType = ...
    TRIPLE: BondType = ...
    TWOANDAHALF: BondType = ...
    UNSPECIFIED: BondType = ...
    ZERO: BondType = ...
    names: dict[str, BondType] = ...
    values: dict[int, BondType] = ...
    __slots__: ClassVar[tuple] = ...

class ChiralType(Boost.Python.enum):
    CHI_ALLENE: ChiralType = ...
    CHI_OCTAHEDRAL: ChiralType = ...
    CHI_OTHER: ChiralType = ...
    CHI_SQUAREPLANAR: ChiralType = ...
    CHI_TETRAHEDRAL: ChiralType = ...
    CHI_TETRAHEDRAL_CCW: ChiralType = ...
    CHI_TETRAHEDRAL_CW: ChiralType = ...
    CHI_TRIGONALBIPYRAMIDAL: ChiralType = ...
    CHI_UNSPECIFIED: ChiralType = ...
    names: dict[str, ChiralType] = ...
    values: dict[int, ChiralType] = ...
    __slots__: ClassVar[tuple] = ...

class CompositeQueryType(Boost.Python.enum):
    COMPOSITE_AND: CompositeQueryType = ...
    COMPOSITE_OR: CompositeQueryType = ...
    COMPOSITE_XOR: CompositeQueryType = ...
    names: dict[str, CompositeQueryType] = ...
    values: dict[int, CompositeQueryType] = ...
    __slots__: ClassVar[tuple] = ...

class Conformer(Boost.Python.instance):
    """
    The class to store 2D or 3D conformation of a molecule

    C++ signature :void __init__(_object*)

    __init__( (object)arg1, (int)arg2) -> None :Constructor with the number of atoms specified

    C++ signature :void __init__(_object*,unsigned int)

    __init__( (object)arg1, (Conformer)arg2) -> None :

    C++ signature :void __init__(_object*,RDKit::Conformer)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def ClearComputedProps(self, arg1: Conformer) -> None:
        """
        Removes all computed properties from the conformer.

        C++ signature :void ClearComputedProps(RDKit::Conformer)"""
        ...
    def ClearProp(self, arg1: Conformer, arg2: str) -> None:
        """
        Removes a property from the conformer.

        ARGUMENTS:
        key: the name of the property to clear (a string).

        C++ signature :void ClearProp(RDKit::Conformer,char const*)"""
        ...
    def GetAtomPosition(self, arg1: Conformer, arg2: int) -> Point3D:
        """
        Get the posistion of an atom

        C++ signature :RDGeom::Point3D GetAtomPosition(RDKit::Conformer const*,unsigned int)
        """
        ...
    def GetBoolProp(self, arg1: Conformer, arg2: str) -> bool:
        """
        Returns the Bool value of the property if possible.

        ARGUMENTS:
        key: the name of the property to return (a string).

        RETURNS: a bool

        NOTE:
        If the property has not been set, a KeyError exception will be raised.

        C++ signature :bool GetBoolProp(RDKit::Conformer*,char const*)"""
        ...
    def GetDoubleProp(self, arg1: Conformer, arg2: str) -> float:
        """
        Returns the double value of the property if possible.

        ARGUMENTS:
        key: the name of the property to return (a string).

        RETURNS: a double

        NOTE:
        If the property has not been set, a KeyError exception will be raised.

        C++ signature :double GetDoubleProp(RDKit::Conformer*,char const*)"""
        ...
    def GetId(self, arg1: Conformer) -> int:
        """
        Get the ID of the conformer

        C++ signature :unsigned int GetId(RDKit::Conformer {lvalue})"""
        ...
    def GetIntProp(self, arg1: Conformer, arg2: str) -> int:
        """
        Returns the integer value of the property if possible.

        ARGUMENTS:
        key: the name of the property to return (a string).

        RETURNS: an integer

        NOTE:
        If the property has not been set, a KeyError exception will be raised.

        C++ signature :int GetIntProp(RDKit::Conformer*,char const*)"""
        ...
    def GetNumAtoms(self, arg1: Conformer) -> int:
        """
        Get the number of atoms in the conformer

        C++ signature :unsigned int GetNumAtoms(RDKit::Conformer {lvalue})"""
        ...
    def GetOwningMol(self, arg1: Conformer) -> Mol:
        """
        Get the owning molecule

        C++ signature :RDKit::ROMol {lvalue} GetOwningMol(RDKit::Conformer {lvalue})"""
        ...
    def GetPositions(self, arg1: Conformer) -> object:
        """
        Get positions of all the atoms

        C++ signature :_object* GetPositions(RDKit::Conformer const*)"""
        ...
    def GetProp(self, arg1: Conformer, arg2: str) -> str:
        """
        Returns the value of the property.

        ARGUMENTS:
        key: the name of the property to return (a string).

        RETURNS: a string

        NOTE:
        If the property has not been set, a KeyError exception will be raised.

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetProp(RDKit::Conformer*,char const*)
        """
        ...
    def GetPropNames(
        self: Conformer, includePrivate: bool = False, includeComputed: bool = False
    ) -> _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE:
        """
        Returns a tuple with all property names for this conformer.

        ARGUMENTS:

        includePrivate: (optional) toggles inclusion of private properties in the result set.Defaults to 0.

        includeComputed: (optional) toggles inclusion of computed properties in the result set.Defaults to 0.

        RETURNS: a tuple of strings

        C++ signature :std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > GetPropNames(RDKit::Conformer {lvalue} [,bool=False [,bool=False]])
        """
        ...
    def GetPropsAsDict(
        self: Conformer, includePrivate: bool = False, includeComputed: bool = False
    ) -> dict:
        """
        Returns a dictionary populated with the conformer’s properties.n.b. Some properties are not able to be converted to python types.

        ARGUMENTS:

        includePrivate: (optional) toggles inclusion of private properties in the result set.Defaults to False.

        includeComputed: (optional) toggles inclusion of computed properties in the result set.Defaults to False.

        RETURNS: a dictionary

        C++ signature :boost::python::dict GetPropsAsDict(RDKit::Conformer [,bool=False [,bool=False]])
        """
        ...
    def GetUnsignedProp(self, arg1: Conformer, arg2: str) -> int:
        """
        Returns the unsigned int value of the property if possible.

        ARGUMENTS:
        key: the name of the property to return (a string).

        RETURNS: an unsigned integer

        NOTE:
        If the property has not been set, a KeyError exception will be raised.

        C++ signature :unsigned int GetUnsignedProp(RDKit::Conformer*,char const*)"""
        ...
    def HasOwningMol(self, arg1: Conformer) -> bool:
        """
        Returns whether or not this instance belongs to a molecule.

        C++ signature :bool HasOwningMol(RDKit::Conformer {lvalue})"""
        ...
    def HasProp(self, arg1: Conformer, arg2: str) -> int:
        """
        Queries a conformer to see if a particular property has been assigned.

        ARGUMENTS:
        key: the name of the property to check for (a string).

        C++ signature :int HasProp(RDKit::Conformer,char const*)"""
        ...
    def Is3D(self, arg1: Conformer) -> bool:
        """
        returns the 3D flag of the conformer

        C++ signature :bool Is3D(RDKit::Conformer {lvalue})"""
        ...
    def Set3D(self, arg1: Conformer, arg2: bool) -> None:
        """
        Set the 3D flag of the conformer

        C++ signature :void Set3D(RDKit::Conformer {lvalue},bool)"""
        ...
    @overload
    def SetAtomPosition(
        self, arg1: Conformer, arg2: int, arg3: AtomPairsParameters
    ) -> None:
        """
        Set the position of the specified atom

        C++ signature :void SetAtomPosition(RDKit::Conformer {lvalue},unsigned int,RDGeom::Point3D)
        """
        ...
    @overload
    def SetAtomPosition(self, arg1: Conformer, arg2: int, arg3: Point3D) -> None: ...
    def SetBoolProp(
        self: Conformer, key: str, val: bool, computed: bool = False
    ) -> None:
        """
        Sets a boolean valued molecular property

        ARGUMENTS:
        key: the name of the property to be set (a string).
        value: the property value as a bool.

        computed: (optional) marks the property as being computed.Defaults to False.

        C++ signature :void SetBoolProp(RDKit::Conformer,char const*,bool [,bool=False])
        """
        ...
    def SetDoubleProp(
        self: Conformer, key: str, val: float, computed: bool = False
    ) -> None:
        """
        Sets a double valued molecular property

        ARGUMENTS:
        key: the name of the property to be set (a string).
        value: the property value as a double.

        computed: (optional) marks the property as being computed.Defaults to 0.

        C++ signature :void SetDoubleProp(RDKit::Conformer,char const*,double [,bool=False])
        """
        ...
    def SetId(self, arg1: Conformer, arg2: int) -> None:
        """
        Set the ID of the conformer

        C++ signature :void SetId(RDKit::Conformer {lvalue},unsigned int)"""
        ...
    def SetIntProp(self: Conformer, key: str, val: int, computed: bool = False) -> None:
        """
        Sets an integer valued molecular property

        ARGUMENTS:
        key: the name of the property to be set (an unsigned number).
        value: the property value as an integer.

        computed: (optional) marks the property as being computed.Defaults to False.

        C++ signature :void SetIntProp(RDKit::Conformer,char const*,int [,bool=False])
        """
        ...
    def SetProp(self: Conformer, key: str, val: str, computed: bool = False) -> None:
        """
        Sets a molecular property

        ARGUMENTS:
        key: the name of the property to be set (a string).
        value: the property value (a string).

        computed: (optional) marks the property as being computed.Defaults to False.

        C++ signature :void SetProp(RDKit::Conformer,char const*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])
        """
        ...
    def SetUnsignedProp(
        self: Conformer, key: str, val: int, computed: bool = False
    ) -> None:
        """
        Sets an unsigned integer valued molecular property

        ARGUMENTS:
        key: the name of the property to be set (a string).
        value: the property value as an unsigned integer.

        computed: (optional) marks the property as being computed.Defaults to False.

        C++ signature :void SetUnsignedProp(RDKit::Conformer,char const*,unsigned int [,bool=False])
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class EditableMol(Boost.Python.instance):
    """
    an editable molecule class
    Construct from a Mol

    C++ signature :void __init__(_object*,RDKit::ROMol)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def AddAtom(self, arg1: EditableMol, atom: Atom) -> int:
        """
        add an atom, returns the index of the newly added atom

        C++ signature :int AddAtom(RDKit::(anonymous namespace)::EditableMol {lvalue},RDKit::Atom*)
        """
        ...
    def AddBond(
        self,
        arg1: EditableMol,
        beginAtomIdx: int,
        endAtomIdx: int,
        order: BondType = BondType.UNSPECIFIED,
    ) -> int:
        """
        add a bond, returns the total number of bonds

        C++ signature :int AddBond(RDKit::(anonymous namespace)::EditableMol {lvalue},unsigned int,unsigned int [,RDKit::Bond::BondType=rdkit.Chem.rdchem.BondType.UNSPECIFIED])
        """
        ...
    def BeginBatchEdit(self, arg1: EditableMol) -> None:
        """
        starts batch editing

        C++ signature :void BeginBatchEdit(RDKit::(anonymous namespace)::EditableMol {lvalue})
        """
        ...
    def CommitBatchEdit(self, arg1: EditableMol) -> None:
        """
        finishes batch editing and makes the actual edits

        C++ signature :void CommitBatchEdit(RDKit::(anonymous namespace)::EditableMol {lvalue})
        """
        ...
    def GetMol(self, arg1: EditableMol) -> Mol:
        """
        Returns a Mol (a normal molecule)

        C++ signature :RDKit::ROMol* GetMol(RDKit::(anonymous namespace)::EditableMol {lvalue})
        """
        ...
    def RemoveAtom(self, arg1: EditableMol, arg2: int) -> None:
        """
        Remove the specified atom from the molecule

        C++ signature :void RemoveAtom(RDKit::(anonymous namespace)::EditableMol {lvalue},unsigned int)
        """
        ...
    def RemoveBond(self, arg1: EditableMol, arg2: int, arg3: int) -> None:
        """
        Remove the specified bond from the molecule

        C++ signature :void RemoveBond(RDKit::(anonymous namespace)::EditableMol {lvalue},unsigned int,unsigned int)
        """
        ...
    def ReplaceAtom(
        self,
        arg1: EditableMol,
        index: int,
        newAtom: Atom,
        updateLabel: bool = False,
        preserveProps: bool = False,
    ) -> None:
        """
        replaces the specified atom with the provided one
        If updateLabel is True, the new atom becomes the active atom
        If preserveProps is True preserve keep the existing props unless explicit set on the new atom

        C++ signature :void ReplaceAtom(RDKit::(anonymous namespace)::EditableMol {lvalue},unsigned int,RDKit::Atom* [,bool=False [,bool=False]])
        """
        ...
    def ReplaceBond(
        self, arg1: EditableMol, index: int, newBond: Bond, preserveProps: bool = False
    ) -> None:
        """
        replaces the specified bond with the provided one.
        If preserveProps is True preserve keep the existing props unless explicit set on the new bond

        C++ signature :void ReplaceBond(RDKit::(anonymous namespace)::EditableMol {lvalue},unsigned int,RDKit::Bond* [,bool=False])
        """
        ...
    def RollbackBatchEdit(self, arg1: EditableMol) -> None:
        """
        cancels batch editing

        C++ signature :void RollbackBatchEdit(RDKit::(anonymous namespace)::EditableMol {lvalue})
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class FixedMolSizeMolBundle(MolBundle):
    """
    A class for storing groups of related molecules.
    Here related means that the molecules have to have the same number of atoms.

    C++ signature :void __init__(_object*)"""

    ...
    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class HybridizationType(Boost.Python.enum):
    OTHER: HybridizationType = ...
    S: HybridizationType = ...
    SP: HybridizationType = ...
    SP2: HybridizationType = ...
    SP2D: HybridizationType = ...
    SP3: HybridizationType = ...
    SP3D: HybridizationType = ...
    SP3D2: HybridizationType = ...
    UNSPECIFIED: HybridizationType = ...
    names: dict[str, HybridizationType] = ...
    values: dict[int, HybridizationType] = ...
    __slots__: ClassVar[tuple] = ...

class KekulizeException(MolSanitizeException):
    ...
    ...

class Mol(Boost.Python.instance):
    """
    The Molecule class.

    In addition to the expected Atoms and Bonds, molecules contain:

    a collection of Atom and Bond bookmarks indexed with integersthat can be used to flag and retrieve particular Atoms or Bonds
    using the {get|set}{Atom|Bond}Bookmark() methods.

    a set of string-valued properties. These can have arbitrary stringlabels and can be set and retrieved using the {set|get}Prop() methods
    Molecular properties can be tagged as being computed, in which case

    they will be automatically cleared under certain circumstances (when the
    molecule itself is modified, for example).

    Molecules also have the concept of private properties, which are taggedby beginning the property name with an underscore (_).

    Constructor, takes no arguments

    C++ signature :void __init__(_object*)

    __init__( (object)arg1, (str)pklString) -> None :

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)

    __init__( (object)arg1, (str)pklString, (int)propertyFlags) -> None :

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,unsigned int)

    __init__( (object)arg1, (Mol)mol [, (bool)quickCopy=False [, (int)confId=-1]]) -> None :

    C++ signature :void __init__(_object*,RDKit::ROMol [,bool=False [,int=-1]])"""

    __getstate_manages_dict__: ClassVar[bool] = ...
    __instance_size__: ClassVar[int] = ...
    __safe_for_unpickling__: ClassVar[bool] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def AddConformer(self: Mol, conf: Conformer, assignId: bool = False) -> int:
        """
        Add a conformer to the molecule and return the conformer ID

        C++ signature :unsigned int AddConformer(RDKit::ROMol {lvalue},RDKit::Conformer* [,bool=False])
        """
        ...
    def ClearComputedProps(self: Mol, includeRings: bool = True) -> None:
        """
        Removes all computed properties from the molecule.

        C++ signature :void ClearComputedProps(RDKit::ROMol [,bool=True])"""
        ...
    def ClearProp(self, arg1: Mol, arg2: str) -> None:
        """
        Removes a property from the molecule.

        ARGUMENTS:
        key: the name of the property to clear (a string).

        C++ signature :void ClearProp(RDKit::ROMol,char const*)"""
        ...
    def Compute2DCoords(
        self,
        mol: Mol,
        canonOrient: bool = True,
        clearConfs: bool = True,
        coordMap: dict = {},
        nFlipsPerSample: int = 0,
        nSample: int = 0,
        sampleSeed: int = 0,
        permuteDeg4Nodes: bool = False,
        bondLength: float = -1.0,
        forceRDKit: bool = False,
        useRingTemplates: bool = False,
    ) -> int:
        """
        Compute 2D coordinates for a molecule. The resulting coordinates are stored on each atom of the molecule
        ARGUMENTS:

        mol - the molecule of interest
        canonOrient - orient the molecule in a canonical way
        clearConfs - if true, all existing conformations on the molecule

        will be cleared

        coordMap - a dictionary mapping atom Ids -> Point2D objects with starting coordinates for atoms that should
        have their positions locked.

        nFlipsPerSample - number of rotatable bonds that areflipped at random at a time.

        nSample - Number of random samplings of rotatable bonds.
        sampleSeed - seed for the random sampling process.
        permuteDeg4Nodes - allow permutation of bonds at a degree 4

        node during the sampling process

        bondLength - change the default bond length for depiction
        forceRDKit - use RDKit to generate coordinates even if

        preferCoordGen is set to true

        useRingTemplates - use templates to generate coordinates of complexring systems

        RETURNS:

        ID of the conformation added to the molecule

        C++ signature :unsigned int Compute2DCoords(RDKit::ROMol {lvalue} [,bool=True [,bool=True [,boost::python::dict {lvalue}={} [,unsigned int=0 [,unsigned int=0 [,int=0 [,bool=False [,double=-1.0 [,bool=False [,bool=False]]]]]]]]]])
        """
        ...
    def ComputeGasteigerCharges(
        self, mol: Mol, nIter: int = 12, throwOnParamFailure: bool = False
    ) -> None:
        """
        Compute Gasteiger partial charges for molecule

        The charges are computed using an iterative procedure presented in
        Ref : J.Gasteiger, M. Marseli, Iterative Equalization of Oribital Electronegatiity
        A Rapid Access to Atomic Charges, Tetrahedron Vol 36 p3219 1980
        The computed charges are stored on each atom are stored a computed property ( under the name
        _GasteigerCharge). In addition, each atom also stored the total charge for the implicit hydrogens
        on the atom (under the property name _GasteigerHCharge)
        ARGUMENTS:

        mol : the molecule of interrest
        nIter : number of iteration (defaults to 12)
        throwOnParamFailure : toggles whether or not an exception should be raised if parameters
        for an atom cannot be found.  If this is false (the default), all parameters for unknown
        atoms will be set to zero.  This has the effect of removing that atom from the iteration.

        C++ signature :void ComputeGasteigerCharges(RDKit::ROMol [,int=12 [,bool=False]])
        """
        ...
    def Debug(self, useStdout=False): ...
    def GetAromaticAtoms(self, arg1: Mol) -> _ROQAtomSeq:
        """
        Returns a read-only sequence containing all of the molecule’s aromatic Atoms.

        C++ signature :RDKit::ReadOnlySeq<RDKit::QueryAtomIterator_<RDKit::Atom, RDKit::ROMol>, RDKit::Atom*, RDKit::AtomCountFunctor>* GetAromaticAtoms(boost::shared_ptr<RDKit::ROMol>)
        """
        ...
    def GetAtomWithIdx(self, arg1: Mol, arg2: int) -> Atom:
        """
        Returns a particular Atom.

        ARGUMENTS:
        idx: which Atom to return

        NOTE: atom indices start at 0

        C++ signature :RDKit::Atom* GetAtomWithIdx(RDKit::ROMol {lvalue},unsigned int)
        """
        ...
    def GetAtoms(self, arg1: Mol) -> _ROAtomSeq:
        """
        Returns a read-only sequence containing all of the molecule’s Atoms.

        C++ signature :RDKit::ReadOnlySeq<RDKit::AtomIterator_<RDKit::Atom, RDKit::ROMol>, RDKit::Atom*, RDKit::AtomCountFunctor>* GetAtoms(boost::shared_ptr<RDKit::ROMol>)
        """
        ...
    def GetAtomsMatchingQuery(self, arg1: Mol, arg2: QueryAtom) -> _ROQAtomSeq:
        """
        Returns a read-only sequence containing all of the atoms in a molecule that match the query atom.

        C++ signature :RDKit::ReadOnlySeq<RDKit::QueryAtomIterator_<RDKit::Atom, RDKit::ROMol>, RDKit::Atom*, RDKit::AtomCountFunctor>* GetAtomsMatchingQuery(boost::shared_ptr<RDKit::ROMol>,RDKit::QueryAtom*)
        """
        ...
    def GetBondBetweenAtoms(self, arg1: Mol, arg2: int, arg3: int) -> Bond:
        """
        Returns the bond between two atoms, if there is one.

        ARGUMENTS:
        idx1,idx2: the Atom indices

        Returns:The Bond between the two atoms, if such a bond exists.
        If there is no Bond between the atoms, None is returned instead.

        NOTE: bond indices start at 0

        C++ signature :RDKit::Bond* GetBondBetweenAtoms(RDKit::ROMol {lvalue},unsigned int,unsigned int)
        """
        ...
    def GetBondWithIdx(self, arg1: Mol, arg2: int) -> Bond:
        """
        Returns a particular Bond.

        ARGUMENTS:
        idx: which Bond to return

        NOTE: bond indices start at 0

        C++ signature :RDKit::Bond* GetBondWithIdx(RDKit::ROMol {lvalue},unsigned int)
        """
        ...
    def GetBonds(self, arg1: Mol) -> _ROBondSeq:
        """
        Returns a read-only sequence containing all of the molecule’s Bonds.

        C++ signature :RDKit::ReadOnlySeq<RDKit::BondIterator_, RDKit::Bond*, RDKit::BondCountFunctor>* GetBonds(boost::shared_ptr<RDKit::ROMol>)
        """
        ...
    def GetBoolProp(self, arg1: Mol, arg2: str) -> bool:
        """
        Returns the Bool value of the property if possible.

        ARGUMENTS:
        key: the name of the property to return (a string).

        RETURNS: a bool

        NOTE:
        If the property has not been set, a KeyError exception will be raised.

        C++ signature :bool GetBoolProp(RDKit::ROMol*,char const*)"""
        ...
    def GetConformer(self: Mol, id: int = -1) -> Conformer:
        """
        Get the conformer with a specified ID

        C++ signature :RDKit::Conformer* GetConformer(RDKit::ROMol {lvalue} [,int=-1])
        """
        ...
    def GetConformers(self, arg1: Mol) -> _ROConformerSeq:
        """
        Returns a read-only sequence containing all of the molecule’s Conformers.

        C++ signature :RDKit::ReadOnlySeq<std::_List_iterator<boost::shared_ptr<RDKit::Conformer> >, boost::shared_ptr<RDKit::Conformer>&, RDKit::ConformerCountFunctor>* GetConformers(boost::shared_ptr<RDKit::ROMol>)
        """
        ...
    def GetDoubleProp(self, arg1: Mol, arg2: str) -> float:
        """
        Returns the double value of the property if possible.

        ARGUMENTS:
        key: the name of the property to return (a string).

        RETURNS: a double

        NOTE:
        If the property has not been set, a KeyError exception will be raised.

        C++ signature :double GetDoubleProp(RDKit::ROMol*,char const*)"""
        ...
    def GetIntProp(self, arg1: Mol, arg2: str) -> int:
        """
        Returns the integer value of the property if possible.

        ARGUMENTS:
        key: the name of the property to return (a string).

        RETURNS: an integer

        NOTE:
        If the property has not been set, a KeyError exception will be raised.

        C++ signature :int GetIntProp(RDKit::ROMol*,char const*)"""
        ...
    def GetNumAtoms(
        self, arg1: Mol, onlyHeavy: int = -1, onlyExplicit: bool = True
    ) -> int:
        """
        Returns the number of atoms in the molecule.

        ARGUMENTS:

        onlyExplicit: (optional) include only explicit atoms (atoms in the molecular graph)defaults to 1.

        NOTE: the onlyHeavy argument is deprecated

        C++ signature :int GetNumAtoms(RDKit::ROMol [,int=-1 [,bool=True]])"""
        ...
    def GetNumBonds(self, arg1: Mol, onlyHeavy: bool = True) -> int:
        """
        Returns the number of Bonds in the molecule.

        ARGUMENTS:

        onlyHeavy: (optional) include only bonds to heavy atoms (not Hs)defaults to 1.

        C++ signature :unsigned int GetNumBonds(RDKit::ROMol {lvalue} [,bool=True])"""
        ...
    def GetNumConformers(self, arg1: Mol) -> int:
        """
        Return the number of conformations on the molecule

        C++ signature :unsigned int GetNumConformers(RDKit::ROMol {lvalue})"""
        ...
    def GetNumHeavyAtoms(self, arg1: Mol) -> int:
        """
        Returns the number of heavy atoms (atomic number >1) in the molecule.

        C++ signature :unsigned int GetNumHeavyAtoms(RDKit::ROMol {lvalue})"""
        ...
    def GetProp(self, arg1: Mol, arg2: str) -> str:
        """
        Returns the value of the property.

        ARGUMENTS:
        key: the name of the property to return (a string).

        RETURNS: a string

        NOTE:
        If the property has not been set, a KeyError exception will be raised.

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetProp(RDKit::ROMol*,char const*)
        """
        ...
    def GetPropNames(
        self: Mol, includePrivate: bool = False, includeComputed: bool = False
    ) -> _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE:
        """
        Returns a tuple with all property names for this molecule.

        ARGUMENTS:

        includePrivate: (optional) toggles inclusion of private properties in the result set.Defaults to 0.

        includeComputed: (optional) toggles inclusion of computed properties in the result set.Defaults to 0.

        RETURNS: a tuple of strings

        C++ signature :std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > GetPropNames(RDKit::ROMol {lvalue} [,bool=False [,bool=False]])
        """
        ...
    def GetPropsAsDict(
        self: Mol, includePrivate: bool = False, includeComputed: bool = False
    ) -> dict:
        """
        Returns a dictionary populated with the molecules properties.n.b. Some properties are not able to be converted to python types.

        ARGUMENTS:

        includePrivate: (optional) toggles inclusion of private properties in the result set.Defaults to False.

        includeComputed: (optional) toggles inclusion of computed properties in the result set.Defaults to False.

        RETURNS: a dictionary

        C++ signature :boost::python::dict GetPropsAsDict(RDKit::ROMol [,bool=False [,bool=False]])
        """
        ...
    def GetRingInfo(self, arg1: Mol) -> RingInfo:
        """
        Returns the number of molecule’s RingInfo object.

        C++ signature :RDKit::RingInfo* GetRingInfo(RDKit::ROMol {lvalue})"""
        ...
    def GetStereoGroups(self, arg1: Mol) -> StereoGroup_vect:
        """
        Returns a list of StereoGroups defining the relative stereochemistry of the atoms.

        C++ signature :std::vector<RDKit::StereoGroup, std::allocator<RDKit::StereoGroup> > GetStereoGroups(RDKit::ROMol {lvalue})
        """
        ...
    @overload
    def GetSubstructMatch(
        self: Mol,
        query: Mol,
        useChirality: bool = False,
        useQueryQueryMatches: bool = False,
    ) -> object:
        """


        C++ signature :_object* GetSubstructMatch(RDKit::ROMol,RDKit::MolBundle [,bool=False [,bool=False]])

        GetSubstructMatch( (Mol)self, (Mol)query, (SubstructMatchParameters)params) -> object :Returns the indices of the molecule’s atoms that match a substructure query.

        ARGUMENTS:
        query: a Molecule
        params: parameters controlling the substructure match

        RETURNS: a tuple of integers

        NOTES:
        only a single match is returned

        the ordering of the indices corresponds to the atom orderingin the query. For example, the first index is for the atom in
        this molecule that matches the first atom in the query.

        C++ signature :_object* GetSubstructMatch(RDKit::ROMol,RDKit::ROMol,RDKit::SubstructMatchParameters)

        GetSubstructMatch( (Mol)self, (MolBundle)query, (SubstructMatchParameters)params) -> object :

        C++ signature :_object* GetSubstructMatch(RDKit::ROMol,RDKit::MolBundle,RDKit::SubstructMatchParameters)
        """
        ...
    @overload
    def GetSubstructMatch(
        self: Mol,
        query: MolBundle,
        useChirality: bool = False,
        useQueryQueryMatches: bool = False,
    ) -> object:
        """
        Returns the indices of the molecule’s atoms that match a substructure query.

        ARGUMENTS:
        query: a Molecule
        params: parameters controlling the substructure match

        RETURNS: a tuple of integers

        NOTES:
        only a single match is returned

        the ordering of the indices corresponds to the atom orderingin the query. For example, the first index is for the atom in
        this molecule that matches the first atom in the query.

        C++ signature :_object* GetSubstructMatch(RDKit::ROMol,RDKit::ROMol,RDKit::SubstructMatchParameters)

        GetSubstructMatch( (Mol)self, (MolBundle)query, (SubstructMatchParameters)params) -> object :

        C++ signature :_object* GetSubstructMatch(RDKit::ROMol,RDKit::MolBundle,RDKit::SubstructMatchParameters)
        """
        ...
    @overload
    def GetSubstructMatch(
        self: Mol, query: Mol, params: SubstructMatchParameters
    ) -> object:
        """


        C++ signature :_object* GetSubstructMatch(RDKit::ROMol,RDKit::MolBundle,RDKit::SubstructMatchParameters)
        """
        ...
    @overload
    def GetSubstructMatch(
        self: Mol, query: MolBundle, params: SubstructMatchParameters
    ) -> object: ...
    @overload
    def GetSubstructMatches(
        self: Mol,
        query: Mol,
        uniquify: bool = True,
        useChirality: bool = False,
        useQueryQueryMatches: bool = False,
        maxMatches: int = 1000,
    ) -> object:
        """


        C++ signature :_object* GetSubstructMatches(RDKit::ROMol,RDKit::MolBundle [,bool=True [,bool=False [,bool=False [,unsigned int=1000]]]])

        GetSubstructMatches( (Mol)self, (Mol)query, (SubstructMatchParameters)params) -> object :Returns tuples of the indices of the molecule’s atoms that match a substructure query.

        ARGUMENTS:
        query: a Molecule.
        params: parameters controlling the substructure match

        RETURNS: a tuple of tuples of integers

        NOTE:

        the ordering of the indices corresponds to the atom orderingin the query. For example, the first index is for the atom in
        this molecule that matches the first atom in the query.

        C++ signature :_object* GetSubstructMatches(RDKit::ROMol,RDKit::ROMol,RDKit::SubstructMatchParameters)

        GetSubstructMatches( (Mol)self, (MolBundle)query, (SubstructMatchParameters)params) -> object :

        C++ signature :_object* GetSubstructMatches(RDKit::ROMol,RDKit::MolBundle,RDKit::SubstructMatchParameters)
        """
        ...
    @overload
    def GetSubstructMatches(
        self: Mol,
        query: MolBundle,
        uniquify: bool = True,
        useChirality: bool = False,
        useQueryQueryMatches: bool = False,
        maxMatches: int = 1000,
    ) -> object:
        """
        Returns tuples of the indices of the molecule’s atoms that match a substructure query.

        ARGUMENTS:
        query: a Molecule.
        params: parameters controlling the substructure match

        RETURNS: a tuple of tuples of integers

        NOTE:

        the ordering of the indices corresponds to the atom orderingin the query. For example, the first index is for the atom in
        this molecule that matches the first atom in the query.

        C++ signature :_object* GetSubstructMatches(RDKit::ROMol,RDKit::ROMol,RDKit::SubstructMatchParameters)

        GetSubstructMatches( (Mol)self, (MolBundle)query, (SubstructMatchParameters)params) -> object :

        C++ signature :_object* GetSubstructMatches(RDKit::ROMol,RDKit::MolBundle,RDKit::SubstructMatchParameters)
        """
        ...
    @overload
    def GetSubstructMatches(
        self: Mol, query: Mol, params: SubstructMatchParameters
    ) -> object:
        """


        C++ signature :_object* GetSubstructMatches(RDKit::ROMol,RDKit::MolBundle,RDKit::SubstructMatchParameters)
        """
        ...
    @overload
    def GetSubstructMatches(
        self: Mol, query: MolBundle, params: SubstructMatchParameters
    ) -> object: ...
    def GetUnsignedProp(self, arg1: Mol, arg2: str) -> int:
        """
        Returns the unsigned int value of the property if possible.

        ARGUMENTS:
        key: the name of the property to return (a string).

        RETURNS: an unsigned integer

        NOTE:
        If the property has not been set, a KeyError exception will be raised.

        C++ signature :unsigned int GetUnsignedProp(RDKit::ROMol*,char const*)"""
        ...
    def HasProp(self, arg1: Mol, arg2: str) -> int:
        """
        Queries a molecule to see if a particular property has been assigned.

        ARGUMENTS:
        key: the name of the property to check for (a string).

        C++ signature :int HasProp(RDKit::ROMol,char const*)"""
        ...
    @overload
    def HasSubstructMatch(
        self: Mol,
        query: Mol,
        recursionPossible: bool = True,
        useChirality: bool = False,
        useQueryQueryMatches: bool = False,
    ) -> bool:
        """


        C++ signature :bool HasSubstructMatch(RDKit::ROMol,RDKit::MolBundle [,bool=True [,bool=False [,bool=False]]])

        HasSubstructMatch( (Mol)self, (Mol)query, (SubstructMatchParameters)params) -> bool :Queries whether or not the molecule contains a particular substructure.

        ARGUMENTS:
        query: a Molecule
        params: parameters controlling the substructure match

        RETURNS: True or False

        C++ signature :bool HasSubstructMatch(RDKit::ROMol,RDKit::ROMol,RDKit::SubstructMatchParameters)

        HasSubstructMatch( (Mol)self, (MolBundle)query [, (SubstructMatchParameters)params=True]) -> bool :

        C++ signature :bool HasSubstructMatch(RDKit::ROMol,RDKit::MolBundle [,RDKit::SubstructMatchParameters=True])
        """
        ...
    @overload
    def HasSubstructMatch(
        self: Mol,
        query: MolBundle,
        recursionPossible: bool = True,
        useChirality: bool = False,
        useQueryQueryMatches: bool = False,
    ) -> bool:
        """
        Queries whether or not the molecule contains a particular substructure.

        ARGUMENTS:
        query: a Molecule
        params: parameters controlling the substructure match

        RETURNS: True or False

        C++ signature :bool HasSubstructMatch(RDKit::ROMol,RDKit::ROMol,RDKit::SubstructMatchParameters)

        HasSubstructMatch( (Mol)self, (MolBundle)query [, (SubstructMatchParameters)params=True]) -> bool :

        C++ signature :bool HasSubstructMatch(RDKit::ROMol,RDKit::MolBundle [,RDKit::SubstructMatchParameters=True])
        """
        ...
    @overload
    def HasSubstructMatch(
        self: Mol, query: Mol, params: SubstructMatchParameters
    ) -> bool:
        """


        C++ signature :bool HasSubstructMatch(RDKit::ROMol,RDKit::MolBundle [,RDKit::SubstructMatchParameters=True])
        """
        ...
    @overload
    def HasSubstructMatch(
        self: Mol, query: MolBundle, params: SubstructMatchParameters = True
    ) -> bool: ...
    def NeedsUpdatePropertyCache(self: Mol) -> bool:
        """
        Returns true or false depending on whether implicit and explicit valence of the molecule have already been calculated.

        C++ signature :bool NeedsUpdatePropertyCache(RDKit::ROMol {lvalue})"""
        ...
    def RemoveAllConformers(self, arg1: Mol) -> None:
        """
        Remove all the conformations on the molecule

        C++ signature :void RemoveAllConformers(RDKit::ROMol {lvalue})"""
        ...
    def RemoveConformer(self, arg1: Mol, arg2: int) -> None:
        """
        Remove the conformer with the specified ID

        C++ signature :void RemoveConformer(RDKit::ROMol {lvalue},unsigned int)"""
        ...
    def SetBoolProp(self: Mol, key: str, val: bool, computed: bool = False) -> None:
        """
        Sets a boolean valued molecular property

        ARGUMENTS:
        key: the name of the property to be set (a string).
        value: the property value as a bool.

        computed: (optional) marks the property as being computed.Defaults to False.

        C++ signature :void SetBoolProp(RDKit::ROMol,char const*,bool [,bool=False])"""
        ...
    def SetDoubleProp(self: Mol, key: str, val: float, computed: bool = False) -> None:
        """
        Sets a double valued molecular property

        ARGUMENTS:
        key: the name of the property to be set (a string).
        value: the property value as a double.

        computed: (optional) marks the property as being computed.Defaults to 0.

        C++ signature :void SetDoubleProp(RDKit::ROMol,char const*,double [,bool=False])
        """
        ...
    def SetIntProp(self: Mol, key: str, val: int, computed: bool = False) -> None:
        """
        Sets an integer valued molecular property

        ARGUMENTS:
        key: the name of the property to be set (an unsigned number).
        value: the property value as an integer.

        computed: (optional) marks the property as being computed.Defaults to False.

        C++ signature :void SetIntProp(RDKit::ROMol,char const*,int [,bool=False])"""
        ...
    def SetProp(self: Mol, key: str, val: str, computed: bool = False) -> None:
        """
        Sets a molecular property

        ARGUMENTS:
        key: the name of the property to be set (a string).
        value: the property value (a string).

        computed: (optional) marks the property as being computed.Defaults to False.

        C++ signature :void SetProp(RDKit::ROMol,char const*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])
        """
        ...
    def SetUnsignedProp(self: Mol, key: str, val: int, computed: bool = False) -> None:
        """
        Sets an unsigned integer valued molecular property

        ARGUMENTS:
        key: the name of the property to be set (a string).
        value: the property value as an unsigned integer.

        computed: (optional) marks the property as being computed.Defaults to False.

        C++ signature :void SetUnsignedProp(RDKit::ROMol,char const*,unsigned int [,bool=False])
        """
        ...
    @overload
    def ToBinary(self, arg1: Mol) -> object:
        """
        Returns a binary string representation of the molecule pickling the specified properties.

        C++ signature :boost::python::api::object ToBinary(RDKit::ROMol,unsigned int)"""
        ...
    @overload
    def ToBinary(self, mol: Mol, propertyFlags: int) -> object: ...
    def UpdatePropertyCache(self: Mol, strict: bool = True) -> None:
        """
        Regenerates computed properties like implicit valence and ring information.

        C++ signature :void UpdatePropertyCache(RDKit::ROMol {lvalue} [,bool=True])"""
        ...
    @classmethod
    def __copy__(cls, boost) -> Any: ...
    @classmethod
    def __deepcopy__(cls) -> Any: ...
    @classmethod
    def __getinitargs__(cls, RDKit) -> Any: ...
    @classmethod
    def __getstate__(cls, boost) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @classmethod
    def __setstate__(cls, state) -> Any: ...

class MolBundle(Boost.Python.instance):
    """
    A class for storing groups of related molecules.

    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def AddMol(self, arg1: MolBundle, arg2: Mol) -> int:
        """
        C++ signature :unsigned long AddMol(RDKit::MolBundle {lvalue},boost::shared_ptr<RDKit::ROMol>)
        """
        ...
    def GetMol(self, arg1: MolBundle, arg2: int) -> Mol:
        """
        C++ signature :boost::shared_ptr<RDKit::ROMol> GetMol(RDKit::MolBundle {lvalue},unsigned long)
        """
        ...
    @overload
    def GetSubstructMatch(
        self: MolBundle,
        query: Mol,
        useChirality: bool = False,
        useQueryQueryMatches: bool = False,
    ) -> object:
        """
        Returns the indices of the atoms from the first molecule in a bundle that matches a substructure query from a bundle.

        ARGUMENTS:
        query: a MolBundle
        useChirality: enables the use of stereochemistry in the matching
        useQueryQueryMatches: use query-query matching logic

        RETURNS: a tuple of integers

        NOTES:
        only a single match is returned

        the ordering of the indices corresponds to the atom orderingin the query. For example, the first index is for the atom in
        this molecule that matches the first atom in the query.

        C++ signature :_object* GetSubstructMatch(RDKit::MolBundle,RDKit::MolBundle [,bool=False [,bool=False]])

        GetSubstructMatch( (MolBundle)self, (Mol)query, (SubstructMatchParameters)params) -> object :Returns the indices of the atoms from the first molecule in a bundle that matches a substructure query.

        ARGUMENTS:
        query: a Molecule
        params: parameters controlling the substructure match

        RETURNS: a tuple of integers

        NOTES:
        only a single match is returned

        the ordering of the indices corresponds to the atom orderingin the query. For example, the first index is for the atom in
        this molecule that matches the first atom in the query.

        C++ signature :_object* GetSubstructMatch(RDKit::MolBundle,RDKit::ROMol,RDKit::SubstructMatchParameters)

        GetSubstructMatch( (MolBundle)self, (MolBundle)query, (SubstructMatchParameters)params) -> object :Returns the indices of the atoms from the first molecule in a bundle that matches a substructure query from a bundle.

        ARGUMENTS:
        query: a MolBundle
        params: parameters controlling the substructure match

        RETURNS: a tuple of integers

        NOTES:
        only a single match is returned

        the ordering of the indices corresponds to the atom orderingin the query. For example, the first index is for the atom in
        this molecule that matches the first atom in the query.

        C++ signature :_object* GetSubstructMatch(RDKit::MolBundle,RDKit::MolBundle,RDKit::SubstructMatchParameters)
        """
        ...
    @overload
    def GetSubstructMatch(
        self: MolBundle,
        query: MolBundle,
        useChirality: bool = False,
        useQueryQueryMatches: bool = False,
    ) -> object:
        """
        Returns the indices of the atoms from the first molecule in a bundle that matches a substructure query.

        ARGUMENTS:
        query: a Molecule
        params: parameters controlling the substructure match

        RETURNS: a tuple of integers

        NOTES:
        only a single match is returned

        the ordering of the indices corresponds to the atom orderingin the query. For example, the first index is for the atom in
        this molecule that matches the first atom in the query.

        C++ signature :_object* GetSubstructMatch(RDKit::MolBundle,RDKit::ROMol,RDKit::SubstructMatchParameters)

        GetSubstructMatch( (MolBundle)self, (MolBundle)query, (SubstructMatchParameters)params) -> object :Returns the indices of the atoms from the first molecule in a bundle that matches a substructure query from a bundle.

        ARGUMENTS:
        query: a MolBundle
        params: parameters controlling the substructure match

        RETURNS: a tuple of integers

        NOTES:
        only a single match is returned

        the ordering of the indices corresponds to the atom orderingin the query. For example, the first index is for the atom in
        this molecule that matches the first atom in the query.

        C++ signature :_object* GetSubstructMatch(RDKit::MolBundle,RDKit::MolBundle,RDKit::SubstructMatchParameters)
        """
        ...
    @overload
    def GetSubstructMatch(
        self: MolBundle, query: Mol, params: SubstructMatchParameters
    ) -> object:
        """
        Returns the indices of the atoms from the first molecule in a bundle that matches a substructure query from a bundle.

        ARGUMENTS:
        query: a MolBundle
        params: parameters controlling the substructure match

        RETURNS: a tuple of integers

        NOTES:
        only a single match is returned

        the ordering of the indices corresponds to the atom orderingin the query. For example, the first index is for the atom in
        this molecule that matches the first atom in the query.

        C++ signature :_object* GetSubstructMatch(RDKit::MolBundle,RDKit::MolBundle,RDKit::SubstructMatchParameters)
        """
        ...
    @overload
    def GetSubstructMatch(
        self: MolBundle, query: MolBundle, params: SubstructMatchParameters
    ) -> object: ...
    @overload
    def GetSubstructMatches(
        self: MolBundle,
        query: Mol,
        uniquify: bool = True,
        useChirality: bool = False,
        useQueryQueryMatches: bool = False,
        maxMatches: int = 1000,
    ) -> object:
        """
        Returns tuple of all indices of the atoms from the first molecule in a bundle that matches a substructure query from the second bundle.

        ARGUMENTS:
        query: a MolBundle.

        uniquify: (optional) determines whether or not the matches are uniquified.Defaults to 1.

        useChirality: enables the use of stereochemistry in the matching
        useQueryQueryMatches: use query-query matching logic

        maxMatches: The maximum number of matches that will be returned.In high-symmetry cases with medium-sized molecules, it is
        very easy to end up with a combinatorial explosion in the
        number of possible matches. This argument prevents that from
        having unintended consequences

        RETURNS: a tuple of tuples of integers

        NOTE:

        the ordering of the indices corresponds to the atom orderingin the query. For example, the first index is for the atom in
        this molecule that matches the first atom in the query.

        C++ signature :_object* GetSubstructMatches(RDKit::MolBundle,RDKit::MolBundle [,bool=True [,bool=False [,bool=False [,unsigned int=1000]]]])

        GetSubstructMatches( (MolBundle)self, (Mol)query, (SubstructMatchParameters)params) -> object :Returns tuple of all indices of the atoms from the first molecule in a bundle that matches a substructure query.

        ARGUMENTS:
        query: a molecule.
        params: parameters controlling the substructure match

        RETURNS: a tuple of tuples of integers

        NOTE:

        the ordering of the indices corresponds to the atom orderingin the query. For example, the first index is for the atom in
        this molecule that matches the first atom in the query.

        C++ signature :_object* GetSubstructMatches(RDKit::MolBundle,RDKit::ROMol,RDKit::SubstructMatchParameters)

        GetSubstructMatches( (MolBundle)self, (MolBundle)query, (SubstructMatchParameters)params) -> object :Returns tuple of all indices of the atoms from the first molecule in a bundle that matches a substructure query from the second bundle.

        ARGUMENTS:
        query: a MolBundle.
        params: parameters controlling the substructure match

        RETURNS: a tuple of tuples of integers

        NOTE:

        the ordering of the indices corresponds to the atom orderingin the query. For example, the first index is for the atom in
        this molecule that matches the first atom in the query.

        C++ signature :_object* GetSubstructMatches(RDKit::MolBundle,RDKit::MolBundle,RDKit::SubstructMatchParameters)
        """
        ...
    @overload
    def GetSubstructMatches(
        self: MolBundle,
        query: MolBundle,
        uniquify: bool = True,
        useChirality: bool = False,
        useQueryQueryMatches: bool = False,
        maxMatches: int = 1000,
    ) -> object:
        """
        Returns tuple of all indices of the atoms from the first molecule in a bundle that matches a substructure query.

        ARGUMENTS:
        query: a molecule.
        params: parameters controlling the substructure match

        RETURNS: a tuple of tuples of integers

        NOTE:

        the ordering of the indices corresponds to the atom orderingin the query. For example, the first index is for the atom in
        this molecule that matches the first atom in the query.

        C++ signature :_object* GetSubstructMatches(RDKit::MolBundle,RDKit::ROMol,RDKit::SubstructMatchParameters)

        GetSubstructMatches( (MolBundle)self, (MolBundle)query, (SubstructMatchParameters)params) -> object :Returns tuple of all indices of the atoms from the first molecule in a bundle that matches a substructure query from the second bundle.

        ARGUMENTS:
        query: a MolBundle.
        params: parameters controlling the substructure match

        RETURNS: a tuple of tuples of integers

        NOTE:

        the ordering of the indices corresponds to the atom orderingin the query. For example, the first index is for the atom in
        this molecule that matches the first atom in the query.

        C++ signature :_object* GetSubstructMatches(RDKit::MolBundle,RDKit::MolBundle,RDKit::SubstructMatchParameters)
        """
        ...
    @overload
    def GetSubstructMatches(
        self: MolBundle, query: Mol, params: SubstructMatchParameters
    ) -> object:
        """
        Returns tuple of all indices of the atoms from the first molecule in a bundle that matches a substructure query from the second bundle.

        ARGUMENTS:
        query: a MolBundle.
        params: parameters controlling the substructure match

        RETURNS: a tuple of tuples of integers

        NOTE:

        the ordering of the indices corresponds to the atom orderingin the query. For example, the first index is for the atom in
        this molecule that matches the first atom in the query.

        C++ signature :_object* GetSubstructMatches(RDKit::MolBundle,RDKit::MolBundle,RDKit::SubstructMatchParameters)
        """
        ...
    @overload
    def GetSubstructMatches(
        self: MolBundle, query: MolBundle, params: SubstructMatchParameters
    ) -> object: ...
    @overload
    def HasSubstructMatch(
        self: MolBundle,
        query: Mol,
        recursionPossible: bool = True,
        useChirality: bool = False,
        useQueryQueryMatches: bool = False,
    ) -> bool:
        """
        Queries whether or not any molecule in the first bundle matches any molecule in the second bundle.

        ARGUMENTS:
        query: a MolBundle
        recursionPossible: (optional)
        useChirality: enables the use of stereochemistry in the matching
        useQueryQueryMatches: use query-query matching logic

        RETURNS: True or False

        C++ signature :bool HasSubstructMatch(RDKit::MolBundle,RDKit::MolBundle [,bool=True [,bool=False [,bool=False]]])

        HasSubstructMatch( (MolBundle)self, (Mol)query, (SubstructMatchParameters)params) -> bool :Queries whether or not any molecule in the bundle contains a particular substructure.

        ARGUMENTS:
        query: a Molecule
        params: parameters controlling the substructure match

        matching

        useQueryQueryMatches: use query-query matching logic

        RETURNS: True or False

        C++ signature :bool HasSubstructMatch(RDKit::MolBundle,RDKit::ROMol,RDKit::SubstructMatchParameters)

        HasSubstructMatch( (MolBundle)self, (MolBundle)query, (SubstructMatchParameters)params) -> bool :Queries whether or not any molecule in the first bundle matches any molecule in the second bundle.

        ARGUMENTS:
        query: a MolBundle
        params: parameters controlling the substructure match

        RETURNS: True or False

        C++ signature :bool HasSubstructMatch(RDKit::MolBundle,RDKit::MolBundle,RDKit::SubstructMatchParameters)
        """
        ...
    @overload
    def HasSubstructMatch(
        self: MolBundle,
        query: MolBundle,
        recursionPossible: bool = True,
        useChirality: bool = False,
        useQueryQueryMatches: bool = False,
    ) -> bool:
        """
        Queries whether or not any molecule in the bundle contains a particular substructure.

        ARGUMENTS:
        query: a Molecule
        params: parameters controlling the substructure match

        matching

        useQueryQueryMatches: use query-query matching logic

        RETURNS: True or False

        C++ signature :bool HasSubstructMatch(RDKit::MolBundle,RDKit::ROMol,RDKit::SubstructMatchParameters)

        HasSubstructMatch( (MolBundle)self, (MolBundle)query, (SubstructMatchParameters)params) -> bool :Queries whether or not any molecule in the first bundle matches any molecule in the second bundle.

        ARGUMENTS:
        query: a MolBundle
        params: parameters controlling the substructure match

        RETURNS: True or False

        C++ signature :bool HasSubstructMatch(RDKit::MolBundle,RDKit::MolBundle,RDKit::SubstructMatchParameters)
        """
        ...
    @overload
    def HasSubstructMatch(
        self: MolBundle, query: Mol, params: SubstructMatchParameters
    ) -> bool:
        """
        Queries whether or not any molecule in the first bundle matches any molecule in the second bundle.

        ARGUMENTS:
        query: a MolBundle
        params: parameters controlling the substructure match

        RETURNS: True or False

        C++ signature :bool HasSubstructMatch(RDKit::MolBundle,RDKit::MolBundle,RDKit::SubstructMatchParameters)
        """
        ...
    @overload
    def HasSubstructMatch(
        self: MolBundle, query: MolBundle, params: SubstructMatchParameters
    ) -> bool: ...
    def Size(self, arg1: MolBundle) -> int:
        """
        C++ signature :unsigned long Size(RDKit::MolBundle {lvalue})"""
        ...
    @classmethod
    def __getitem__(cls, RDKit, unsignedlong) -> Any: ...
    @classmethod
    def __len__(cls, RDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class MolSanitizeException(ValueError):
    ...
    ...

class PeriodicTable(Boost.Python.instance):
    """
    A class which stores information from the Periodic Table.
    It is not possible to create a PeriodicTable object directly from Python,
    use GetPeriodicTable() to get the global table.
    The PeriodicTable object can be queried for a variety of properties:

    GetAtomicWeight
    GetAtomicNumber
    GetElementSymbol
    GetElementName
    GetRvdw (van der Waals radius)
    GetRCovalent (covalent radius)
    GetDefaultValence
    GetValenceList
    GetNOuterElecs (number of valence electrons)
    GetMostCommonIsotope
    GetMostCommonIsotopeMass
    GetRb0
    GetAbundanceForIsotope
    GetMassForIsotope

    When it makes sense, these can be queried using either an atomic number (integer)
    or an atomic symbol (string)
    Raises an exception
    This class cannot be instantiated from Python"""

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @overload
    def GetAbundanceForIsotope(
        self, arg1: PeriodicTable, arg2: int, arg3: int
    ) -> float:
        """


        C++ signature :double GetAbundanceForIsotope(RDKit::PeriodicTable {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,unsigned int)
        """
        ...
    @overload
    def GetAbundanceForIsotope(
        self, arg1: PeriodicTable, arg2: str, arg3: int
    ) -> float: ...
    def GetAtomicNumber(self, arg1: PeriodicTable, arg2: str) -> int:
        """
        C++ signature :int GetAtomicNumber(RDKit::PeriodicTable {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    @overload
    def GetAtomicWeight(self, arg1: PeriodicTable, arg2: int) -> float:
        """


        C++ signature :double GetAtomicWeight(RDKit::PeriodicTable {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    @overload
    def GetAtomicWeight(self, arg1: PeriodicTable, arg2: str) -> float: ...
    @overload
    def GetDefaultValence(self, arg1: PeriodicTable, arg2: int) -> int:
        """


        C++ signature :int GetDefaultValence(RDKit::PeriodicTable {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    @overload
    def GetDefaultValence(self, arg1: PeriodicTable, arg2: str) -> int: ...
    def GetElementName(self, arg1: PeriodicTable, arg2: int) -> str:
        """
        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetElementName(RDKit::PeriodicTable {lvalue},unsigned int)
        """
        ...
    def GetElementSymbol(self, arg1: PeriodicTable, arg2: int) -> str:
        """
        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetElementSymbol(RDKit::PeriodicTable {lvalue},unsigned int)
        """
        ...
    @overload
    def GetMassForIsotope(self, arg1: PeriodicTable, arg2: int, arg3: int) -> float:
        """


        C++ signature :double GetMassForIsotope(RDKit::PeriodicTable {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,unsigned int)
        """
        ...
    @overload
    def GetMassForIsotope(self, arg1: PeriodicTable, arg2: str, arg3: int) -> float: ...
    @overload
    def GetMostCommonIsotope(self, arg1: PeriodicTable, arg2: int) -> int:
        """


        C++ signature :int GetMostCommonIsotope(RDKit::PeriodicTable {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    @overload
    def GetMostCommonIsotope(self, arg1: PeriodicTable, arg2: str) -> int: ...
    @overload
    def GetMostCommonIsotopeMass(self, arg1: PeriodicTable, arg2: int) -> float:
        """


        C++ signature :double GetMostCommonIsotopeMass(RDKit::PeriodicTable {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    @overload
    def GetMostCommonIsotopeMass(self, arg1: PeriodicTable, arg2: str) -> float: ...
    @overload
    def GetNOuterElecs(self, arg1: PeriodicTable, arg2: int) -> int:
        """


        C++ signature :int GetNOuterElecs(RDKit::PeriodicTable {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    @overload
    def GetNOuterElecs(self, arg1: PeriodicTable, arg2: str) -> int: ...
    @overload
    def GetRb0(self, arg1: PeriodicTable, arg2: int) -> float:
        """


        C++ signature :double GetRb0(RDKit::PeriodicTable {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    @overload
    def GetRb0(self, arg1: PeriodicTable, arg2: str) -> float: ...
    @overload
    def GetRcovalent(self, arg1: PeriodicTable, arg2: int) -> float:
        """


        C++ signature :double GetRcovalent(RDKit::PeriodicTable {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    @overload
    def GetRcovalent(self, arg1: PeriodicTable, arg2: str) -> float: ...
    @overload
    def GetRvdw(self, arg1: PeriodicTable, arg2: int) -> float:
        """


        C++ signature :double GetRvdw(RDKit::PeriodicTable {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    @overload
    def GetRvdw(self, arg1: PeriodicTable, arg2: str) -> float: ...
    @overload
    def GetValenceList(self, arg1: PeriodicTable, arg2: int) -> _vecti:
        """


        C++ signature :std::vector<int, std::allocator<int> > GetValenceList(RDKit::PeriodicTable {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    @overload
    def GetValenceList(self, arg1: PeriodicTable, arg2: str) -> _vecti: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class PropertyPickleOptions(Boost.Python.enum):
    AllProps: PropertyPickleOptions = ...
    AtomProps: PropertyPickleOptions = ...
    BondProps: PropertyPickleOptions = ...
    ComputedProps: PropertyPickleOptions = ...
    CoordsAsDouble: PropertyPickleOptions = ...
    MolProps: PropertyPickleOptions = ...
    NoConformers: PropertyPickleOptions = ...
    NoProps: PropertyPickleOptions = ...
    PrivateProps: PropertyPickleOptions = ...
    QueryAtomData: PropertyPickleOptions = ...
    names: dict[str, PropertyPickleOptions] = ...
    values: dict[int, PropertyPickleOptions] = ...
    __slots__: ClassVar[tuple] = ...

class QueryAtom(Atom):
    """
    The class to store QueryAtoms.
    These cannot currently be constructed directly from Python
    Raises an exception
    This class cannot be instantiated from Python"""

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def ExpandQuery(
        self: QueryAtom,
        other: QueryAtom,
        how: CompositeQueryType = CompositeQueryType.COMPOSITE_AND,
        maintainOrder: bool = True,
    ) -> None:
        """
        combines the query from other with ours

        C++ signature :void ExpandQuery(RDKit::QueryAtom*,RDKit::QueryAtom const* [,Queries::CompositeQueryType=rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_AND [,bool=True]])
        """
        ...
    def SetQuery(self: QueryAtom, other: QueryAtom) -> None:
        """
        Replace our query with a copy of the other query

        C++ signature :void SetQuery(RDKit::QueryAtom*,RDKit::QueryAtom const*)"""
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class QueryBond(Bond):
    """
    The class to store QueryBonds.
    These cannot currently be constructed directly from Python
    Raises an exception
    This class cannot be instantiated from Python"""

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def ExpandQuery(
        self: QueryBond,
        other: QueryBond,
        how: CompositeQueryType = CompositeQueryType.COMPOSITE_AND,
        maintainOrder: bool = True,
    ) -> None:
        """
        combines the query from other with ours

        C++ signature :void ExpandQuery(RDKit::QueryBond*,RDKit::QueryBond const* [,Queries::CompositeQueryType=rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_AND [,bool=True]])
        """
        ...
    def SetQuery(self: QueryBond, other: QueryBond) -> None:
        """
        Replace our query with a copy of the other query

        C++ signature :void SetQuery(RDKit::QueryBond*,RDKit::QueryBond const*)"""
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class RWMol(Mol):
    """
    The RW molecule class (read/write)
    This class is a more-performant version of the EditableMolecule class in that
    it is a ‘live’ molecule and shares the interface from the Mol class.
    All changes are performed without the need to create a copy of the
    molecule using GetMol() (this is still available, however).
    n.b. Eventually this class may become a direct replacement for EditableMol
    Construct from a Mol

    C++ signature :void __init__(_object*,RDKit::ROMol)

    __init__( (object)arg1) -> None :

    C++ signature :void __init__(_object*)

    __init__( (object)arg1, (str)pklString) -> None :

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)

    __init__( (object)arg1, (str)pklString, (int)propertyFlags) -> None :

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,unsigned int)

    __init__( (object)arg1, (Mol)mol [, (bool)quickCopy=False [, (int)confId=-1]]) -> None :

    C++ signature :void __init__(_object*,RDKit::ROMol [,bool=False [,int=-1]])"""

    __getstate_manages_dict__: ClassVar[bool] = ...
    __instance_size__: ClassVar[int] = ...
    __safe_for_unpickling__: ClassVar[bool] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def AddAtom(self, arg1: RWMol, atom: Atom) -> int:
        """
        add an atom, returns the index of the newly added atom

        C++ signature :int AddAtom(RDKit::ReadWriteMol {lvalue},RDKit::Atom*)"""
        ...
    def AddBond(
        self,
        arg1: RWMol,
        beginAtomIdx: int,
        endAtomIdx: int,
        order: BondType = BondType.UNSPECIFIED,
    ) -> int:
        """
        add a bond, returns the new number of bonds

        C++ signature :int AddBond(RDKit::ReadWriteMol {lvalue},unsigned int,unsigned int [,RDKit::Bond::BondType=rdkit.Chem.rdchem.BondType.UNSPECIFIED])
        """
        ...
    def BeginBatchEdit(self, arg1: RWMol) -> None:
        """
        starts batch editing

        C++ signature :void BeginBatchEdit(RDKit::ReadWriteMol {lvalue})"""
        ...
    def CommitBatchEdit(self, arg1: RWMol) -> None:
        """
        finishes batch editing and makes the actual changes

        C++ signature :void CommitBatchEdit(RDKit::ReadWriteMol {lvalue})"""
        ...
    def GetMol(self, arg1: RWMol) -> Mol:
        """
        Returns a Mol (a normal molecule)

        C++ signature :RDKit::ROMol* GetMol(RDKit::ReadWriteMol {lvalue})"""
        ...
    def InsertMol(self, arg1: RWMol, mol: Mol) -> None:
        """
        Insert (add) the given molecule into this one

        C++ signature :void InsertMol(RDKit::ReadWriteMol {lvalue},RDKit::ROMol)"""
        ...
    def RemoveAtom(self, arg1: RWMol, arg2: int) -> None:
        """
        Remove the specified atom from the molecule

        C++ signature :void RemoveAtom(RDKit::ReadWriteMol {lvalue},unsigned int)"""
        ...
    def RemoveBond(self, arg1: RWMol, arg2: int, arg3: int) -> None:
        """
        Remove the specified bond from the molecule

        C++ signature :void RemoveBond(RDKit::ReadWriteMol {lvalue},unsigned int,unsigned int)
        """
        ...
    def ReplaceAtom(
        self,
        arg1: RWMol,
        index: int,
        newAtom: Atom,
        updateLabel: bool = False,
        preserveProps: bool = False,
    ) -> None:
        """
        replaces the specified atom with the provided one
        If updateLabel is True, the new atom becomes the active atom
        If preserveProps is True preserve keep the existing props unless explicit set on the new atom

        C++ signature :void ReplaceAtom(RDKit::ReadWriteMol {lvalue},unsigned int,RDKit::Atom* [,bool=False [,bool=False]])
        """
        ...
    def ReplaceBond(
        self,
        arg1: RWMol,
        index: int,
        newBond: Bond,
        preserveProps: bool = False,
        keepSGroups: bool = True,
    ) -> None:
        """
        replaces the specified bond with the provided one.
        If preserveProps is True preserve keep the existing props unless explicit set on the new bond. If keepSGroups is False, allSubstance Groups referencing the bond will be dropped.

        C++ signature :void ReplaceBond(RDKit::ReadWriteMol {lvalue},unsigned int,RDKit::Bond* [,bool=False [,bool=True]])
        """
        ...
    def RollbackBatchEdit(self, arg1: RWMol) -> None:
        """
        cancels batch editing

        C++ signature :void RollbackBatchEdit(RDKit::ReadWriteMol {lvalue})"""
        ...
    def SetStereoGroups(self, arg1: RWMol, stereo_groups: list) -> None:
        """
        Set the stereo groups

        C++ signature :void SetStereoGroups(RDKit::ReadWriteMol {lvalue},boost::python::list {lvalue})
        """
        ...
    @classmethod
    def __copy__(cls, boost) -> Any: ...
    @classmethod
    def __deepcopy__(cls) -> Any: ...
    @classmethod
    def __enter__(cls, RDKit) -> Any: ...
    @classmethod
    def __exit__(cls, type, value, traceback) -> Any: ...
    @classmethod
    def __getinitargs__(cls, RDKit) -> Any: ...
    @classmethod
    def __getstate__(cls, boost) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @classmethod
    def __setstate__(cls, state) -> Any: ...

class ResonanceFlags(Boost.Python.enum):
    ALLOW_CHARGE_SEPARATION: ResonanceFlags = ...
    ALLOW_INCOMPLETE_OCTETS: ResonanceFlags = ...
    KEKULE_ALL: ResonanceFlags = ...
    UNCONSTRAINED_ANIONS: ResonanceFlags = ...
    UNCONSTRAINED_CATIONS: ResonanceFlags = ...
    names: dict[str, ResonanceFlags] = ...
    values: dict[int, ResonanceFlags] = ...
    __slots__: ClassVar[tuple] = ...

class ResonanceMolSupplier(Boost.Python.instance):
    """
    A class which supplies resonance structures (as mols) from a mol.
    Usage examples:

    Lazy evaluation: the resonance structures are not constructed
    until we ask for them:
    >>> suppl = ResonanceMolSupplier(mol)
    >>> for resMol in suppl:
    ...    resMol.GetNumAtoms()

    Lazy evaluation 2:
    >>> suppl = ResonanceMolSupplier(mol)
    >>> resMol1 = next(suppl)
    >>> resMol2 = next(suppl)
    >>> suppl.reset()
    >>> resMol3 = next(suppl)
    # resMol3 and resMol1 are the same:
    >>> MolToSmiles(resMol3)==MolToSmiles(resMol1)

    Random Access:
    >>> suppl = ResonanceMolSupplier(mol)
    >>> resMol1 = suppl[0]
    >>> resMol2 = suppl[1]

    NOTE: this will generate an IndexError if the supplier doesn’t have that many
    molecules.

    Random Access 2: looping over all resonance structures
    >>> suppl = ResonanceMolSupplier(mol)
    >>> nResMols = len(suppl)
    >>> for i in range(nResMols):
    …   suppl[i].GetNumAtoms()

    C++ signature :void __init__(_object*,RDKit::ROMol {lvalue} [,unsigned int=0 [,unsigned int=1000]])
    """

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def Enumerate(self, arg1: ResonanceMolSupplier) -> None:
        """
        Ask ResonanceMolSupplier to enumerate resonance structures(automatically done as soon as any attempt to access them is made).

        C++ signature :void Enumerate(RDKit::ResonanceMolSupplier {lvalue})"""
        ...
    def GetAtomConjGrpIdx(self, arg1: ResonanceMolSupplier, arg2: int) -> int:
        """
        Given an atom index, it returns the index of the conjugated groupthe atom belongs to, or -1 if it is not conjugated.

        C++ signature :unsigned int GetAtomConjGrpIdx(RDKit::ResonanceMolSupplier {lvalue},unsigned int)
        """
        ...
    def GetBondConjGrpIdx(self, arg1: ResonanceMolSupplier, arg2: int) -> int:
        """
        Given a bond index, it returns the index of the conjugated groupthe bond belongs to, or -1 if it is not conjugated.

        C++ signature :unsigned int GetBondConjGrpIdx(RDKit::ResonanceMolSupplier {lvalue},unsigned int)
        """
        ...
    def GetIsEnumerated(self, arg1: ResonanceMolSupplier) -> bool:
        """
        Returns true if resonance structure enumeration has already happened.

        C++ signature :bool GetIsEnumerated(RDKit::ResonanceMolSupplier {lvalue})"""
        ...
    def GetNumConjGrps(self, arg1: ResonanceMolSupplier) -> int:
        """
        Returns the number of individual conjugated groups in the molecule.

        C++ signature :unsigned int GetNumConjGrps(RDKit::ResonanceMolSupplier {lvalue})
        """
        ...
    def GetProgressCallback(self, arg1: ResonanceMolSupplier) -> object:
        """
        Get the ResonanceMolSupplierCallback subclass instance,
        or None if none was set.

        C++ signature :boost::python::api::object GetProgressCallback(RDKit::ResonanceMolSupplier)
        """
        ...
    def GetSubstructMatch(
        self: ResonanceMolSupplier,
        query: Mol,
        useChirality: bool = False,
        useQueryQueryMatches: bool = False,
    ) -> object:
        """
        Returns the indices of the molecule’s atoms that match a substructure query,
        taking into account all resonance structures in ResonanceMolSupplier.

        ARGUMENTS:
        query: a Molecule
        useChirality: enables the use of stereochemistry in the matching
        useQueryQueryMatches: use query-query matching logic

        RETURNS: a tuple of integers

        NOTES:
        only a single match is returned

        the ordering of the indices corresponds to the atom orderingin the query. For example, the first index is for the atom in
        this molecule that matches the first atom in the query.

        C++ signature :_object* GetSubstructMatch(RDKit::ResonanceMolSupplier {lvalue},RDKit::ROMol [,bool=False [,bool=False]])
        """
        ...
    def GetSubstructMatches(
        self: ResonanceMolSupplier,
        query: Mol,
        uniquify: bool = False,
        useChirality: bool = False,
        useQueryQueryMatches: bool = False,
        maxMatches: int = 1000,
        numThreads: int = 1,
    ) -> object:
        """
        Returns tuples of the indices of the molecule’s atoms that match a substructure query,
        taking into account all resonance structures in ResonanceMolSupplier.

        ARGUMENTS:
        query: a Molecule.

        uniquify: (optional) determines whether or not the matches are uniquified.Defaults to 1.

        useChirality: enables the use of stereochemistry in the matching
        useQueryQueryMatches: use query-query matching logic

        maxMatches: The maximum number of matches that will be returned.In high-symmetry cases with medium-sized molecules, it is
        very easy to end up with a combinatorial explosion in the
        number of possible matches. This argument prevents that from
        having unintended consequences

        numThreads: The number of threads to be used (defaults to 1; 0 selects thenumber of concurrent threads supported by the hardware; negative
        values are added to the number of concurrent threads supported
        by the hardware).

        RETURNS: a tuple of tuples of integers

        NOTE:

        the ordering of the indices corresponds to the atom orderingin the query. For example, the first index is for the atom in
        this molecule that matches the first atom in the query.

        C++ signature :_object* GetSubstructMatches(RDKit::ResonanceMolSupplier {lvalue},RDKit::ROMol [,bool=False [,bool=False [,bool=False [,unsigned int=1000 [,int=1]]]]])
        """
        ...
    def SetNumThreads(self, arg1: ResonanceMolSupplier, arg2: int) -> None:
        """
        Sets the number of threads to be used to enumerate resonance
        structures (defaults to 1; 0 selects the number of concurrent
        threads supported by the hardware; negative values are added
        to the number of concurrent threads supported by the hardware).

        C++ signature :void SetNumThreads(RDKit::ResonanceMolSupplier {lvalue},unsigned int)
        """
        ...
    def SetProgressCallback(self, arg1: ResonanceMolSupplier, arg2: object) -> None:
        """
        Pass an instance of a class derived from
        ResonanceMolSupplierCallback, which must implement the
        __call__() method.

        C++ signature :void SetProgressCallback(RDKit::ResonanceMolSupplier {lvalue},_object*)
        """
        ...
    def WasCanceled(self, arg1: ResonanceMolSupplier) -> bool:
        """
        Returns True if the resonance structure generation was canceled.

        C++ signature :bool WasCanceled(RDKit::ResonanceMolSupplier {lvalue})"""
        ...
    def atEnd(self, arg1: ResonanceMolSupplier) -> bool:
        """
        Returns whether or not we have hit the end of the resonance structure supplier.

        C++ signature :bool atEnd(RDKit::ResonanceMolSupplier {lvalue})"""
        ...
    def reset(self, arg1: ResonanceMolSupplier) -> None:
        """
        Resets our position in the resonance structure supplier to the beginning.

        C++ signature :void reset(RDKit::ResonanceMolSupplier {lvalue})"""
        ...
    @classmethod
    def __getitem__(cls, RDKit, int) -> Any: ...
    @classmethod
    def __iter__(cls, RDKit) -> Any: ...
    @classmethod
    def __len__(cls, RDKit) -> Any: ...
    @classmethod
    def __next__(cls, RDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class ResonanceMolSupplierCallback(Boost.Python.instance):
    """
    Create a derived class from this abstract base class and
    implement the __call__() method.
    The __call__() method is called at each iteration of the
    algorithm, and provides a mechanism to monitor or stop
    its progress.
    To have your callback called, pass an instance of your
    derived class to ResonanceMolSupplier.SetProgressCallback()

    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def GetMaxStructures(self, arg1: ResonanceMolSupplierCallback) -> int:
        """
        Get the number of conjugated groups this molecule has.

        C++ signature :unsigned long GetMaxStructures(RDKit::PyResonanceMolSupplierCallback {lvalue})
        """
        ...
    def GetNumConjGrps(self, arg1: ResonanceMolSupplierCallback) -> int:
        """
        Returns the number of individual conjugated groups in the molecule.

        C++ signature :unsigned int GetNumConjGrps(RDKit::PyResonanceMolSupplierCallback {lvalue})
        """
        ...
    def GetNumDiverseStructures(
        self, arg1: ResonanceMolSupplierCallback, arg2: int
    ) -> int:
        """
        Get the number of non-degenrate resonance structures generated so far for the passed conjugated group index.

        C++ signature :unsigned long GetNumDiverseStructures(RDKit::PyResonanceMolSupplierCallback {lvalue},unsigned int)
        """
        ...
    def GetNumStructures(self, arg1: ResonanceMolSupplierCallback, arg2: int) -> int:
        """
        Get the number of resonance structures generated so far for the passed conjugated group index.

        C++ signature :unsigned long GetNumStructures(RDKit::PyResonanceMolSupplierCallback {lvalue},unsigned int)
        """
        ...
    @overload
    @classmethod
    def __call__(cls, RDKit) -> Any: ...
    @overload
    @classmethod
    def __call__(cls, RDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class RingInfo(Boost.Python.instance):
    """
    contains information about a molecule’s rings
    Raises an exception
    This class cannot be instantiated from Python"""

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def AddRing(
        self: RingInfo, atomIds: AtomPairsParameters, bondIds: AtomPairsParameters
    ) -> None:
        """
        Adds a ring to the set. Be very careful with this operation.

        C++ signature :void AddRing(RDKit::RingInfo*,boost::python::api::object,boost::python::api::object)
        """
        ...
    def AreAtomsInSameRing(self, arg1: RingInfo, arg2: int, arg3: int) -> bool:
        """
        C++ signature :bool AreAtomsInSameRing(RDKit::RingInfo {lvalue},unsigned int,unsigned int)
        """
        ...
    def AreAtomsInSameRingOfSize(
        self, arg1: RingInfo, arg2: int, arg3: int, arg4: int
    ) -> bool:
        """
        C++ signature :bool AreAtomsInSameRingOfSize(RDKit::RingInfo {lvalue},unsigned int,unsigned int,unsigned int)
        """
        ...
    def AreBondsInSameRing(self, arg1: RingInfo, arg2: int, arg3: int) -> bool:
        """
        C++ signature :bool AreBondsInSameRing(RDKit::RingInfo {lvalue},unsigned int,unsigned int)
        """
        ...
    def AreBondsInSameRingOfSize(
        self, arg1: RingInfo, arg2: int, arg3: int, arg4: int
    ) -> bool:
        """
        C++ signature :bool AreBondsInSameRingOfSize(RDKit::RingInfo {lvalue},unsigned int,unsigned int,unsigned int)
        """
        ...
    def AreRingFamiliesInitialized(self, arg1: RingInfo) -> bool:
        """
        C++ signature :bool AreRingFamiliesInitialized(RDKit::RingInfo {lvalue})"""
        ...
    def AtomMembers(self, arg1: RingInfo, arg2: int) -> object:
        """
        C++ signature :boost::python::api::object AtomMembers(RDKit::RingInfo const*,unsigned int)
        """
        ...
    def AtomRingFamilies(self, arg1: RingInfo) -> object:
        """
        C++ signature :boost::python::api::object AtomRingFamilies(RDKit::RingInfo const*)
        """
        ...
    def AtomRingSizes(self, arg1: RingInfo, arg2: int) -> object:
        """
        C++ signature :boost::python::api::object AtomRingSizes(RDKit::RingInfo const*,unsigned int)
        """
        ...
    def AtomRings(self, arg1: RingInfo) -> object:
        """
        C++ signature :boost::python::api::object AtomRings(RDKit::RingInfo const*)"""
        ...
    def BondMembers(self, arg1: RingInfo, arg2: int) -> object:
        """
        C++ signature :boost::python::api::object BondMembers(RDKit::RingInfo const*,unsigned int)
        """
        ...
    def BondRingFamilies(self, arg1: RingInfo) -> object:
        """
        C++ signature :boost::python::api::object BondRingFamilies(RDKit::RingInfo const*)
        """
        ...
    def BondRingSizes(self, arg1: RingInfo, arg2: int) -> object:
        """
        C++ signature :boost::python::api::object BondRingSizes(RDKit::RingInfo const*,unsigned int)
        """
        ...
    def BondRings(self, arg1: RingInfo) -> object:
        """
        C++ signature :boost::python::api::object BondRings(RDKit::RingInfo const*)"""
        ...
    def IsAtomInRingOfSize(self, arg1: RingInfo, arg2: int, arg3: int) -> bool:
        """
        C++ signature :bool IsAtomInRingOfSize(RDKit::RingInfo {lvalue},unsigned int,unsigned int)
        """
        ...
    def IsBondInRingOfSize(self, arg1: RingInfo, arg2: int, arg3: int) -> bool:
        """
        C++ signature :bool IsBondInRingOfSize(RDKit::RingInfo {lvalue},unsigned int,unsigned int)
        """
        ...
    def MinAtomRingSize(self, arg1: RingInfo, arg2: int) -> int:
        """
        C++ signature :unsigned int MinAtomRingSize(RDKit::RingInfo {lvalue},unsigned int)
        """
        ...
    def MinBondRingSize(self, arg1: RingInfo, arg2: int) -> int:
        """
        C++ signature :unsigned int MinBondRingSize(RDKit::RingInfo {lvalue},unsigned int)
        """
        ...
    def NumAtomRings(self, arg1: RingInfo, arg2: int) -> int:
        """
        C++ signature :unsigned int NumAtomRings(RDKit::RingInfo {lvalue},unsigned int)
        """
        ...
    def NumBondRings(self, arg1: RingInfo, arg2: int) -> int:
        """
        C++ signature :unsigned int NumBondRings(RDKit::RingInfo {lvalue},unsigned int)
        """
        ...
    def NumRelevantCycles(self, arg1: RingInfo) -> int:
        """
        C++ signature :unsigned int NumRelevantCycles(RDKit::RingInfo {lvalue})"""
        ...
    def NumRingFamilies(self, arg1: RingInfo) -> int:
        """
        C++ signature :unsigned int NumRingFamilies(RDKit::RingInfo {lvalue})"""
        ...
    def NumRings(self, arg1: RingInfo) -> int:
        """
        C++ signature :unsigned int NumRings(RDKit::RingInfo {lvalue})"""
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class StereoDescriptor(Boost.Python.enum):
    Bond_Cis: StereoDescriptor = ...
    Bond_Trans: StereoDescriptor = ...
    NoValue: StereoDescriptor = ...
    Tet_CCW: StereoDescriptor = ...
    Tet_CW: StereoDescriptor = ...
    names: dict[str, StereoDescriptor] = ...
    values: dict[int, StereoDescriptor] = ...
    __slots__: ClassVar[tuple] = ...

class StereoGroup(Boost.Python.instance):
    """
    A collection of atoms with a defined stereochemical relationship.
    Used to help represent a sample with unknown stereochemistry, or that is a mix
    of diastereomers.
    Raises an exception
    This class cannot be instantiated from Python"""

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def GetAtoms(self, arg1: StereoGroup) -> object:
        """
        access the atoms in the StereoGroup.

        C++ signature :boost::python::api::object GetAtoms(RDKit::StereoGroup {lvalue})
        """
        ...
    def GetGroupType(self, arg1: StereoGroup) -> StereoGroupType:
        """
        Returns the StereoGroupType.

        C++ signature :RDKit::StereoGroupType GetGroupType(RDKit::StereoGroup {lvalue})
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class StereoGroupType(Boost.Python.enum):
    STEREO_ABSOLUTE: StereoGroupType = ...
    STEREO_AND: StereoGroupType = ...
    STEREO_OR: StereoGroupType = ...
    names: dict[str, StereoGroupType] = ...
    values: dict[int, StereoGroupType] = ...
    __slots__: ClassVar[tuple] = ...

class StereoGroup_vect(Boost.Python.instance):
    """
    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def append(self, arg1: StereoGroup_vect, arg2: AtomPairsParameters) -> None:
        """
        C++ signature :void append(std::vector<RDKit::StereoGroup, std::allocator<RDKit::StereoGroup> > {lvalue},boost::python::api::object)
        """
        ...
    def extend(self, arg1: StereoGroup_vect, arg2: AtomPairsParameters) -> None:
        """
        C++ signature :void extend(std::vector<RDKit::StereoGroup, std::allocator<RDKit::StereoGroup> > {lvalue},boost::python::api::object)
        """
        ...
    @classmethod
    def __contains__(cls, other) -> Any: ...
    @classmethod
    def __delitem__(cls, other) -> Any: ...
    @classmethod
    def __getitem__(cls, index) -> Any: ...
    @classmethod
    def __iter__(cls, boost, std) -> Any: ...
    @classmethod
    def __len__(cls) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @classmethod
    def __setitem__(cls, index, object) -> Any: ...

class StereoInfo(Boost.Python.instance):
    """
    Class describing stereochemistry

    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...
    centeredOn: Any
    descriptor: Any
    permutation: Any
    specified: Any
    type: Any

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    NOATOM: property[int] = ...

    @property
    def controllingAtoms(self) -> Any: ...

class StereoSpecified(Boost.Python.enum):
    Specified: StereoSpecified = ...
    Unknown: StereoSpecified = ...
    Unspecified: StereoSpecified = ...
    names: dict[str, StereoSpecified] = ...
    values: dict[int, StereoSpecified] = ...
    __slots__: ClassVar[tuple] = ...

class StereoType(Boost.Python.enum):
    Atom_Octahedral: StereoType = ...
    Atom_SquarePlanar: StereoType = ...
    Atom_Tetrahedral: StereoType = ...
    Atom_TrigonalBipyramidal: StereoType = ...
    Bond_Atropisomer: StereoType = ...
    Bond_Cumulene_Even: StereoType = ...
    Bond_Double: StereoType = ...
    Unspecified: StereoType = ...
    names: dict[str, StereoType] = ...
    values: dict[int, StereoType] = ...
    __slots__: ClassVar[tuple] = ...

class SubstanceGroup(Boost.Python.instance):
    """
    A collection of atoms and bonds with associated properties
    Raises an exception
    This class cannot be instantiated from Python"""

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def AddAtomWithBookmark(self, arg1: SubstanceGroup, arg2: int) -> None:
        """
        C++ signature :void AddAtomWithBookmark(RDKit::SubstanceGroup {lvalue},int)"""
        ...
    def AddAtomWithIdx(self, arg1: SubstanceGroup, arg2: int) -> None:
        """
        C++ signature :void AddAtomWithIdx(RDKit::SubstanceGroup {lvalue},unsigned int)
        """
        ...
    def AddAttachPoint(
        self, arg1: SubstanceGroup, arg2: int, arg3: int, arg4: str
    ) -> None:
        """
        C++ signature :void AddAttachPoint(RDKit::SubstanceGroup {lvalue},unsigned int,int,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    def AddBondWithBookmark(self, arg1: SubstanceGroup, arg2: int) -> None:
        """
        C++ signature :void AddBondWithBookmark(RDKit::SubstanceGroup {lvalue},int)"""
        ...
    def AddBondWithIdx(self, arg1: SubstanceGroup, arg2: int) -> None:
        """
        C++ signature :void AddBondWithIdx(RDKit::SubstanceGroup {lvalue},unsigned int)
        """
        ...
    def AddBracket(self, arg1: SubstanceGroup, arg2: AtomPairsParameters) -> None:
        """
        C++ signature :void AddBracket(RDKit::SubstanceGroup {lvalue},boost::python::api::object)
        """
        ...
    def AddCState(self, arg1: SubstanceGroup, arg2: int, arg3: Point3D) -> None:
        """
        C++ signature :void AddCState(RDKit::SubstanceGroup {lvalue},unsigned int,RDGeom::Point3D)
        """
        ...
    def AddParentAtomWithBookmark(self, arg1: SubstanceGroup, arg2: int) -> None:
        """
        C++ signature :void AddParentAtomWithBookmark(RDKit::SubstanceGroup {lvalue},int)
        """
        ...
    def AddParentAtomWithIdx(self, arg1: SubstanceGroup, arg2: int) -> None:
        """
        C++ signature :void AddParentAtomWithIdx(RDKit::SubstanceGroup {lvalue},unsigned int)
        """
        ...
    def ClearAttachPoints(self, arg1: SubstanceGroup) -> None:
        """
        C++ signature :void ClearAttachPoints(RDKit::SubstanceGroup {lvalue})"""
        ...
    def ClearBrackets(self, arg1: SubstanceGroup) -> None:
        """
        C++ signature :void ClearBrackets(RDKit::SubstanceGroup {lvalue})"""
        ...
    def ClearCStates(self, arg1: SubstanceGroup) -> None:
        """
        C++ signature :void ClearCStates(RDKit::SubstanceGroup {lvalue})"""
        ...
    def ClearProp(self, arg1: SubstanceGroup, arg2: str) -> None:
        """
        Removes a particular property (does nothing if not set).

        C++ signature :void ClearProp(RDKit::SubstanceGroup {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    def GetAtoms(self, arg1: SubstanceGroup) -> _vectj:
        """
        returns a list of the indices of the atoms in this SubstanceGroup

        C++ signature :std::vector<unsigned int, std::allocator<unsigned int> > GetAtoms(RDKit::SubstanceGroup {lvalue})
        """
        ...
    def GetAttachPoints(self, arg1: SubstanceGroup) -> tuple:
        """
        C++ signature :boost::python::tuple GetAttachPoints(RDKit::SubstanceGroup)"""
        ...
    def GetBonds(self, arg1: SubstanceGroup) -> _vectj:
        """
        returns a list of the indices of the bonds in this SubstanceGroup

        C++ signature :std::vector<unsigned int, std::allocator<unsigned int> > GetBonds(RDKit::SubstanceGroup {lvalue})
        """
        ...
    def GetBoolProp(self, arg1: SubstanceGroup, arg2: str) -> bool:
        """
        returns the value of a particular property

        C++ signature :bool GetBoolProp(RDKit::SubstanceGroup {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    def GetBrackets(self, arg1: SubstanceGroup) -> tuple:
        """
        C++ signature :boost::python::tuple GetBrackets(RDKit::SubstanceGroup)"""
        ...
    def GetCStates(self, arg1: SubstanceGroup) -> tuple:
        """
        C++ signature :boost::python::tuple GetCStates(RDKit::SubstanceGroup)"""
        ...
    def GetDoubleProp(self, arg1: SubstanceGroup, arg2: str) -> float:
        """
        returns the value of a particular property

        C++ signature :double GetDoubleProp(RDKit::SubstanceGroup {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    def GetIndexInMol(self, arg1: SubstanceGroup) -> int:
        """
        returns the index of this SubstanceGroup in the owning molecule’s list.

        C++ signature :unsigned int GetIndexInMol(RDKit::SubstanceGroup {lvalue})"""
        ...
    def GetIntProp(self, arg1: SubstanceGroup, arg2: str) -> int:
        """
        returns the value of a particular property

        C++ signature :int GetIntProp(RDKit::SubstanceGroup {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    def GetOwningMol(self, arg1: SubstanceGroup) -> Mol:
        """
        returns the molecule owning this SubstanceGroup

        C++ signature :RDKit::ROMol {lvalue} GetOwningMol(RDKit::SubstanceGroup {lvalue})
        """
        ...
    def GetParentAtoms(self, arg1: SubstanceGroup) -> _vectj:
        """
        returns a list of the indices of the parent atoms in this SubstanceGroup

        C++ signature :std::vector<unsigned int, std::allocator<unsigned int> > GetParentAtoms(RDKit::SubstanceGroup {lvalue})
        """
        ...
    def GetProp(self, arg1: SubstanceGroup, arg2: str) -> str:
        """
        returns the value of a particular property

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetProp(RDKit::SubstanceGroup {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    def GetPropNames(
        self: SubstanceGroup,
        includePrivate: bool = False,
        includeComputed: bool = False,
    ) -> _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE:
        """
        Returns a list of the properties set on the SubstanceGroup.

        C++ signature :std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > GetPropNames(RDKit::SubstanceGroup {lvalue} [,bool=False [,bool=False]])
        """
        ...
    def GetPropsAsDict(
        self: SubstanceGroup, includePrivate: bool = True, includeComputed: bool = True
    ) -> dict:
        """
        Returns a dictionary of the properties set on the SubstanceGroup.n.b. some properties cannot be converted to python types.

        C++ signature :boost::python::dict GetPropsAsDict(RDKit::SubstanceGroup [,bool=True [,bool=True]])
        """
        ...
    def GetStringVectProp(
        self, arg1: SubstanceGroup, arg2: str
    ) -> _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE:
        """
        returns the value of a particular property

        C++ signature :std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > GetStringVectProp(RDKit::SubstanceGroup {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    def GetUnsignedProp(self, arg1: SubstanceGroup, arg2: str) -> int:
        """
        returns the value of a particular property

        C++ signature :unsigned int GetUnsignedProp(RDKit::SubstanceGroup {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    def GetUnsignedVectProp(self, arg1: SubstanceGroup, arg2: str) -> _vectj:
        """
        returns the value of a particular property

        C++ signature :std::vector<unsigned int, std::allocator<unsigned int> > GetUnsignedVectProp(RDKit::SubstanceGroup {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    def HasProp(self, arg1: SubstanceGroup, arg2: str) -> bool:
        """
        returns whether or not a particular property exists

        C++ signature :bool HasProp(RDKit::SubstanceGroup {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    def SetAtoms(self, arg1: SubstanceGroup, arg2: AtomPairsParameters) -> None:
        """
        Set the list of the indices of the atoms in this SubstanceGroup.
        Note that this does not update properties, CStates or Attachment Points.

        C++ signature :void SetAtoms(RDKit::SubstanceGroup {lvalue},boost::python::api::object)
        """
        ...
    def SetBonds(self, arg1: SubstanceGroup, arg2: AtomPairsParameters) -> None:
        """
        Set the list of the indices of the bonds in this SubstanceGroup.
        Note that this does not update properties, CStates or Attachment Points.

        C++ signature :void SetBonds(RDKit::SubstanceGroup {lvalue},boost::python::api::object)
        """
        ...
    def SetBoolProp(
        self: SubstanceGroup, key: str, val: bool, computed: bool = False
    ) -> None:
        """
        sets the value of a particular property

        C++ signature :void SetBoolProp(RDKit::SubstanceGroup {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,bool [,bool=False])
        """
        ...
    def SetDoubleProp(
        self: SubstanceGroup, key: str, val: float, computed: bool = False
    ) -> None:
        """
        sets the value of a particular property

        C++ signature :void SetDoubleProp(RDKit::SubstanceGroup {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,double [,bool=False])
        """
        ...
    def SetIntProp(
        self: SubstanceGroup, key: str, val: int, computed: bool = False
    ) -> None:
        """
        sets the value of a particular property

        C++ signature :void SetIntProp(RDKit::SubstanceGroup {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,int [,bool=False])
        """
        ...
    def SetParentAtoms(self, arg1: SubstanceGroup, arg2: AtomPairsParameters) -> None:
        """
        Set the list of the indices of the parent atoms in this SubstanceGroup.
        Note that this does not update properties, CStates or Attachment Points.

        C++ signature :void SetParentAtoms(RDKit::SubstanceGroup {lvalue},boost::python::api::object)
        """
        ...
    def SetProp(
        self: SubstanceGroup, key: str, val: str, computed: bool = False
    ) -> None:
        """
        sets the value of a particular property

        C++ signature :void SetProp(RDKit::SubstanceGroup {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])
        """
        ...
    def SetUnsignedProp(
        self: SubstanceGroup, key: str, val: int, computed: bool = False
    ) -> None:
        """
        sets the value of a particular property

        C++ signature :void SetUnsignedProp(RDKit::SubstanceGroup {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,unsigned int [,bool=False])
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class SubstanceGroupAttach(Boost.Python.instance):
    """
    AttachPoint for a SubstanceGroup

    C++ signature :void __init__(_object*)

    property aIdx¶
    attachment index

    property id¶
    attachment id

    property lvIdx¶
    leaving atom or index (0 for implied)"""

    ...
    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @property
    def aIdx(self) -> Any: ...
    @property
    def id(self) -> Any: ...
    @property
    def lvIdx(self) -> Any: ...

class SubstanceGroupCState(Boost.Python.instance):
    """
    CSTATE for a SubstanceGroup

    C++ signature :void __init__(_object*)

    property bondIdx¶

    property vector¶"""

    ...
    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @property
    def bondIdx(self) -> Any: ...
    @property
    def vector(self) -> Any: ...

class SubstanceGroup_VECT(Boost.Python.instance):
    """
    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def append(self, arg1: SubstanceGroup_VECT, arg2: AtomPairsParameters) -> None:
        """
        C++ signature :void append(std::vector<RDKit::SubstanceGroup, std::allocator<RDKit::SubstanceGroup> > {lvalue},boost::python::api::object)
        """
        ...
    def extend(self, arg1: SubstanceGroup_VECT, arg2: AtomPairsParameters) -> None:
        """
        C++ signature :void extend(std::vector<RDKit::SubstanceGroup, std::allocator<RDKit::SubstanceGroup> > {lvalue},boost::python::api::object)
        """
        ...
    @classmethod
    def __contains__(cls, other) -> Any: ...
    @classmethod
    def __delitem__(cls, other) -> Any: ...
    @classmethod
    def __getitem__(cls, index) -> Any: ...
    @classmethod
    def __iter__(cls, boost, std) -> Any: ...
    @classmethod
    def __len__(cls) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @classmethod
    def __setitem__(cls, index, object) -> Any: ...

class SubstructMatchParameters(Boost.Python.instance):
    """
    Parameters controlling substructure matching

    C++ signature :void __init__(_object*)

    property aromaticMatchesConjugated¶
    aromatic and conjugated bonds match each other

    property maxMatches¶
    maximum number of matches to return

    property numThreads¶
    number of threads to use when multi-threading is possible.0 selects the number of concurrent threads supported by thehardware. negative values are added to the number of concurrentthreads supported by the hardware.

    property recursionPossible¶
    Allow recursive queries"""

    __instance_size__: ClassVar[int] = ...
    aromaticMatchesConjugated: Any
    maxMatches: Any
    numThreads: Any
    recursionPossible: Any
    uniquify: Any
    useChirality: Any
    useEnhancedStereo: Any
    useGenericMatchers: Any
    useQueryQueryMatches: Any

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def setExtraFinalCheck(
        self, arg1: SubstructMatchParameters, arg2: AtomPairsParameters
    ) -> None:
        """
        allows you to provide a function that will be called
        with the molecule

        and a vector of atom IDs containing a potential match.
        The function should return true or false indicating whether or not
        that match should be accepted.

        C++ signature :void setExtraFinalCheck(RDKit::SubstructMatchParameters {lvalue},boost::python::api::object)
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class AtomKekulizeException(AtomSanitizeException):
    ...
    ...

class AtomSanitizeException(MolSanitizeException):
    ...
    ...

class AtomValenceException(AtomSanitizeException):
    ...
    ...

class KekulizeException(MolSanitizeException):
    ...
    ...

class MolSanitizeException(ValueError):
    ...
    ...

class _ROAtomSeq(Boost.Python.instance):
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __getitem__(cls, index) -> Any: ...
    @classmethod
    def __iter__(cls) -> Any: ...
    @classmethod
    def __len__(cls) -> Any: ...
    @classmethod
    def __next__(cls) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class _ROBondSeq(Boost.Python.instance):
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __getitem__(cls, index) -> Any: ...
    @classmethod
    def __iter__(cls) -> Any: ...
    @classmethod
    def __len__(cls) -> Any: ...
    @classmethod
    def __next__(cls) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class _ROConformerSeq(Boost.Python.instance):
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __getitem__(cls, index) -> Any: ...
    @classmethod
    def __iter__(cls) -> Any: ...
    @classmethod
    def __len__(cls) -> Any: ...
    @classmethod
    def __next__(cls) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class _ROQAtomSeq(Boost.Python.instance):
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __getitem__(cls, index) -> Any: ...
    @classmethod
    def __iter__(cls) -> Any: ...
    @classmethod
    def __len__(cls) -> Any: ...
    @classmethod
    def __next__(cls) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class _cppAtomKekulizeException(_cppMolSanitizeException):
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def GetAtomIndices(cls, RDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class _cppAtomSanitizeException(_cppMolSanitizeException):
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def GetAtomIdx(cls, RDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class _cppAtomValenceException(_cppAtomSanitizeException):
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class _cppMolSanitizeException(Boost.Python.instance):
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def GetType(cls, RDKit) -> Any: ...
    @classmethod
    def Message(cls, RDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class _listN5boost10shared_ptrIN5RDKit9ConformerEEE(Boost.Python.instance):
    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __contains__(cls, other) -> Any: ...
    @classmethod
    def __delitem__(cls, other) -> Any: ...
    @classmethod
    def __getitem__(cls, index) -> Any: ...
    @classmethod
    def __iter__(cls, boost, std) -> Any: ...
    @classmethod
    def __len__(cls) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @classmethod
    def __setitem__(cls, index, object) -> Any: ...

class _listPN5RDKit4AtomE(Boost.Python.instance):
    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __contains__(cls, other) -> Any: ...
    @classmethod
    def __delitem__(cls, other) -> Any: ...
    @classmethod
    def __getitem__(cls, index) -> Any: ...
    @classmethod
    def __iter__(cls, boost, std) -> Any: ...
    @classmethod
    def __len__(cls) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @classmethod
    def __setitem__(cls, index, object) -> Any: ...

class _listPN5RDKit4BondE(Boost.Python.instance):
    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __contains__(cls, other) -> Any: ...
    @classmethod
    def __delitem__(cls, other) -> Any: ...
    @classmethod
    def __getitem__(cls, index) -> Any: ...
    @classmethod
    def __iter__(cls, boost, std) -> Any: ...
    @classmethod
    def __len__(cls) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @classmethod
    def __setitem__(cls, index, object) -> Any: ...

def AddMolSubstanceGroup(self, mol: Mol, sgroup: SubstanceGroup) -> SubstanceGroup:
    """
    adds a copy of a SubstanceGroup to a molecule, returns the new SubstanceGroup

    C++ signature :RDKit::SubstanceGroup* AddMolSubstanceGroup(RDKit::ROMol {lvalue},RDKit::SubstanceGroup)
    """
    ...

def ClearMolSubstanceGroups(self, arg1: Mol) -> None:
    """
    removes all SubstanceGroups from a molecule (if any)

    C++ signature :void ClearMolSubstanceGroups(RDKit::ROMol {lvalue})"""
    ...

def CreateMolDataSubstanceGroup(
    self, mol: Mol, fieldName: str, value: str
) -> SubstanceGroup:
    """
    creates a new DATA SubstanceGroup associated with a molecule, returns the new SubstanceGroup

    C++ signature :RDKit::SubstanceGroup* CreateMolDataSubstanceGroup(RDKit::ROMol {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
    ...

def CreateMolSubstanceGroup(self, mol: Mol, type: str) -> SubstanceGroup:
    """
    creates a new SubstanceGroup associated with a molecule, returns the new SubstanceGroup

    C++ signature :RDKit::SubstanceGroup* CreateMolSubstanceGroup(RDKit::ROMol {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
    ...

def CreateStereoGroup(
    self, stereoGroupType: StereoGroupType, mol: Mol, atomIds: AtomPairsParameters
) -> StereoGroup:
    """
    creates a StereoGroup associated with a molecule from a list of atom Ids

    C++ signature :RDKit::StereoGroup* CreateStereoGroup(RDKit::StereoGroupType,RDKit::ROMol {lvalue},boost::python::api::object)
    """
    ...

def GetAtomAlias(self, atom: Atom) -> str:
    """
    Returns the atom’s MDL alias text

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetAtomAlias(RDKit::Atom const*)
    """
    ...

def GetAtomRLabel(self, atom: Atom) -> int:
    """
    Returns the atom’s MDL AtomRLabel (this is an integer from 0 to 99)

    C++ signature :int GetAtomRLabel(RDKit::Atom const*)"""
    ...

def GetAtomValue(self, atom: Atom) -> str:
    """
    Returns the atom’s MDL alias text

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetAtomValue(RDKit::Atom const*)
    """
    ...

def GetDefaultPickleProperties(self) -> int:
    """
    Get the current global mol pickler options.

    C++ signature :unsigned int GetDefaultPickleProperties()"""
    ...

def GetMolSubstanceGroupWithIdx(self, arg1: Mol, arg2: int) -> SubstanceGroup:
    """
    returns a particular SubstanceGroup from the molecule

    C++ signature :RDKit::SubstanceGroup* GetMolSubstanceGroupWithIdx(RDKit::ROMol {lvalue},unsigned int)
    """
    ...

def GetMolSubstanceGroups(self, arg1: Mol) -> SubstanceGroup_VECT:
    """
    returns a copy of the molecule’s SubstanceGroups (if any)

    C++ signature :std::vector<RDKit::SubstanceGroup, std::allocator<RDKit::SubstanceGroup> > GetMolSubstanceGroups(RDKit::ROMol {lvalue})
    """
    ...

def GetPeriodicTable(self) -> PeriodicTable:
    """
    Returns the application’s PeriodicTable instance.

    C++ signature :RDKit::PeriodicTable* GetPeriodicTable()"""
    ...

def GetSupplementalSmilesLabel(self, atom: Atom) -> str:
    """
    Gets the supplemental smiles label on an atom, returns an empty string if not present.

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetSupplementalSmilesLabel(RDKit::Atom const*)
    """
    ...

def SetAtomAlias(self, atom: Atom, rlabel: str) -> None:
    """
    Sets the atom’s MDL alias text.
    Setting to an empty string clears the alias.

    C++ signature :void SetAtomAlias(RDKit::Atom*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
    ...

def SetAtomRLabel(self, atom: Atom, rlabel: int) -> None:
    """
    Sets the atom’s MDL RLabel (this is an integer from 0 to 99).
    Setting to 0 clears the rlabel.

    C++ signature :void SetAtomRLabel(RDKit::Atom*,int)"""
    ...

def SetAtomValue(self, atom: Atom, rlabel: str) -> None:
    """
    Sets the atom’s MDL alias text.
    Setting to an empty string clears the alias.

    C++ signature :void SetAtomValue(RDKit::Atom*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
    ...

def SetDefaultPickleProperties(self, arg1: int) -> None:
    """
    Set the current global mol pickler options.

    C++ signature :void SetDefaultPickleProperties(unsigned int)"""
    ...

def SetSupplementalSmilesLabel(self, atom: Atom, label: str) -> None:
    """
    Sets a supplemental label on an atom that is written to the smiles string.
    >>> m = Chem.MolFromSmiles("C")
    >>> Chem.SetSupplementalSmilesLabel(m.GetAtomWithIdx(0), '<xxx>')
    >>> Chem.MolToSmiles(m)
    'C<xxx>'

    C++ signature :void SetSupplementalSmilesLabel(RDKit::Atom*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
    ...

def _HasSubstructMatchStr(std, RDKit) -> Any: ...
def tossit(self) -> None:
    """
    C++ signature :void tossit()"""
    ...
