"""
rdkit.Chem.rdqueries module¶
Module containing RDKit functionality for querying molecules.
"""
from typing import Any, overload

from rdkit.Chem.rdchem import QueryAtom, QueryBond
from rdkit.DataStructs.cDataStructs import ExplicitBitVect

def AAtomQueryAtom(self, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms when AAtom is True.

    C++ signature :RDKit::QueryAtom* AAtomQueryAtom([ bool=False])"""
    ...

def AHAtomQueryAtom(self, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms when AHAtom is True.

    C++ signature :RDKit::QueryAtom* AHAtomQueryAtom([ bool=False])"""
    ...

def AtomNumEqualsQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where AtomNum is equal to the target value.

    C++ signature :RDKit::QueryAtom* AtomNumEqualsQueryAtom(int [,bool=False])"""
    ...

def AtomNumGreaterQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where AtomNum is equal to the target value.
    NOTE: the direction of comparison is reversed relative to the C++ API

    C++ signature :RDKit::QueryAtom* AtomNumGreaterQueryAtom(int [,bool=False])"""
    ...

def AtomNumLessQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where AtomNum is less than the target value.
    NOTE: the direction of comparison is reversed relative to the C++ API

    C++ signature :RDKit::QueryAtom* AtomNumLessQueryAtom(int [,bool=False])"""
    ...

def ExplicitDegreeEqualsQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where ExplicitDegree is equal to the target value.

    C++ signature :RDKit::QueryAtom* ExplicitDegreeEqualsQueryAtom(int [,bool=False])"""
    ...

def ExplicitDegreeGreaterQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where ExplicitDegree is equal to the target value.
    NOTE: the direction of comparison is reversed relative to the C++ API

    C++ signature :RDKit::QueryAtom* ExplicitDegreeGreaterQueryAtom(int [,bool=False])
    """
    ...

def ExplicitDegreeLessQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where ExplicitDegree is less than the target value.
    NOTE: the direction of comparison is reversed relative to the C++ API

    C++ signature :RDKit::QueryAtom* ExplicitDegreeLessQueryAtom(int [,bool=False])"""
    ...

def ExplicitValenceEqualsQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where ExplicitValence is equal to the target value.

    C++ signature :RDKit::QueryAtom* ExplicitValenceEqualsQueryAtom(int [,bool=False])
    """
    ...

def ExplicitValenceGreaterQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where ExplicitValence is equal to the target value.
    NOTE: the direction of comparison is reversed relative to the C++ API

    C++ signature :RDKit::QueryAtom* ExplicitValenceGreaterQueryAtom(int [,bool=False])
    """
    ...

def ExplicitValenceLessQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where ExplicitValence is less than the target value.
    NOTE: the direction of comparison is reversed relative to the C++ API

    C++ signature :RDKit::QueryAtom* ExplicitValenceLessQueryAtom(int [,bool=False])"""
    ...

def FormalChargeEqualsQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where FormalCharge is equal to the target value.

    C++ signature :RDKit::QueryAtom* FormalChargeEqualsQueryAtom(int [,bool=False])"""
    ...

def FormalChargeGreaterQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where FormalCharge is equal to the target value.
    NOTE: the direction of comparison is reversed relative to the C++ API

    C++ signature :RDKit::QueryAtom* FormalChargeGreaterQueryAtom(int [,bool=False])"""
    ...

def FormalChargeLessQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where FormalCharge is less than the target value.
    NOTE: the direction of comparison is reversed relative to the C++ API

    C++ signature :RDKit::QueryAtom* FormalChargeLessQueryAtom(int [,bool=False])"""
    ...

def HCountEqualsQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where HCount is equal to the target value.

    C++ signature :RDKit::QueryAtom* HCountEqualsQueryAtom(int [,bool=False])"""
    ...

def HCountGreaterQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where HCount is equal to the target value.
    NOTE: the direction of comparison is reversed relative to the C++ API

    C++ signature :RDKit::QueryAtom* HCountGreaterQueryAtom(int [,bool=False])"""
    ...

def HCountLessQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where HCount is less than the target value.
    NOTE: the direction of comparison is reversed relative to the C++ API

    C++ signature :RDKit::QueryAtom* HCountLessQueryAtom(int [,bool=False])"""
    ...

def HasBitVectPropWithValueQueryAtom(
    self,
    propname: str,
    val: ExplicitBitVect,
    negate: bool = False,
    tolerance: float = 0,
) -> QueryAtom:
    """
    Returns a QueryAtom that matches when the propery ‘propname’ has the specified explicit bit vector value.  The Tolerance is the allowed Tanimoto difference

    C++ signature :RDKit::QueryAtom* HasBitVectPropWithValueQueryAtom(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,ExplicitBitVect [,bool=False [,float=0]])
    """
    ...

def HasBoolPropWithValueQueryAtom(
    self, propname: str, val: bool, negate: bool = False
) -> QueryAtom:
    """
    Returns a QueryAtom that matches when the propery ‘propname’ has the specified boolean value.

    C++ signature :RDKit::QueryAtom* HasBoolPropWithValueQueryAtom(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,bool [,bool=False])
    """
    ...

def HasBoolPropWithValueQueryBond(
    self, propname: str, val: bool, negate: bool = False
) -> QueryBond:
    """
    Returns a QueryBond that matches when the propery ‘propname’ has the specified boolean value.

    C++ signature :RDKit::QueryBond* HasBoolPropWithValueQueryBond(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,bool [,bool=False])
    """
    ...

def HasChiralTagQueryAtom(self, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms when HasChiralTag is True.

    C++ signature :RDKit::QueryAtom* HasChiralTagQueryAtom([ bool=False])"""
    ...

def HasDoublePropWithValueQueryAtom(
    self, propname: str, val: float, negate: bool = False, tolerance: float = 0.0
) -> QueryAtom:
    """
    Returns a QueryAtom that matches when the propery ‘propname’ has the specified value +- tolerance

    C++ signature :RDKit::QueryAtom* HasDoublePropWithValueQueryAtom(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,double [,bool=False [,double=0.0]])
    """
    ...

def HasDoublePropWithValueQueryBond(
    self, propname: str, val: float, negate: bool = False, tolerance: float = 0.0
) -> QueryBond:
    """
    Returns a QueryBond that matches when the propery ‘propname’ has the specified value +- tolerance

    C++ signature :RDKit::QueryBond* HasDoublePropWithValueQueryBond(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,double [,bool=False [,double=0.0]])
    """
    ...

def HasIntPropWithValueQueryAtom(
    self, propname: str, val: int, negate: bool = False, tolerance: int = 0
) -> QueryAtom:
    """
    Returns a QueryAtom that matches when the propery ‘propname’ has the specified int value.

    C++ signature :RDKit::QueryAtom* HasIntPropWithValueQueryAtom(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,int [,bool=False [,int=0]])
    """
    ...

def HasIntPropWithValueQueryBond(
    self, propname: str, val: int, negate: bool = False, tolerance: int = 0
) -> QueryBond:
    """
    Returns a QueryBond that matches when the propery ‘propname’ has the specified int value.

    C++ signature :RDKit::QueryBond* HasIntPropWithValueQueryBond(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,int [,bool=False [,int=0]])
    """
    ...

def HasPropQueryAtom(self, propname: str, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches when the propery ‘propname’ exists in the atom.

    C++ signature :RDKit::QueryAtom* HasPropQueryAtom(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])
    """
    ...

@overload
def HasPropQueryBond(self, propname: str, negate: bool = False) -> QueryBond:
    """
    Returns a QueryBond that matches when the propery ‘propname’ exists in the bond.

    C++ signature :RDKit::QueryBond* HasPropQueryBond(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])

    HasPropQueryBond( (str)propname [, (bool)negate=False]) -> QueryBond :Returns a QueryBond that matches when the propery ‘propname’ exists in the bond.

    C++ signature :RDKit::QueryBond* HasPropQueryBond(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])
    """
    ...

@overload
def HasPropQueryBond(self, propname: str, negate: bool = False) -> QueryBond:
    """
    Returns a QueryBond that matches when the propery ‘propname’ exists in the bond.

    C++ signature :RDKit::QueryBond* HasPropQueryBond(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])
    """
    ...

@overload
def HasPropQueryBond(self, propname: str, negate: bool = False) -> QueryBond: ...
def HasStringPropWithValueQueryAtom(
    self, propname: str, val: str, negate: bool = False
) -> QueryAtom:
    """
    Returns a QueryAtom that matches when the propery ‘propname’ has the specified string value.

    C++ signature :RDKit::QueryAtom* HasStringPropWithValueQueryAtom(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])
    """
    ...

def HasStringPropWithValueQueryBond(
    self, propname: str, val: str, negate: bool = False
) -> QueryBond:
    """
    Returns a QueryBond that matches when the propery ‘propname’ has the specified string value.

    C++ signature :RDKit::QueryBond* HasStringPropWithValueQueryBond(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])
    """
    ...

def HybridizationEqualsQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where Hybridization is equal to the target value.

    C++ signature :RDKit::QueryAtom* HybridizationEqualsQueryAtom(int [,bool=False])"""
    ...

def HybridizationGreaterQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where Hybridization is equal to the target value.
    NOTE: the direction of comparison is reversed relative to the C++ API

    C++ signature :RDKit::QueryAtom* HybridizationGreaterQueryAtom(int [,bool=False])"""
    ...

def HybridizationLessQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where Hybridization is less than the target value.
    NOTE: the direction of comparison is reversed relative to the C++ API

    C++ signature :RDKit::QueryAtom* HybridizationLessQueryAtom(int [,bool=False])"""
    ...

def InNRingsEqualsQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where InNRings is equal to the target value.

    C++ signature :RDKit::QueryAtom* InNRingsEqualsQueryAtom(int [,bool=False])"""
    ...

def InNRingsGreaterQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where InNRings is equal to the target value.
    NOTE: the direction of comparison is reversed relative to the C++ API

    C++ signature :RDKit::QueryAtom* InNRingsGreaterQueryAtom(int [,bool=False])"""
    ...

def InNRingsLessQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where InNRings is less than the target value.
    NOTE: the direction of comparison is reversed relative to the C++ API

    C++ signature :RDKit::QueryAtom* InNRingsLessQueryAtom(int [,bool=False])"""
    ...

def IsAliphaticQueryAtom(self, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms when IsAliphatic is True.

    C++ signature :RDKit::QueryAtom* IsAliphaticQueryAtom([ bool=False])"""
    ...

def IsAromaticQueryAtom(self, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms when IsAromatic is True.

    C++ signature :RDKit::QueryAtom* IsAromaticQueryAtom([ bool=False])"""
    ...

def IsBridgeheadQueryAtom(self, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms when IsBridgehead is True.

    C++ signature :RDKit::QueryAtom* IsBridgeheadQueryAtom([ bool=False])"""
    ...

def IsInRingQueryAtom(self, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms when IsInRing is True.

    C++ signature :RDKit::QueryAtom* IsInRingQueryAtom([ bool=False])"""
    ...

def IsUnsaturatedQueryAtom(self, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms when IsUnsaturated is True.

    C++ signature :RDKit::QueryAtom* IsUnsaturatedQueryAtom([ bool=False])"""
    ...

def IsotopeEqualsQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where Isotope is equal to the target value.

    C++ signature :RDKit::QueryAtom* IsotopeEqualsQueryAtom(int [,bool=False])"""
    ...

def IsotopeGreaterQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where Isotope is equal to the target value.
    NOTE: the direction of comparison is reversed relative to the C++ API

    C++ signature :RDKit::QueryAtom* IsotopeGreaterQueryAtom(int [,bool=False])"""
    ...

def IsotopeLessQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where Isotope is less than the target value.
    NOTE: the direction of comparison is reversed relative to the C++ API

    C++ signature :RDKit::QueryAtom* IsotopeLessQueryAtom(int [,bool=False])"""
    ...

def MAtomQueryAtom(self, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms when MAtom is True.

    C++ signature :RDKit::QueryAtom* MAtomQueryAtom([ bool=False])"""
    ...

def MHAtomQueryAtom(self, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms when MHAtom is True.

    C++ signature :RDKit::QueryAtom* MHAtomQueryAtom([ bool=False])"""
    ...

def MassEqualsQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where Mass is equal to the target value.

    C++ signature :RDKit::QueryAtom* MassEqualsQueryAtom(int [,bool=False])"""
    ...

def MassGreaterQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where Mass is equal to the target value.
    NOTE: the direction of comparison is reversed relative to the C++ API

    C++ signature :RDKit::QueryAtom* MassGreaterQueryAtom(int [,bool=False])"""
    ...

def MassLessQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where Mass is less than the target value.
    NOTE: the direction of comparison is reversed relative to the C++ API

    C++ signature :RDKit::QueryAtom* MassLessQueryAtom(int [,bool=False])"""
    ...

def MinRingSizeEqualsQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where MinRingSize is equal to the target value.

    C++ signature :RDKit::QueryAtom* MinRingSizeEqualsQueryAtom(int [,bool=False])"""
    ...

def MinRingSizeGreaterQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where MinRingSize is equal to the target value.
    NOTE: the direction of comparison is reversed relative to the C++ API

    C++ signature :RDKit::QueryAtom* MinRingSizeGreaterQueryAtom(int [,bool=False])"""
    ...

def MinRingSizeLessQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where MinRingSize is less than the target value.
    NOTE: the direction of comparison is reversed relative to the C++ API

    C++ signature :RDKit::QueryAtom* MinRingSizeLessQueryAtom(int [,bool=False])"""
    ...

def MissingChiralTagQueryAtom(self, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms when MissingChiralTag is True.

    C++ signature :RDKit::QueryAtom* MissingChiralTagQueryAtom([ bool=False])"""
    ...

def NonHydrogenDegreeEqualsQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where NonHydrogenDegree is equal to the target value.

    C++ signature :RDKit::QueryAtom* NonHydrogenDegreeEqualsQueryAtom(int [,bool=False])
    """
    ...

def NonHydrogenDegreeGreaterQueryAtom(
    self, val: int, negate: bool = False
) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where NonHydrogenDegree is equal to the target value.
    NOTE: the direction of comparison is reversed relative to the C++ API

    C++ signature :RDKit::QueryAtom* NonHydrogenDegreeGreaterQueryAtom(int [,bool=False])
    """
    ...

def NonHydrogenDegreeLessQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where NonHydrogenDegree is less than the target value.
    NOTE: the direction of comparison is reversed relative to the C++ API

    C++ signature :RDKit::QueryAtom* NonHydrogenDegreeLessQueryAtom(int [,bool=False])
    """
    ...

def NumAliphaticHeteroatomNeighborsEqualsQueryAtom(
    self, val: int, negate: bool = False
) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where NumAliphaticHeteroatomNeighbors is equal to the target value.

    C++ signature :RDKit::QueryAtom* NumAliphaticHeteroatomNeighborsEqualsQueryAtom(int [,bool=False])
    """
    ...

def NumAliphaticHeteroatomNeighborsGreaterQueryAtom(
    self, val: int, negate: bool = False
) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where NumAliphaticHeteroatomNeighbors is equal to the target value.
    NOTE: the direction of comparison is reversed relative to the C++ API

    C++ signature :RDKit::QueryAtom* NumAliphaticHeteroatomNeighborsGreaterQueryAtom(int [,bool=False])
    """
    ...

def NumAliphaticHeteroatomNeighborsLessQueryAtom(
    self, val: int, negate: bool = False
) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where NumAliphaticHeteroatomNeighbors is less than the target value.
    NOTE: the direction of comparison is reversed relative to the C++ API

    C++ signature :RDKit::QueryAtom* NumAliphaticHeteroatomNeighborsLessQueryAtom(int [,bool=False])
    """
    ...

def NumHeteroatomNeighborsEqualsQueryAtom(
    self, val: int, negate: bool = False
) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where NumHeteroatomNeighbors is equal to the target value.

    C++ signature :RDKit::QueryAtom* NumHeteroatomNeighborsEqualsQueryAtom(int [,bool=False])
    """
    ...

def NumHeteroatomNeighborsGreaterQueryAtom(
    self, val: int, negate: bool = False
) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where NumHeteroatomNeighbors is equal to the target value.
    NOTE: the direction of comparison is reversed relative to the C++ API

    C++ signature :RDKit::QueryAtom* NumHeteroatomNeighborsGreaterQueryAtom(int [,bool=False])
    """
    ...

def NumHeteroatomNeighborsLessQueryAtom(
    self, val: int, negate: bool = False
) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where NumHeteroatomNeighbors is less than the target value.
    NOTE: the direction of comparison is reversed relative to the C++ API

    C++ signature :RDKit::QueryAtom* NumHeteroatomNeighborsLessQueryAtom(int [,bool=False])
    """
    ...

def NumRadicalElectronsEqualsQueryAtom(
    self, val: int, negate: bool = False
) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where NumRadicalElectrons is equal to the target value.

    C++ signature :RDKit::QueryAtom* NumRadicalElectronsEqualsQueryAtom(int [,bool=False])
    """
    ...

def NumRadicalElectronsGreaterQueryAtom(
    self, val: int, negate: bool = False
) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where NumRadicalElectrons is equal to the target value.
    NOTE: the direction of comparison is reversed relative to the C++ API

    C++ signature :RDKit::QueryAtom* NumRadicalElectronsGreaterQueryAtom(int [,bool=False])
    """
    ...

def NumRadicalElectronsLessQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where NumRadicalElectrons is less than the target value.
    NOTE: the direction of comparison is reversed relative to the C++ API

    C++ signature :RDKit::QueryAtom* NumRadicalElectronsLessQueryAtom(int [,bool=False])
    """
    ...

def QAtomQueryAtom(self, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms when QAtom is True.

    C++ signature :RDKit::QueryAtom* QAtomQueryAtom([ bool=False])"""
    ...

def QHAtomQueryAtom(self, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms when QHAtom is True.

    C++ signature :RDKit::QueryAtom* QHAtomQueryAtom([ bool=False])"""
    ...

def RingBondCountEqualsQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where RingBondCount is equal to the target value.

    C++ signature :RDKit::QueryAtom* RingBondCountEqualsQueryAtom(int [,bool=False])"""
    ...

def RingBondCountGreaterQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where RingBondCount is equal to the target value.
    NOTE: the direction of comparison is reversed relative to the C++ API

    C++ signature :RDKit::QueryAtom* RingBondCountGreaterQueryAtom(int [,bool=False])"""
    ...

def RingBondCountLessQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where RingBondCount is less than the target value.
    NOTE: the direction of comparison is reversed relative to the C++ API

    C++ signature :RDKit::QueryAtom* RingBondCountLessQueryAtom(int [,bool=False])"""
    ...

def TotalDegreeEqualsQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where TotalDegree is equal to the target value.

    C++ signature :RDKit::QueryAtom* TotalDegreeEqualsQueryAtom(int [,bool=False])"""
    ...

def TotalDegreeGreaterQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where TotalDegree is equal to the target value.
    NOTE: the direction of comparison is reversed relative to the C++ API

    C++ signature :RDKit::QueryAtom* TotalDegreeGreaterQueryAtom(int [,bool=False])"""
    ...

def TotalDegreeLessQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where TotalDegree is less than the target value.
    NOTE: the direction of comparison is reversed relative to the C++ API

    C++ signature :RDKit::QueryAtom* TotalDegreeLessQueryAtom(int [,bool=False])"""
    ...

def TotalValenceEqualsQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where TotalValence is equal to the target value.

    C++ signature :RDKit::QueryAtom* TotalValenceEqualsQueryAtom(int [,bool=False])"""
    ...

def TotalValenceGreaterQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where TotalValence is equal to the target value.
    NOTE: the direction of comparison is reversed relative to the C++ API

    C++ signature :RDKit::QueryAtom* TotalValenceGreaterQueryAtom(int [,bool=False])"""
    ...

def TotalValenceLessQueryAtom(self, val: int, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms where TotalValence is less than the target value.
    NOTE: the direction of comparison is reversed relative to the C++ API

    C++ signature :RDKit::QueryAtom* TotalValenceLessQueryAtom(int [,bool=False])"""
    ...

def XAtomQueryAtom(self, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms when XAtom is True.

    C++ signature :RDKit::QueryAtom* XAtomQueryAtom([ bool=False])"""
    ...

def XHAtomQueryAtom(self, negate: bool = False) -> QueryAtom:
    """
    Returns a QueryAtom that matches atoms when XHAtom is True.

    C++ signature :RDKit::QueryAtom* XHAtomQueryAtom([ bool=False])"""
    ...
