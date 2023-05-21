"""
Module containing functions to compute molecular descriptors
"""
from typing import Any, ClassVar, overload

import Boost.Python
import rdkit.rdBase
from rdkit import rdBase
from rdkit.Chem.rdchem import Atom, Mol
from rdkit.DataStructs.cDataStructs import (
    ExplicitBitVect,
    IntSparseIntVect,
    LongSparseIntVect,
    UIntSparseIntVect,
)
from rdkit.rdBase import (
    _vectd,
    _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE,
)

_BCUT2D_version: str
_CalcAUTOCORR2D_version: str
_CalcAUTOCORR3D_version: str
_CalcAsphericity_version: str
_CalcChi0n_version: str
_CalcChi0v_version: str
_CalcChi1n_version: str
_CalcChi1v_version: str
_CalcChi2n_version: str
_CalcChi2v_version: str
_CalcChi3n_version: str
_CalcChi3v_version: str
_CalcChi4n_version: str
_CalcChi4v_version: str
_CalcChiNn_version: str
_CalcChiNv_version: str
_CalcCoulombMat_version: str
_CalcCrippenDescriptors_version: str
_CalcEMMcharges_version: str
_CalcEccentricity_version: str
_CalcExactMolWt_version: str
_CalcFractionCSP3_version: str
_CalcGETAWAY_version: str
_CalcHallKierAlpha_version: str
_CalcInertialShapeFactor_version: str
_CalcKappa1_version: str
_CalcKappa2_version: str
_CalcKappa3_version: str
_CalcLabuteASA_version: str
_CalcMORSE_version: str
_CalcMolFormula_version: str
_CalcMolWt_version: str
_CalcNPR1_version: str
_CalcNPR2_version: str
_CalcNumAliphaticCarbocycles_version: str
_CalcNumAliphaticHeterocycles_version: str
_CalcNumAliphaticRings_version: str
_CalcNumAmideBonds_version: str
_CalcNumAromaticCarbocycles_version: str
_CalcNumAromaticHeterocycles_version: str
_CalcNumAromaticRings_version: str
_CalcNumAtomStereoCenters_version: str
_CalcNumAtoms_version: str
_CalcNumBridgeheadAtoms_version: str
_CalcNumHBA_version: str
_CalcNumHBD_version: str
_CalcNumHeavyAtoms_version: str
_CalcNumHeteroatoms_version: str
_CalcNumHeterocycles_version: str
_CalcNumLipinskiHBA_version: str
_CalcNumLipinskiHBD_version: str
_CalcNumRings_version: str
_CalcNumRotatableBonds_version: str
_CalcNumSaturatedCarbocycles_version: str
_CalcNumSaturatedHeterocycles_version: str
_CalcNumSaturatedRings_version: str
_CalcNumSpiroAtoms_version: str
_CalcNumUnspecifiedAtomStereoCenters_version: str
_CalcPBF_version: str
_CalcPMI1_version: str
_CalcPMI2_version: str
_CalcPMI3_version: str
_CalcPhi_version: str
_CalcRDF_version: str
_CalcRadiusOfGyration_version: str
_CalcSpherocityIndex_version: str
_CalcTPSA_version: str
_CalcWHIM_version: str
_ConnectivityInvariants_version: str
_FeatureInvariants_version: str
_GetAtomFeatures_version: str
_MorganFingerprint_version: str

class AtomPairsParameters(Boost.Python.instance):
    """
    C++ signature :void __init__(_object*)"""

    atomTypes: rdBase._vectj = ...
    codeSize: int = ...
    numAtomPairFingerprintBits: int = ...
    numBranchBits: int = ...
    numChiralBits: int = ...
    numPathBits: int = ...
    numPiBits: int = ...
    numTypeBits: int = ...
    version: str = ...
    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class NumRotatableBondsOptions(Boost.Python.enum):
    """
    Options for generating rotatable bonds
    NonStrict - standard loose definitions
    Strict - stricter definition excluding amides, esters, etc
    StrictLinkages - adds rotors between rotatable bonds
    Default - Current RDKit default"""

    Default: NumRotatableBondsOptions = ...
    NonStrict: NumRotatableBondsOptions = ...
    Strict: NumRotatableBondsOptions = ...
    StrictLinkages: NumRotatableBondsOptions = ...
    names: dict[str, NumRotatableBondsOptions] = ...
    values: dict[int, NumRotatableBondsOptions] = ...
    __slots__: ClassVar[tuple] = ...

class Properties(Boost.Python.instance):
    """
    Property computation and registry system.  To compute all registered properties:
    mol = Chem.MolFromSmiles(‘c1ccccc1’)
    properties = rdMolDescriptors.Properties()
    for name, value in zip(properties.GetPropertyNames(), properties.ComputeProperties(mol)):

    print(name, value)

    To compute a subset
    properties = rdMolDescriptors.Properties([‘exactmw’, ‘lipinskiHBA’])
    for name, value in zip(properties.GetPropertyNames(), properties.ComputeProperties(mol)):

    print(name, value)

    C++ signature :void __init__(_object*)

    __init__( (object)arg1, (_vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE)arg2) -> None :

    C++ signature :void __init__(_object*,std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >)
    """

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def AnnotateProperties(self, arg1: Properties, mol: Mol) -> None:
        """
        Annotate the molecule with the computed properties.  These properties will be available as SDData or from mol.GetProp(prop)

        C++ signature :void AnnotateProperties(RDKit::Descriptors::Properties {lvalue},RDKit::ROMol {lvalue})
        """
        ...
    def ComputeProperties(
        self, arg1: Properties, mol: Mol, annotateMol: bool = False
    ) -> _vectd:
        """
        Return a list of computed properties, if annotateMol==True, annotate the molecule with the computed properties.

        C++ signature :std::vector<double, std::allocator<double> > ComputeProperties(RDKit::Descriptors::Properties {lvalue},RDKit::ROMol [,bool=False])
        """
        ...
    def GetAvailableProperties(
        self,
    ) -> _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE:
        """
        Return all available property names that can be computed

        C++ signature :std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > GetAvailableProperties()
        """
        ...
    def GetProperty(self, propName: str) -> PropertyFunctor:
        """
        Return the named property if it exists

        C++ signature :boost::shared_ptr<RDKit::Descriptors::PropertyFunctor> GetProperty(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    def GetPropertyNames(
        self, arg1: Properties
    ) -> _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE:
        """
        Return the property names computed by this instance

        C++ signature :std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > GetPropertyNames(RDKit::Descriptors::Properties {lvalue})
        """
        ...
    def RegisterProperty(self, propertyFunctor: PropertyFunctor) -> int:
        """
        Register a new property object (not thread safe)

        C++ signature :int RegisterProperty(RDKit::Descriptors::PropertyFunctor*)"""
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class PropertyFunctor(Boost.Python.instance):
    """
    Property computation class stored in the property registry.
    See rdkit.Chem.rdMolDescriptor.Properties.GetProperty and
    rdkit.Chem.Descriptor.Properties.PropertyFunctor for creating new ones
    Raises an exception
    This class cannot be instantiated from Python"""

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def GetName(self, arg1: PropertyFunctor) -> str:
        """
        Return the name of the property to calculate

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetName(RDKit::Descriptors::PropertyFunctor {lvalue})
        """
        ...
    def GetVersion(self, arg1: PropertyFunctor) -> str:
        """
        Return the version of the calculated property

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetVersion(RDKit::Descriptors::PropertyFunctor {lvalue})
        """
        ...
    @classmethod
    def __call__(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class PropertyRangeQuery(Boost.Python.instance):
    """
    Property Range Query for a molecule.  Match(mol) -> true if in range
    Raises an exception
    This class cannot be instantiated from Python"""

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def Match(self, arg1: PropertyRangeQuery, arg2: Mol) -> bool:
        """
        C++ signature :bool Match(Queries::RangeQuery<double, RDKit::ROMol const&, true> {lvalue},RDKit::ROMol)
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class PythonPropertyFunctor(PropertyFunctor):
    """
    C++ signature :void __init__(_object*,_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """

    ...
    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __call__(cls, anonymousnamespace, RDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

@overload
def BCUT2D(self, mol: Mol) -> list:
    """

    Returns a 2D BCUT (eigen value hi, eigenvalue low) given the molecule and the specified atom propsatom_props must be a list or tuple of floats equal in size to the number of atoms in mol

    C++ signature :std::pair<double, double> BCUT2D(RDKit::ROMol,boost::python::list)

    BCUT2D( (Mol)mol, (tuple)atom_props) -> tuple :
    Returns a 2D BCUT (eigen value hi, eigenvalue low) given the molecule and the specified atom propsatom_props must be a list or tuple of floats equal in size to the number of atoms in mol

    C++ signature :std::pair<double, double> BCUT2D(RDKit::ROMol,boost::python::tuple)

    BCUT2D( (Mol)mol, (str)atom_propname) -> tuple :Returns a 2D BCUT (eigen value high, eigen value low) given the molecule and the specified atom prop name
    atom_propname must exist on each atom and be convertible to a float

    C++ signature :std::pair<double, double> BCUT2D(RDKit::ROMol,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
    ...

@overload
def BCUT2D(self, mol: Mol, atom_props: list) -> tuple:
    """

    Returns a 2D BCUT (eigen value hi, eigenvalue low) given the molecule and the specified atom propsatom_props must be a list or tuple of floats equal in size to the number of atoms in mol

    C++ signature :std::pair<double, double> BCUT2D(RDKit::ROMol,boost::python::tuple)

    BCUT2D( (Mol)mol, (str)atom_propname) -> tuple :Returns a 2D BCUT (eigen value high, eigen value low) given the molecule and the specified atom prop name
    atom_propname must exist on each atom and be convertible to a float

    C++ signature :std::pair<double, double> BCUT2D(RDKit::ROMol,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
    ...

@overload
def BCUT2D(self, mol: Mol, atom_props: tuple) -> tuple:
    """
    Returns a 2D BCUT (eigen value high, eigen value low) given the molecule and the specified atom prop name
    atom_propname must exist on each atom and be convertible to a float

    C++ signature :std::pair<double, double> BCUT2D(RDKit::ROMol,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
    ...

@overload
def BCUT2D(self, mol: Mol, atom_propname: str) -> tuple: ...
def CalcAUTOCORR2D(self, mol: Mol, CustomAtomProperty: str = "") -> list:
    """
    Returns 2D Autocorrelation descriptors vector

    C++ signature :boost::python::list CalcAUTOCORR2D(RDKit::ROMol [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’’])
    """
    ...

def CalcAUTOCORR3D(
    self, mol: Mol, confId: int = -1, CustomAtomProperty: str = ""
) -> list:
    """
    Returns 3D Autocorrelation descriptors vector

    C++ signature :boost::python::list CalcAUTOCORR3D(RDKit::ROMol [,int=-1 [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’’]])
    """
    ...

def CalcAsphericity(
    self, mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True
) -> float:
    """
    C++ signature :double CalcAsphericity(RDKit::ROMol [,int=-1 [,bool=True [,bool=True]]])
    """
    ...

def CalcChi0n(self, mol: Mol, force: bool = False) -> float:
    """
    C++ signature :double CalcChi0n(RDKit::ROMol [,bool=False])"""
    ...

def CalcChi0v(self, mol: Mol, force: bool = False) -> float:
    """
    C++ signature :double CalcChi0v(RDKit::ROMol [,bool=False])"""
    ...

def CalcChi1n(self, mol: Mol, force: bool = False) -> float:
    """
    C++ signature :double CalcChi1n(RDKit::ROMol [,bool=False])"""
    ...

def CalcChi1v(self, mol: Mol, force: bool = False) -> float:
    """
    C++ signature :double CalcChi1v(RDKit::ROMol [,bool=False])"""
    ...

def CalcChi2n(self, mol: Mol, force: bool = False) -> float:
    """
    C++ signature :double CalcChi2n(RDKit::ROMol [,bool=False])"""
    ...

def CalcChi2v(self, mol: Mol, force: bool = False) -> float:
    """
    C++ signature :double CalcChi2v(RDKit::ROMol [,bool=False])"""
    ...

def CalcChi3n(self, mol: Mol, force: bool = False) -> float:
    """
    C++ signature :double CalcChi3n(RDKit::ROMol [,bool=False])"""
    ...

def CalcChi3v(self, mol: Mol, force: bool = False) -> float:
    """
    C++ signature :double CalcChi3v(RDKit::ROMol [,bool=False])"""
    ...

def CalcChi4n(self, mol: Mol, force: bool = False) -> float:
    """
    C++ signature :double CalcChi4n(RDKit::ROMol [,bool=False])"""
    ...

def CalcChi4v(self, mol: Mol, force: bool = False) -> float:
    """
    C++ signature :double CalcChi4v(RDKit::ROMol [,bool=False])"""
    ...

def CalcChiNn(self, mol: Mol, n: int, force: bool = False) -> float:
    """
    C++ signature :double CalcChiNn(RDKit::ROMol,unsigned int [,bool=False])"""
    ...

def CalcChiNv(self, mol: Mol, n: int, force: bool = False) -> float:
    """
    C++ signature :double CalcChiNv(RDKit::ROMol,unsigned int [,bool=False])"""
    ...

def CalcCoulombMat(self, mol: Mol, confId: int = -1) -> tuple:
    """
    Returns severals Coulomb randomized matrices

    C++ signature :boost::python::tuple CalcCoulombMat(RDKit::ROMol [,int=-1])"""
    ...

def CalcCrippenDescriptors(
    self, mol: Mol, includeHs: bool = True, force: bool = False
) -> tuple:
    """
    returns a 2-tuple with the Wildman-Crippen logp,mr values

    C++ signature :boost::python::tuple CalcCrippenDescriptors(RDKit::ROMol [,bool=True [,bool=False]])
    """
    ...

def CalcEEMcharges(self, mol: Mol, confId: int = -1) -> list:
    """
    Returns EEM atomic partial charges

    C++ signature :boost::python::list CalcEEMcharges(RDKit::ROMol {lvalue} [,int=-1])
    """
    ...

def CalcEccentricity(
    self, mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True
) -> float:
    """
    C++ signature :double CalcEccentricity(RDKit::ROMol [,int=-1 [,bool=True [,bool=True]]])
    """
    ...

def CalcExactMolWt(self, mol: Mol, onlyHeavy: bool = False) -> float:
    """
    returns the molecule’s exact molecular weight

    C++ signature :double CalcExactMolWt(RDKit::ROMol [,bool=False])"""
    ...

def CalcFractionCSP3(self, mol: Mol) -> float:
    """
    returns the fraction of C atoms that are SP3 hybridized

    C++ signature :double CalcFractionCSP3(RDKit::ROMol)"""
    ...

def CalcGETAWAY(
    self, mol: Mol, confId: int = -1, precision: float = 2, CustomAtomProperty: str = ""
) -> list:
    """
    Returns the GETAWAY descriptors vector

    C++ signature :boost::python::list CalcGETAWAY(RDKit::ROMol [,int=-1 [,double=2 [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’’]]])
    """
    ...

def CalcHallKierAlpha(
    self, mol: Mol, atomContribs: AtomPairsParameters = None
) -> float:
    """
    C++ signature :double CalcHallKierAlpha(RDKit::ROMol [,boost::python::api::object=None])
    """
    ...

def CalcInertialShapeFactor(
    self, mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True
) -> float:
    """
    C++ signature :double CalcInertialShapeFactor(RDKit::ROMol [,int=-1 [,bool=True [,bool=True]]])
    """
    ...

def CalcKappa1(self, mol: Mol) -> float:
    """
    C++ signature :double CalcKappa1(RDKit::ROMol)"""
    ...

def CalcKappa2(self, mol: Mol) -> float:
    """
    C++ signature :double CalcKappa2(RDKit::ROMol)"""
    ...

def CalcKappa3(self, mol: Mol) -> float:
    """
    C++ signature :double CalcKappa3(RDKit::ROMol)"""
    ...

def CalcLabuteASA(self, mol: Mol, includeHs: bool = True, force: bool = False) -> float:
    """
    returns the Labute ASA value for a molecule

    C++ signature :double CalcLabuteASA(RDKit::ROMol [,bool=True [,bool=False]])"""
    ...

def CalcMORSE(self, mol: Mol, confId: int = -1, CustomAtomProperty: str = "") -> list:
    """
    Returns Molecule Representation of Structures based on Electron diffraction descriptors

    C++ signature :boost::python::list CalcMORSE(RDKit::ROMol [,int=-1 [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’’]])
    """
    ...

def CalcMolFormula(
    self, mol: Mol, separateIsotopes: bool = False, abbreviateHIsotopes: bool = True
) -> str:
    """
    returns the molecule’s formula

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > CalcMolFormula(RDKit::ROMol [,bool=False [,bool=True]])
    """
    ...

def CalcNPR1(
    self, mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True
) -> float:
    """
    C++ signature :double CalcNPR1(RDKit::ROMol [,int=-1 [,bool=True [,bool=True]]])"""
    ...

def CalcNPR2(
    self, mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True
) -> float:
    """
    C++ signature :double CalcNPR2(RDKit::ROMol [,int=-1 [,bool=True [,bool=True]]])"""
    ...

def CalcNumAliphaticCarbocycles(self, mol: Mol) -> int:
    """
    returns the number of aliphatic (containing at least one non-aromatic bond) carbocycles for a molecule

    C++ signature :unsigned int CalcNumAliphaticCarbocycles(RDKit::ROMol)"""
    ...

def CalcNumAliphaticHeterocycles(self, mol: Mol) -> int:
    """
    returns the number of aliphatic (containing at least one non-aromatic bond) heterocycles for a molecule

    C++ signature :unsigned int CalcNumAliphaticHeterocycles(RDKit::ROMol)"""
    ...

def CalcNumAliphaticRings(self, mol: Mol) -> int:
    """
    returns the number of aliphatic (containing at least one non-aromatic bond) rings for a molecule

    C++ signature :unsigned int CalcNumAliphaticRings(RDKit::ROMol)"""
    ...

def CalcNumAmideBonds(self, mol: Mol) -> int:
    """
    returns the number of amide bonds in a molecule

    C++ signature :unsigned int CalcNumAmideBonds(RDKit::ROMol)"""
    ...

def CalcNumAromaticCarbocycles(self, mol: Mol) -> int:
    """
    returns the number of aromatic carbocycles for a molecule

    C++ signature :unsigned int CalcNumAromaticCarbocycles(RDKit::ROMol)"""
    ...

def CalcNumAromaticHeterocycles(self, mol: Mol) -> int:
    """
    returns the number of aromatic heterocycles for a molecule

    C++ signature :unsigned int CalcNumAromaticHeterocycles(RDKit::ROMol)"""
    ...

def CalcNumAromaticRings(self, mol: Mol) -> int:
    """
    returns the number of aromatic rings for a molecule

    C++ signature :unsigned int CalcNumAromaticRings(RDKit::ROMol)"""
    ...

def CalcNumAtomStereoCenters(self, mol: Mol) -> int:
    """
    Returns the total number of atomic stereocenters (specified and unspecified)

    C++ signature :unsigned int CalcNumAtomStereoCenters(RDKit::ROMol)"""
    ...

def CalcNumAtoms(self, mol: Mol) -> int:
    """
    returns the total number of atoms for a molecule

    C++ signature :unsigned int CalcNumAtoms(RDKit::ROMol)"""
    ...

def CalcNumBridgeheadAtoms(self, mol: Mol, atoms: AtomPairsParameters = None) -> int:
    """
    Returns the number of bridgehead atoms (atoms shared between rings that share at least two bonds)

    C++ signature :unsigned int CalcNumBridgeheadAtoms(RDKit::ROMol [,boost::python::api::object=None])
    """
    ...

def CalcNumHBA(self, mol: Mol) -> int:
    """
    returns the number of H-bond acceptors for a molecule

    C++ signature :unsigned int CalcNumHBA(RDKit::ROMol)"""
    ...

def CalcNumHBD(self, mol: Mol) -> int:
    """
    returns the number of H-bond donors for a molecule

    C++ signature :unsigned int CalcNumHBD(RDKit::ROMol)"""
    ...

def CalcNumHeavyAtoms(self, mol: Mol) -> int:
    """
    returns the number of heavy atoms for a molecule

    C++ signature :unsigned int CalcNumHeavyAtoms(RDKit::ROMol)"""
    ...

def CalcNumHeteroatoms(self, mol: Mol) -> int:
    """
    returns the number of heteroatoms for a molecule

    C++ signature :unsigned int CalcNumHeteroatoms(RDKit::ROMol)"""
    ...

def CalcNumHeterocycles(self, mol: Mol) -> int:
    """
    returns the number of heterocycles for a molecule

    C++ signature :unsigned int CalcNumHeterocycles(RDKit::ROMol)"""
    ...

def CalcNumLipinskiHBA(self, mol: Mol) -> int:
    """
    returns the number of Lipinski H-bond acceptors for a molecule

    C++ signature :unsigned int CalcNumLipinskiHBA(RDKit::ROMol)"""
    ...

def CalcNumLipinskiHBD(self, mol: Mol) -> int:
    """
    returns the number of Lipinski H-bond donors for a molecule

    C++ signature :unsigned int CalcNumLipinskiHBD(RDKit::ROMol)"""
    ...

def CalcNumRings(self, mol: Mol) -> int:
    """
    returns the number of rings for a molecule

    C++ signature :unsigned int CalcNumRings(RDKit::ROMol)"""
    ...

@overload
def CalcNumRotatableBonds(self, mol: Mol, strict: bool) -> int:
    """

    returns the number of rotatable bonds for a molecule.strict = NumRotatableBondsOptions.NonStrict - Simple rotatable bond definition.
    strict = NumRotatableBondsOptions.Strict - (default) does not count things like

    amide or ester bonds

    strict = NumRotatableBondsOptions.StrictLinkages - handles linkages between ringsystems.
    - Single bonds between aliphatic ring Cs are always rotatable. This

    means that the central bond in CC1CCCC(C)C1-C1C(C)CCCC1C is now
    considered rotatable; it was not before

    Heteroatoms in the linked rings no longer affect whether or not
    the linking bond is rotatable

    the linking bond in systems like Cc1cccc(C)c1-c1c(C)cccc1 is nowconsidered non-rotatable

    C++ signature :unsigned int CalcNumRotatableBonds(RDKit::ROMol [,RDKit::Descriptors::NumRotatableBondsOptions=rdkit.Chem.rdMolDescriptors.NumRotatableBondsOptions.Default])
    """
    ...

@overload
def CalcNumRotatableBonds(
    self, mol: Mol, strict: NumRotatableBondsOptions = NumRotatableBondsOptions.Default
) -> int: ...
def CalcNumSaturatedCarbocycles(self, mol: Mol) -> int:
    """
    returns the number of saturated carbocycles for a molecule

    C++ signature :unsigned int CalcNumSaturatedCarbocycles(RDKit::ROMol)"""
    ...

def CalcNumSaturatedHeterocycles(self, mol: Mol) -> int:
    """
    returns the number of saturated heterocycles for a molecule

    C++ signature :unsigned int CalcNumSaturatedHeterocycles(RDKit::ROMol)"""
    ...

def CalcNumSaturatedRings(self, mol: Mol) -> int:
    """
    returns the number of saturated rings for a molecule

    C++ signature :unsigned int CalcNumSaturatedRings(RDKit::ROMol)"""
    ...

def CalcNumSpiroAtoms(self, mol: Mol, atoms: AtomPairsParameters = None) -> int:
    """
    Returns the number of spiro atoms (atoms shared between rings that share exactly one atom)

    C++ signature :unsigned int CalcNumSpiroAtoms(RDKit::ROMol [,boost::python::api::object=None])
    """
    ...

def CalcNumUnspecifiedAtomStereoCenters(self, mol: Mol) -> int:
    """
    Returns the number of unspecified atomic stereocenters

    C++ signature :unsigned int CalcNumUnspecifiedAtomStereoCenters(RDKit::ROMol)"""
    ...

def CalcOxidationNumbers(self, mol: Mol) -> None:
    """
    Adds the oxidation number/state to the atoms of a molecule as property OxidationNumber on each atom.  Use Pauling electronegativities.  This is experimental code, still under development.

    C++ signature :void CalcOxidationNumbers(RDKit::ROMol)"""
    ...

def CalcPBF(self, mol: Mol, confId: int = -1) -> float:
    """
    Returns the PBF (plane of best fit) descriptor (https://doi.org/10.1021/ci300293f)

    C++ signature :double CalcPBF(RDKit::ROMol [,int=-1])"""
    ...

def CalcPMI1(
    self, mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True
) -> float:
    """
    C++ signature :double CalcPMI1(RDKit::ROMol [,int=-1 [,bool=True [,bool=True]]])"""
    ...

def CalcPMI2(
    self, mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True
) -> float:
    """
    C++ signature :double CalcPMI2(RDKit::ROMol [,int=-1 [,bool=True [,bool=True]]])"""
    ...

def CalcPMI3(
    self, mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True
) -> float:
    """
    C++ signature :double CalcPMI3(RDKit::ROMol [,int=-1 [,bool=True [,bool=True]]])"""
    ...

def CalcPhi(self, mol: Mol) -> float:
    """
    C++ signature :double CalcPhi(RDKit::ROMol)"""
    ...

def CalcRDF(self, mol: Mol, confId: int = -1, CustomAtomProperty: str = "") -> list:
    """
    Returns radial distribution fonction descriptors (RDF)

    C++ signature :boost::python::list CalcRDF(RDKit::ROMol [,int=-1 [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’’]])
    """
    ...

def CalcRadiusOfGyration(
    self, mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True
) -> float:
    """
    C++ signature :double CalcRadiusOfGyration(RDKit::ROMol [,int=-1 [,bool=True [,bool=True]]])
    """
    ...

def CalcSpherocityIndex(self, mol: Mol, confId: int = -1, force: bool = True) -> float:
    """
    C++ signature :double CalcSpherocityIndex(RDKit::ROMol [,int=-1 [,bool=True]])"""
    ...

def CalcTPSA(self, mol: Mol, force: bool = False, includeSandP: bool = False) -> float:
    """
    returns the TPSA value for a molecule

    C++ signature :double CalcTPSA(RDKit::ROMol [,bool=False [,bool=False]])"""
    ...

def CalcWHIM(
    self,
    mol: Mol,
    confId: int = -1,
    thresh: float = 0.001,
    CustomAtomProperty: str = "",
) -> list:
    """
    Returns the WHIM descriptors vector

    C++ signature :boost::python::list CalcWHIM(RDKit::ROMol [,int=-1 [,double=0.001 [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’’]]])
    """
    ...

def CustomProp_VSA_(
    self, mol: Mol, customPropName: str, bins: AtomPairsParameters, force: bool = False
) -> list:
    """
    C++ signature :boost::python::list CustomProp_VSA_(RDKit::ROMol,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,boost::python::api::object [,bool=False])
    """
    ...

def GetAtomFeatures(self, mol: Mol, atomid: int, addchiral: bool = False) -> list:
    """
    Returns the Atom Features vector

    C++ signature :boost::python::list GetAtomFeatures(RDKit::ROMol,int [,bool=False])
    """
    ...

def GetAtomPairAtomCode(
    self, atom: Atom, branchSubtract: int = 0, includeChirality: bool = False
) -> int:
    """
    Returns the atom code (hash) for an atom

    C++ signature :unsigned int GetAtomPairAtomCode(RDKit::Atom const* [,unsigned int=0 [,bool=False]])
    """
    ...

def GetAtomPairCode(
    self, atom1Code: int, atom2Code: int, distance: int, includeChirality: bool = False
) -> int:
    """
    Returns the atom-pair code (hash) for a pair of atoms separated by a certain number of bonds

    C++ signature :unsigned int GetAtomPairCode(unsigned int,unsigned int,unsigned int [,bool=False])
    """
    ...

def GetAtomPairFingerprint(
    self,
    mol: Mol,
    minLength: int = 1,
    maxLength: int = 30,
    fromAtoms: AtomPairsParameters = 0,
    ignoreAtoms: AtomPairsParameters = 0,
    atomInvariants: AtomPairsParameters = 0,
    includeChirality: bool = False,
    use2D: bool = True,
    confId: int = -1,
) -> IntSparseIntVect:
    """
    Returns the atom-pair fingerprint for a molecule as an IntSparseIntVect

    C++ signature :RDKit::SparseIntVect<int>* GetAtomPairFingerprint(RDKit::ROMol [,unsigned int=1 [,unsigned int=30 [,boost::python::api::object=0 [,boost::python::api::object=0 [,boost::python::api::object=0 [,bool=False [,bool=True [,int=-1]]]]]]]])
    """
    ...

def GetConnectivityInvariants(
    self, mol: Mol, includeRingMembership: bool = True
) -> list:
    """
    Returns connectivity invariants (ECFP-like) for a molecule.

    C++ signature :boost::python::list GetConnectivityInvariants(RDKit::ROMol [,bool=True])
    """
    ...

def GetFeatureInvariants(self, mol: Mol) -> list:
    """
    Returns feature invariants (FCFP-like) for a molecule.

    C++ signature :boost::python::list GetFeatureInvariants(RDKit::ROMol)"""
    ...

def GetHashedAtomPairFingerprint(
    self,
    mol: Mol,
    nBits: int = 2048,
    minLength: int = 1,
    maxLength: int = 30,
    fromAtoms: AtomPairsParameters = 0,
    ignoreAtoms: AtomPairsParameters = 0,
    atomInvariants: AtomPairsParameters = 0,
    includeChirality: bool = False,
    use2D: bool = True,
    confId: int = -1,
) -> IntSparseIntVect:
    """
    Returns the hashed atom-pair fingerprint for a molecule as an IntSparseIntVect

    C++ signature :RDKit::SparseIntVect<int>* GetHashedAtomPairFingerprint(RDKit::ROMol [,unsigned int=2048 [,unsigned int=1 [,unsigned int=30 [,boost::python::api::object=0 [,boost::python::api::object=0 [,boost::python::api::object=0 [,bool=False [,bool=True [,int=-1]]]]]]]]])
    """
    ...

def GetHashedAtomPairFingerprintAsBitVect(
    self,
    mol: Mol,
    nBits: int = 2048,
    minLength: int = 1,
    maxLength: int = 30,
    fromAtoms: AtomPairsParameters = 0,
    ignoreAtoms: AtomPairsParameters = 0,
    atomInvariants: AtomPairsParameters = 0,
    nBitsPerEntry: int = 4,
    includeChirality: bool = False,
    use2D: bool = True,
    confId: int = -1,
) -> ExplicitBitVect:
    """
    Returns the atom-pair fingerprint for a molecule as an ExplicitBitVect

    C++ signature :ExplicitBitVect* GetHashedAtomPairFingerprintAsBitVect(RDKit::ROMol [,unsigned int=2048 [,unsigned int=1 [,unsigned int=30 [,boost::python::api::object=0 [,boost::python::api::object=0 [,boost::python::api::object=0 [,unsigned int=4 [,bool=False [,bool=True [,int=-1]]]]]]]]]])
    """
    ...

def GetHashedMorganFingerprint(
    self,
    mol: Mol,
    radius: int,
    nBits: int = 2048,
    invariants: AtomPairsParameters = [],
    fromAtoms: AtomPairsParameters = [],
    useChirality: bool = False,
    useBondTypes: bool = True,
    useFeatures: bool = False,
    bitInfo: AtomPairsParameters = None,
    includeRedundantEnvironments: bool = False,
) -> UIntSparseIntVect:
    """
    Returns a hashed Morgan fingerprint for a molecule

    C++ signature :RDKit::SparseIntVect<unsigned int>* GetHashedMorganFingerprint(RDKit::ROMol,unsigned int [,unsigned int=2048 [,boost::python::api::object=[] [,boost::python::api::object=[] [,bool=False [,bool=True [,bool=False [,boost::python::api::object=None [,bool=False]]]]]]]])
    """
    ...

def GetHashedTopologicalTorsionFingerprint(
    self,
    mol: Mol,
    nBits: int = 2048,
    targetSize: int = 4,
    fromAtoms: AtomPairsParameters = 0,
    ignoreAtoms: AtomPairsParameters = 0,
    atomInvariants: AtomPairsParameters = 0,
    includeChirality: bool = False,
) -> LongSparseIntVect:
    """
    Returns the hashed topological-torsion fingerprint for a molecule as a LongIntSparseIntVect

    C++ signature :RDKit::SparseIntVect<long>* GetHashedTopologicalTorsionFingerprint(RDKit::ROMol [,unsigned int=2048 [,unsigned int=4 [,boost::python::api::object=0 [,boost::python::api::object=0 [,boost::python::api::object=0 [,bool=False]]]]]])
    """
    ...

def GetHashedTopologicalTorsionFingerprintAsBitVect(
    self,
    mol: Mol,
    nBits: int = 2048,
    targetSize: int = 4,
    fromAtoms: AtomPairsParameters = 0,
    ignoreAtoms: AtomPairsParameters = 0,
    atomInvariants: AtomPairsParameters = 0,
    nBitsPerEntry: int = 4,
    includeChirality: bool = False,
) -> ExplicitBitVect:
    """
    Returns the topological-torsion fingerprint for a molecule as an ExplicitBitVect

    C++ signature :ExplicitBitVect* GetHashedTopologicalTorsionFingerprintAsBitVect(RDKit::ROMol [,unsigned int=2048 [,unsigned int=4 [,boost::python::api::object=0 [,boost::python::api::object=0 [,boost::python::api::object=0 [,unsigned int=4 [,bool=False]]]]]]])
    """
    ...

def GetMACCSKeysFingerprint(self, mol: Mol) -> ExplicitBitVect:
    """
    Returns the MACCS keys for a molecule as an ExplicitBitVect

    C++ signature :ExplicitBitVect* GetMACCSKeysFingerprint(RDKit::ROMol)"""
    ...

def GetMorganFingerprint(
    self,
    mol: Mol,
    radius: int,
    invariants: AtomPairsParameters = [],
    fromAtoms: AtomPairsParameters = [],
    useChirality: bool = False,
    useBondTypes: bool = True,
    useFeatures: bool = False,
    useCounts: bool = True,
    bitInfo: AtomPairsParameters = None,
    includeRedundantEnvironments: bool = False,
) -> UIntSparseIntVect:
    """
    Returns a Morgan fingerprint for a molecule

    C++ signature :RDKit::SparseIntVect<unsigned int>* GetMorganFingerprint(RDKit::ROMol,unsigned int [,boost::python::api::object=[] [,boost::python::api::object=[] [,bool=False [,bool=True [,bool=False [,bool=True [,boost::python::api::object=None [,bool=False]]]]]]]])
    """
    ...

def GetMorganFingerprintAsBitVect(
    self,
    mol: Mol,
    radius: int,
    nBits: int = 2048,
    invariants: AtomPairsParameters = [],
    fromAtoms: AtomPairsParameters = [],
    useChirality: bool = False,
    useBondTypes: bool = True,
    useFeatures: bool = False,
    bitInfo: AtomPairsParameters = None,
    includeRedundantEnvironments: bool = False,
) -> ExplicitBitVect:
    """
    Returns a Morgan fingerprint for a molecule as a bit vector

    C++ signature :ExplicitBitVect* GetMorganFingerprintAsBitVect(RDKit::ROMol,unsigned int [,unsigned int=2048 [,boost::python::api::object=[] [,boost::python::api::object=[] [,bool=False [,bool=True [,bool=False [,boost::python::api::object=None [,bool=False]]]]]]]])
    """
    ...

def GetTopologicalTorsionFingerprint(
    self,
    mol: Mol,
    targetSize: int = 4,
    fromAtoms: AtomPairsParameters = 0,
    ignoreAtoms: AtomPairsParameters = 0,
    atomInvariants: AtomPairsParameters = 0,
    includeChirality: bool = False,
) -> LongSparseIntVect:
    """
    Returns the topological-torsion fingerprint for a molecule as a LongIntSparseIntVect

    C++ signature :RDKit::SparseIntVect<long>* GetTopologicalTorsionFingerprint(RDKit::ROMol [,unsigned int=4 [,boost::python::api::object=0 [,boost::python::api::object=0 [,boost::python::api::object=0 [,bool=False]]]]])
    """
    ...

def GetUSR(self, mol: Mol, confId: int = -1) -> list:
    """
    Returns a USR descriptor for one conformer of a molecule

    C++ signature :boost::python::list GetUSR(RDKit::ROMol [,int=-1])"""
    ...

def GetUSRCAT(
    self, mol: Mol, atomSelections: AtomPairsParameters = None, confId: int = -1
) -> list:
    """
    Returns a USRCAT descriptor for one conformer of a molecule

    C++ signature :boost::python::list GetUSRCAT(RDKit::ROMol [,boost::python::api::object=None [,int=-1]])
    """
    ...

def GetUSRDistributions(
    self, coords: AtomPairsParameters, points: AtomPairsParameters = None
) -> list:
    """
    Returns the four USR distance distributions for a set of coordinates

    C++ signature :boost::python::list GetUSRDistributions(boost::python::api::object [,boost::python::api::object=None])
    """
    ...

def GetUSRDistributionsFromPoints(
    self, coords: AtomPairsParameters, points: AtomPairsParameters
) -> list:
    """
    Returns the USR distance distributions for a set of coordinates and points

    C++ signature :boost::python::list GetUSRDistributionsFromPoints(boost::python::api::object,boost::python::api::object)
    """
    ...

def GetUSRFromDistributions(self, distances: AtomPairsParameters) -> list:
    """
    Returns the USR descriptor from a set of distance distributions

    C++ signature :boost::python::list GetUSRFromDistributions(boost::python::api::object)
    """
    ...

def GetUSRScore(
    self,
    descriptor1: AtomPairsParameters,
    descriptor2: AtomPairsParameters,
    weights: AtomPairsParameters = [],
) -> float:
    """
    Returns the USR score for two USR or USRCAT descriptors

    C++ signature :double GetUSRScore(boost::python::api::object,boost::python::api::object [,boost::python::api::object=[]])
    """
    ...

def MQNs_(self, mol: Mol, force: bool = False) -> list:
    """
    C++ signature :boost::python::list MQNs_(RDKit::ROMol [,bool=False])"""
    ...

def MakePropertyRangeQuery(
    self, name: str, min: float, max: float
) -> PropertyRangeQuery:
    """
    Generates a Range property for the specified property, between min and max
    query = MakePropertyRangeQuery(‘exactmw’, 0, 500)
    query.Match( mol )

    C++ signature :Queries::RangeQuery<double, RDKit::ROMol const&, true>* MakePropertyRangeQuery(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,double,double)
    """
    ...

def PEOE_VSA_(
    self, mol: Mol, bins: AtomPairsParameters = [], force: bool = False
) -> list:
    """
    C++ signature :boost::python::list PEOE_VSA_(RDKit::ROMol [,boost::python::api::object=[] [,bool=False]])
    """
    ...

def SMR_VSA_(
    self, mol: Mol, bins: AtomPairsParameters = [], force: bool = False
) -> list:
    """
    C++ signature :boost::python::list SMR_VSA_(RDKit::ROMol [,boost::python::api::object=[] [,bool=False]])
    """
    ...

def SlogP_VSA_(
    self, mol: Mol, bins: AtomPairsParameters = [], force: bool = False
) -> list:
    """
    C++ signature :boost::python::list SlogP_VSA_(RDKit::ROMol [,boost::python::api::object=[] [,bool=False]])
    """
    ...

def _CalcCrippenContribs(RDKit) -> Any: ...
def _CalcLabuteASAContribs(RDKit) -> Any: ...
def _CalcMolWt(RDKit) -> Any: ...
def _CalcTPSAContribs(RDKit) -> Any: ...
