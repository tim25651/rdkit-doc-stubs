"""
rdkit.Chem.rdFreeSASA module¶
Module containing rdFreeSASA classes and functions.
"""
from typing import Any, ClassVar

import Boost.Python
from rdkit.Chem.rdchem import Atom, Mol, QueryAtom
from rdkit.Chem.rdMolDescriptors import AtomPairsParameters

APolar: SASAClass
LeeRichards: SASAAlgorithm
NACCESS: SASAClassifier
OONS: SASAClassifier
Polar: SASAClass
Protor: SASAClassifier
ShrakeRupley: SASAAlgorithm
Unclassified: SASAClass

class SASAAlgorithm(Boost.Python.enum):
    LeeRichards: SASAAlgorithm = ...
    ShrakeRupley: SASAAlgorithm = ...
    names: dict[str, SASAAlgorithm] = ...
    values: dict[int, SASAAlgorithm] = ...
    __slots__: ClassVar[tuple] = ...

class SASAClass(Boost.Python.enum):
    APolar: SASAClass = ...
    Polar: SASAClass = ...
    Unclassified: SASAClass = ...
    names: dict[str, SASAClass] = ...
    values: dict[int, SASAClass] = ...
    __slots__: ClassVar[tuple] = ...

class SASAClassifier(Boost.Python.enum):
    NACCESS: SASAClassifier = ...
    OONS: SASAClassifier = ...
    Protor: SASAClassifier = ...
    names: dict[str, SASAClassifier] = ...
    values: dict[int, SASAClassifier] = ...
    __slots__: ClassVar[tuple] = ...

class SASAOpts(Boost.Python.instance):
    """
    Constructor takes no arguments

    C++ signature :void __init__(_object*)

    __init__( (object)arg1, (SASAAlgorithm)arg2, (SASAClassifier)arg3) -> None :

    C++ signature :void __init__(_object*,FreeSASA::SASAOpts::Algorithm,FreeSASA::SASAOpts::Classifier)

    __init__( (object)arg1, (SASAAlgorithm)arg2, (SASAClassifier)arg3, (float)arg4) -> None :

    C++ signature :void __init__(_object*,FreeSASA::SASAOpts::Algorithm,FreeSASA::SASAOpts::Classifier,double)

    property algorithm¶

    property classifier¶

    property probeRadius¶"""

    ...
    __instance_size__: ClassVar[int] = ...
    algorithm: Any
    classifier: Any
    probeRadius: Any

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

def CalcSASA(
    self,
    mol: Mol,
    radii: AtomPairsParameters,
    confIdx: int = -1,
    query: Atom = None,
    opts: SASAOpts = ...,
) -> float:
    """
    Compute the Solvent Accessible Surface Area using the FreeSASA library
    ARGUMENTS:

    mol: The molecule to compute.

    radii:  A list of atom raddii where radii[atom.GetIdx()] is the radius of the atomThese can be passed in or calculated with classifyAtoms for some proteins

    confIdx: Specify the conformer to use for the 3D geometry  [default -1]

    query: Pass along a query atom to compute the SASA for a subset of atoms.precanned query atoms can be made with MakeFreeSasaPolarAtomQuery and
    MakeFreeSasaAPolarAtomQuery for classified polar and apolar atoms respectively.

    opts: a SASAOpts class specifying the algorithm to use

    RETURNS:
    The computed solvent accessible surface area.

    C++ signature :double CalcSASA(RDKit::ROMol,boost::python::api::object [,int=-1 [,RDKit::Atom const*=None [,FreeSASA::SASAOpts=<rdkit.Chem.rdFreeSASA.SASAOpts object at 0x7fd4abb89a50>]]])
    """
    ...

def MakeFreeSasaAPolarAtomQuery(self) -> QueryAtom:
    """
    Returns an APolar atom query for use with CalcSASA.  An apolar atom has the SASAClass
    and SASAClassName set to the APOLAR class.  (see classifyAtoms)

    C++ signature :RDKit::QueryAtom const* MakeFreeSasaAPolarAtomQuery()"""
    ...

def MakeFreeSasaPolarAtomQuery(self) -> QueryAtom:
    """
    Returns a polar atom query for use with CalcSASA.  An polar atom has the SASAClass
    and SASAClassName set to the POLAR class.  (see classifyAtoms)

    C++ signature :RDKit::QueryAtom const* MakeFreeSasaPolarAtomQuery()"""
    ...

def classifyAtoms(self, mol: Mol, options: SASAOpts = ...) -> object:
    """
    Classify the atoms in the molecule returning their radii if possible.
    ARGUMENTS:

    mol: molecule to classify

    options: FreeSASA options class specifying the classification method.Current classifiers are Protor, NACCESS and OONS
    classification is stored as atom property ‘SASAClass’ for the integer value

    and ‘SASAClassName’ for the string name of the class, Polar, APolar…

    RETURNS:list of radii where radii[atom.GetIdx()] is the radii of the atom.
    If classification fails, NONE is returned

    C++ signature :boost::python::api::object classifyAtoms(RDKit::ROMol {lvalue} [,FreeSASA::SASAOpts=<rdkit.Chem.rdFreeSASA.SASAOpts object at 0x7fd4abb89970>])
    """
    ...
