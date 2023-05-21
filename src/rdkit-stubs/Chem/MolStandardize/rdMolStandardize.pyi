"""
rdkit.Chem.MolStandardize.rdMolStandardize module¶
Module containing tools for normalizing molecules defined by SMARTS patterns
"""
from typing import Any, ClassVar, overload

import Boost.Python
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdMolDescriptors import AtomPairsParameters

class AllowedAtomsValidation(Boost.Python.instance):
    """
    C++ signature :void* __init__(boost::python::api::object,boost::python::api::object)
    """

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def validate(
        self: AllowedAtomsValidation, mol: Mol, reportAllFailures: bool = False
    ) -> list:
        """
        C++ signature :boost::python::list validate(RDKit::MolStandardize::AllowedAtomsValidation {lvalue},RDKit::ROMol [,bool=False])
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class ChargeCorrection(Boost.Python.instance):
    """
    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,int)

    property Charge¶

    property Name¶

    property Smarts¶"""

    ...
    __instance_size__: ClassVar[int] = ...
    Charge: Any
    Name: Any
    Smarts: Any

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class CleanupParameters(Boost.Python.instance):
    """
    Parameters controlling molecular standardization

    C++ signature :void __init__(_object*)

    property acidbaseFile¶
    file containing the acid and base definitions

    property doCanonical¶
    apply atom-order dependent normalizations (like uncharging) in a canonical order

    property fragmentFile¶
    file containing the acid and base definitions

    property largestFragmentChooserCountHeavyAtomsOnly¶
    whether LargestFragmentChooser should only count heavy atoms (defaults to False)

    property largestFragmentChooserUseAtomCount¶
    Whether LargestFragmentChooser should use atom count as main criterion before MW (defaults to True)

    property maxRestarts¶
    maximum number of restarts

    property maxTautomers¶
    maximum number of tautomers to generate (defaults to 1000)

    property maxTransforms¶
    maximum number of transforms to apply during tautomer enumeration (defaults to 1000)

    property normalizationsFile¶
    file containing the normalization transformations

    property preferOrganic¶
    prefer organic fragments to inorganic ones when deciding what to keep

    property tautomerReassignStereo¶
    call AssignStereochemistry on all generated tautomers (defaults to True)

    property tautomerRemoveBondStereo¶
    remove stereochemistry from double bonds involved in tautomerism (defaults to True)

    property tautomerRemoveIsotopicHs¶
    remove isotopic Hs from centers involved in tautomerism (defaults to True)

    property tautomerRemoveSp3Stereo¶
    remove stereochemistry from sp3 centers involved in tautomerism (defaults to True)

    property tautomerTransformsFile¶
    file containing the tautomer transformations"""

    ...
    __instance_size__: ClassVar[int] = ...
    acidbaseFile: Any
    doCanonical: Any
    fragmentFile: Any
    largestFragmentChooserCountHeavyAtomsOnly: Any
    largestFragmentChooserUseAtomCount: Any
    maxRestarts: Any
    maxTautomers: Any
    maxTransforms: Any
    normalizationsFile: Any
    preferOrganic: Any
    tautomerReassignStereo: Any
    tautomerRemoveBondStereo: Any
    tautomerRemoveIsotopicHs: Any
    tautomerRemoveSp3Stereo: Any
    tautomerTransformsFile: Any

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class DisallowedAtomsValidation(Boost.Python.instance):
    """
    C++ signature :void* __init__(boost::python::api::object,boost::python::api::object)
    """

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def validate(
        self: DisallowedAtomsValidation, mol: Mol, reportAllFailures: bool = False
    ) -> list:
        """
        C++ signature :boost::python::list validate(RDKit::MolStandardize::DisallowedAtomsValidation {lvalue},RDKit::ROMol [,bool=False])
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class FragmentRemover(Boost.Python.instance):
    """
    C++ signature :void __init__(_object*)

    __init__( (object)arg1 [, (str)fragmentFilename=’’ [, (bool)leave_last=True [, (bool)skip_if_all_match=False]]]) -> None :

    C++ signature :void __init__(_object* [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’’ [,bool=True [,bool=False]]])
    """

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def remove(self: FragmentRemover, mol: Mol) -> Mol:
        """
        C++ signature :RDKit::ROMol* remove(RDKit::MolStandardize::FragmentRemover {lvalue},RDKit::ROMol)
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class FragmentValidation(MolVSValidations):
    """
    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def run(
        self: FragmentValidation, mol: Mol, reportAllFailures: bool, errors: object
    ) -> None:
        """
        C++ signature :void run(RDKit::MolStandardize::FragmentValidation {lvalue},RDKit::ROMol,bool,std::vector<RDKit::MolStandardize::ValidationErrorInfo, std::allocator<RDKit::MolStandardize::ValidationErrorInfo> > {lvalue})
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class IsotopeValidation(MolVSValidations):
    """
    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def run(
        self: IsotopeValidation, mol: Mol, reportAllFailures: bool, errors: object
    ) -> None:
        """
        C++ signature :void run(RDKit::MolStandardize::IsotopeValidation {lvalue},RDKit::ROMol,bool,std::vector<RDKit::MolStandardize::ValidationErrorInfo, std::allocator<RDKit::MolStandardize::ValidationErrorInfo> > {lvalue})
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class LargestFragmentChooser(Boost.Python.instance):
    """
    C++ signature :void __init__(_object* [,bool=False])

    __init__( (object)arg1, (CleanupParameters)params) -> None :

    C++ signature :void __init__(_object*,RDKit::MolStandardize::CleanupParameters)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def choose(self: LargestFragmentChooser, mol: Mol) -> Mol:
        """
        C++ signature :RDKit::ROMol* choose(RDKit::MolStandardize::LargestFragmentChooser {lvalue},RDKit::ROMol)
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class MetalDisconnector(Boost.Python.instance):
    """
    a class to disconnect metals that are defined as covalently bonded to non-metals

    C++ signature :void __init__(_object* [,boost::python::api::object=None])"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def Disconnect(self: MetalDisconnector, mol: Mol) -> Mol:
        """
        a class to disconnect metals that are defined as covalently bonded to non-metals

        C++ signature :RDKit::ROMol* Disconnect((anonymous namespace)::MetalDisconnectorWrap {lvalue},RDKit::ROMol)
        """
        ...
    def SetMetalNof(self: MetalDisconnector, mol: Mol) -> None:
        """
        Set the query molecule defining the metals to disconnect if attached to Nitrogen, Oxygen or Fluorine.

        C++ signature :void SetMetalNof((anonymous namespace)::MetalDisconnectorWrap {lvalue},RDKit::ROMol)
        """
        ...
    def SetMetalNon(self: MetalDisconnector, mol: Mol) -> None:
        """
        Set the query molecule defining the metals to disconnect from other inorganic elements.

        C++ signature :void SetMetalNon((anonymous namespace)::MetalDisconnectorWrap {lvalue},RDKit::ROMol)
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @property
    def MetalNof(self) -> Any: ...
    @property
    def MetalNon(self) -> Any: ...

class MetalDisconnectorOptions(Boost.Python.instance):
    """
    Metal Disconnector Options

    C++ signature :void __init__(_object*)

    property adjustCharges¶
    Whether to adjust charges on ligand atoms.  Default true.

    property removeHapticDummies¶
    Whether to remove the dummy atoms representing haptic bonds.  Such dummies are bonded to the metal with a bond that has the MolFileBondEndPts prop set.  Default false.

    property splitAromaticC¶
    Whether to split metal-aromatic C bonds.  Default false.

    property splitGrignards¶
    Whether to split Grignard-type complexes. Default false."""

    ...
    __instance_size__: ClassVar[int] = ...
    adjustCharges: Any
    removeHapticDummies: Any
    splitAromaticC: Any
    splitGrignards: Any

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class MolVSValidation(Boost.Python.instance):
    """
    C++ signature :void __init__(_object*)

    __init__( (AtomPairsParameters)arg1, (AtomPairsParameters)arg2) -> object :

    C++ signature :void* __init__(boost::python::api::object,boost::python::api::object)
    """

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def validate(
        self: MolVSValidation, mol: Mol, reportAllFailures: bool = False
    ) -> list:
        """
        C++ signature :boost::python::list validate(RDKit::MolStandardize::MolVSValidation {lvalue},RDKit::ROMol [,bool=False])
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class MolVSValidations(Boost.Python.instance):
    """
    Raises an exception
    This class cannot be instantiated from Python"""

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def run(
        self: MolVSValidations, mol: Mol, reportAllFailures: bool, errors: object
    ) -> None:
        """
        C++ signature :void run(RDKit::MolStandardize::MolVSValidations {lvalue},RDKit::ROMol,bool,std::vector<RDKit::MolStandardize::ValidationErrorInfo, std::allocator<RDKit::MolStandardize::ValidationErrorInfo> > {lvalue})
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class NeutralValidation(MolVSValidations):
    """
    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def run(
        self: NeutralValidation, mol: Mol, reportAllFailures: bool, errors: object
    ) -> None:
        """
        C++ signature :void run(RDKit::MolStandardize::NeutralValidation {lvalue},RDKit::ROMol,bool,std::vector<RDKit::MolStandardize::ValidationErrorInfo, std::allocator<RDKit::MolStandardize::ValidationErrorInfo> > {lvalue})
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class NoAtomValidation(MolVSValidations):
    """
    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def run(
        self: NoAtomValidation, mol: Mol, reportAllFailures: bool, errors: object
    ) -> None:
        """
        C++ signature :void run(RDKit::MolStandardize::NoAtomValidation {lvalue},RDKit::ROMol,bool,std::vector<RDKit::MolStandardize::ValidationErrorInfo, std::allocator<RDKit::MolStandardize::ValidationErrorInfo> > {lvalue})
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class Normalizer(Boost.Python.instance):
    """
    C++ signature :void __init__(_object*)

    __init__( (object)arg1, (str)normalizeFilename, (int)maxRestarts) -> None :

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,unsigned int)
    """

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def normalize(self: Normalizer, mol: Mol) -> Mol:
        """
        C++ signature :RDKit::ROMol* normalize(RDKit::MolStandardize::Normalizer {lvalue},RDKit::ROMol)
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class RDKitValidation(Boost.Python.instance):
    """
    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def validate(
        self: RDKitValidation, mol: Mol, reportAllFailures: bool = False
    ) -> list:
        """
        C++ signature :boost::python::list validate(RDKit::MolStandardize::RDKitValidation {lvalue},RDKit::ROMol [,bool=False])
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class Reionizer(Boost.Python.instance):
    """
    C++ signature :void __init__(_object*)

    __init__( (object)arg1, (str)arg2) -> None :

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)

    __init__( (object)arg1, (str)arg2, (object)arg3) -> None :

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,std::vector<RDKit::MolStandardize::ChargeCorrection, std::allocator<RDKit::MolStandardize::ChargeCorrection> >)
    """

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def reionize(self: Reionizer, mol: Mol) -> Mol:
        """
        C++ signature :RDKit::ROMol* reionize(RDKit::MolStandardize::Reionizer {lvalue},RDKit::ROMol)
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class SmilesTautomerMap(Boost.Python.instance):
    """
    maps SMILES strings to the respective Tautomer objects
    Raises an exception
    This class cannot be instantiated from Python"""

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def items(self, arg1: SmilesTautomerMap) -> tuple:
        """
        C++ signature :boost::python::tuple items(std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, RDKit::MolStandardize::Tautomer, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, RDKit::MolStandardize::Tautomer> > >)
        """
        ...
    def keys(self, arg1: SmilesTautomerMap) -> tuple:
        """
        C++ signature :boost::python::tuple keys(std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, RDKit::MolStandardize::Tautomer, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, RDKit::MolStandardize::Tautomer> > >)
        """
        ...
    def values(self, arg1: SmilesTautomerMap) -> tuple:
        """
        C++ signature :boost::python::tuple values(std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, RDKit::MolStandardize::Tautomer, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, RDKit::MolStandardize::Tautomer> > >)
        """
        ...
    @classmethod
    def __contains__(cls, other) -> Any: ...
    @classmethod
    def __delitem__(cls, other) -> Any: ...
    @classmethod
    def __getitem__(cls, index) -> Any: ...
    @classmethod
    def __iter__(cls) -> Any: ...
    @classmethod
    def __len__(cls) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @classmethod
    def __setitem__(cls, index, object) -> Any: ...

class Tautomer(Boost.Python.instance):
    """
    used to hold the aromatic and kekulized versions of each tautomer
    Raises an exception
    This class cannot be instantiated from Python

    property kekulized¶
    kekulized version of the tautomer

    property tautomer¶
    aromatic version of the tautomer"""

    ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @property
    def kekulized(self) -> Any: ...
    @property
    def tautomer(self) -> Any: ...

class TautomerEnumerator(Boost.Python.instance):
    """
    C++ signature :void* __init__(boost::python::api::object)

    __init__( (AtomPairsParameters)arg1, (CleanupParameters)arg2) -> object :

    C++ signature :void* __init__(boost::python::api::object,RDKit::MolStandardize::CleanupParameters)

    __init__( (AtomPairsParameters)arg1, (TautomerEnumerator)arg2) -> object :

    C++ signature :void* __init__(boost::python::api::object,RDKit::MolStandardize::TautomerEnumerator)
    """

    tautomerScoreVersion: property[property[str]] = ...

    @overload
    @classmethod
    def __init__(cls, boost) -> Any: ...
    @overload
    @classmethod
    def __init__(cls, boost, RDKit) -> Any: ...
    @overload
    @classmethod
    def __init__(cls, boost, RDKit) -> Any: ...
    @overload
    def Canonicalize(self: TautomerEnumerator, mol: Mol) -> Mol:
        """
        picks the canonical tautomer from an iterable of molecules using a custom scoring function

        C++ signature :RDKit::ROMol* Canonicalize(RDKit::MolStandardize::TautomerEnumerator,RDKit::ROMol,boost::python::api::object)
        """
        ...
    @overload
    def Canonicalize(
        self: TautomerEnumerator, mol: Mol, scoreFunc: AtomPairsParameters
    ) -> Mol: ...
    def Enumerate(self: TautomerEnumerator, mol: Mol) -> TautomerEnumeratorResult:
        """
        Generates the tautomers for a molecule.

        The enumeration rules are inspired by the publication:
        M. Sitzmann et al., “Tautomerism in Large Databases.”, JCAMD 24:521 (2010)
        https://doi.org/10.1007/s10822-010-9346-4
        Note: the definitions used here are that the atoms modified during
        tautomerization are the atoms at the beginning and end of each tautomer
        transform (the H “donor” and H “acceptor” in the transform) and the bonds
        modified during transformation are any bonds whose order is changed during
        the tautomer transform (these are the bonds between the “donor” and the
        “acceptor”).

        C++ signature :(anonymous namespace)::PyTautomerEnumeratorResult* Enumerate(RDKit::MolStandardize::TautomerEnumerator,RDKit::ROMol)
        """
        ...
    def GetCallback(self, arg1: TautomerEnumerator) -> object:
        """
        Get the TautomerEnumeratorCallback subclass instance,
        or None if none was set.

        C++ signature :boost::python::api::object GetCallback(RDKit::MolStandardize::TautomerEnumerator)
        """
        ...
    def GetMaxTautomers(self: TautomerEnumerator) -> int:
        """
        returns the maximum number of tautomers to be generated.

        C++ signature :unsigned int GetMaxTautomers(RDKit::MolStandardize::TautomerEnumerator {lvalue})
        """
        ...
    def GetMaxTransforms(self: TautomerEnumerator) -> int:
        """
        returns the maximum number of transformations to be applied.

        C++ signature :unsigned int GetMaxTransforms(RDKit::MolStandardize::TautomerEnumerator {lvalue})
        """
        ...
    def GetReassignStereo(self: TautomerEnumerator) -> bool:
        """
        returns whether AssignStereochemistry will be called on each tautomer generated by the Enumerate() method.

        C++ signature :bool GetReassignStereo(RDKit::MolStandardize::TautomerEnumerator {lvalue})
        """
        ...
    def GetRemoveBondStereo(self: TautomerEnumerator) -> bool:
        """
        returns whether stereochemistry information will be removed from double bonds involved in tautomerism.

        C++ signature :bool GetRemoveBondStereo(RDKit::MolStandardize::TautomerEnumerator {lvalue})
        """
        ...
    def GetRemoveSp3Stereo(self: TautomerEnumerator) -> bool:
        """
        returns whether stereochemistry information will be removed from sp3 atoms involved in tautomerism.

        C++ signature :bool GetRemoveSp3Stereo(RDKit::MolStandardize::TautomerEnumerator {lvalue})
        """
        ...
    @overload
    def PickCanonical(self: TautomerEnumerator, iterable: AtomPairsParameters) -> Mol:
        """
        returns the canonical tautomer for a molecule using a custom scoring function

        C++ signature :RDKit::ROMol* PickCanonical(RDKit::MolStandardize::TautomerEnumerator,boost::python::api::object,boost::python::api::object)
        """
        ...
    @overload
    def PickCanonical(
        self: TautomerEnumerator,
        iterable: AtomPairsParameters,
        scoreFunc: AtomPairsParameters,
    ) -> Mol: ...
    def ScoreTautomer(self, mol: Mol) -> int:
        """
        returns the score for a tautomer using the default scoring scheme.

        C++ signature :int ScoreTautomer(RDKit::ROMol)"""
        ...
    def SetCallback(self, arg1: TautomerEnumerator, arg2: object) -> None:
        """
        Pass an instance of a class derived from
        TautomerEnumeratorCallback, which must implement the
        __call__() method.

        C++ signature :void SetCallback(RDKit::MolStandardize::TautomerEnumerator {lvalue},_object*)
        """
        ...
    def SetMaxTautomers(self: TautomerEnumerator, maxTautomers: int) -> None:
        """
        set the maximum number of tautomers to be generated.

        C++ signature :void SetMaxTautomers(RDKit::MolStandardize::TautomerEnumerator {lvalue},unsigned int)
        """
        ...
    def SetMaxTransforms(self: TautomerEnumerator, maxTransforms: int) -> None:
        """
        set the maximum number of transformations to be applied. This limit is usually hit earlier than the maxTautomers limit and leads to a more linear scaling of CPU time with increasing number of tautomeric centers (see Sitzmann et al.).

        C++ signature :void SetMaxTransforms(RDKit::MolStandardize::TautomerEnumerator {lvalue},unsigned int)
        """
        ...
    def SetReassignStereo(self: TautomerEnumerator, reassignStereo: bool) -> None:
        """
        set to True if you wish AssignStereochemistry to be called on each tautomer generated by the Enumerate() method. This defaults to True.

        C++ signature :void SetReassignStereo(RDKit::MolStandardize::TautomerEnumerator {lvalue},bool)
        """
        ...
    def SetRemoveBondStereo(self: TautomerEnumerator, removeBondStereo: bool) -> None:
        """
        set to True if you wish stereochemistry information to be removed from double bonds involved in tautomerism. This means that enols will lose their E/Z stereochemistry after going through tautomer enumeration because of the keto-enolic tautomerism. This defaults to True in the RDKit and also in the workflow described by Sitzmann et al.

        C++ signature :void SetRemoveBondStereo(RDKit::MolStandardize::TautomerEnumerator {lvalue},bool)
        """
        ...
    def SetRemoveSp3Stereo(self: TautomerEnumerator, removeSp3Stereo: bool) -> None:
        """
        set to True if you wish stereochemistry information to be removed from sp3 atoms involved in tautomerism. This means that S-aminoacids will lose their stereochemistry after going through tautomer enumeration because of the amido-imidol tautomerism. This defaults to True in RDKit, and to False in the workflow described by Sitzmann et al.

        C++ signature :void SetRemoveSp3Stereo(RDKit::MolStandardize::TautomerEnumerator {lvalue},bool)
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    tautomerScoreVersion: property[property[str]] = ...

class TautomerEnumeratorCallback(Boost.Python.instance):
    """
    Create a derived class from this abstract base class and
    implement the __call__() method.
    The __call__() method is called in the innermost loop of the
    algorithm, and provides a mechanism to monitor or stop
    its progress.
    To have your callback called, pass an instance of your
    derived class to TautomerEnumerator.SetCallback()

    C++ signature :void __init__(_object*)"""

    ...
    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __call__(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class TautomerEnumeratorResult(Boost.Python.instance):
    """
    used to return tautomer enumeration results
    Raises an exception
    This class cannot be instantiated from Python

    property modifiedAtoms¶
    tuple of atom indices modified by the transforms

    property modifiedBonds¶
    tuple of bond indices modified by the transforms

    property smiles¶
    SMILES of tautomers generated by the enumerator

    property smilesTautomerMap¶
    dictionary mapping SMILES strings to the respective Tautomer objects

    property status¶
    whether the enumeration completed or not; see TautomerEnumeratorStatus for possible values

    property tautomers¶
    tautomers generated by the enumerator"""

    ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __call__(cls, anonymousnamespace) -> Any: ...
    @classmethod
    def __getitem__(cls, anonymousnamespace, int) -> Any: ...
    @classmethod
    def __iter__(cls, boost) -> Any: ...
    @classmethod
    def __len__(cls, anonymousnamespace) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @property
    def modifiedAtoms(self) -> Any: ...
    @property
    def modifiedBonds(self) -> Any: ...
    @property
    def smiles(self) -> Any: ...
    @property
    def smilesTautomerMap(self) -> Any: ...
    @property
    def status(self) -> Any: ...
    @property
    def tautomers(self) -> Any: ...

class TautomerEnumeratorStatus(Boost.Python.enum):
    Canceled: TautomerEnumeratorStatus = ...
    Completed: TautomerEnumeratorStatus = ...
    MaxTautomersReached: TautomerEnumeratorStatus = ...
    MaxTransformsReached: TautomerEnumeratorStatus = ...
    names: dict[str, TautomerEnumeratorStatus] = ...
    values: dict[int, TautomerEnumeratorStatus] = ...
    __slots__: ClassVar[tuple] = ...

class Uncharger(Boost.Python.instance):
    """
    C++ signature :void __init__(_object* [,bool=True])"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def uncharge(self: Uncharger, mol: Mol) -> Mol:
        """
        C++ signature :RDKit::ROMol* uncharge(RDKit::MolStandardize::Uncharger {lvalue},RDKit::ROMol)
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class map_indexing_suite_SmilesTautomerMap_entry(Boost.Python.instance):
    """
    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def data(self, arg1: map_indexing_suite_SmilesTautomerMap_entry) -> Tautomer:
        """
        C++ signature :RDKit::MolStandardize::Tautomer data(std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, RDKit::MolStandardize::Tautomer> {lvalue})
        """
        ...
    def key(self, arg1: map_indexing_suite_SmilesTautomerMap_entry) -> str:
        """
        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > key(std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, RDKit::MolStandardize::Tautomer> {lvalue})
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

def CHARGE_CORRECTIONS(self) -> object:
    """
    C++ signature :std::vector<RDKit::MolStandardize::ChargeCorrection, std::allocator<RDKit::MolStandardize::ChargeCorrection> > CHARGE_CORRECTIONS()
    """
    ...

def CanonicalTautomer(self, mol: Mol, params: AtomPairsParameters = None) -> Mol:
    """
    Returns the canonical tautomer for the molecule

    C++ signature :RDKit::ROMol* CanonicalTautomer(RDKit::ROMol const* [,boost::python::api::object=None])
    """
    ...

def ChargeParent(
    self, mol: Mol, params: AtomPairsParameters = None, skipStandardize: bool = False
) -> Mol:
    """
    Returns the uncharged version of the largest fragment

    C++ signature :RDKit::ROMol* ChargeParent(RDKit::ROMol const* [,boost::python::api::object=None [,bool=False]])
    """
    ...

def Cleanup(self, mol: Mol, params: AtomPairsParameters = None) -> Mol:
    """
    Standardizes a molecule

    C++ signature :RDKit::ROMol* Cleanup(RDKit::ROMol const* [,boost::python::api::object=None])
    """
    ...

def DisconnectOrganometallics(
    self, mol: Mol, params: AtomPairsParameters = None
) -> Mol:
    """
    Returns the molecule disconnected using the organometallics rules.

    C++ signature :RDKit::ROMol* DisconnectOrganometallics(RDKit::ROMol {lvalue} [,boost::python::api::object=None])
    """
    ...

def FragmentParent(
    self, mol: Mol, params: AtomPairsParameters = None, skipStandardize: bool = False
) -> Mol:
    """
    Returns the largest fragment after doing a cleanup

    C++ signature :RDKit::ROMol* FragmentParent(RDKit::ROMol const* [,boost::python::api::object=None [,bool=False]])
    """
    ...

def FragmentRemoverFromData(
    self, fragmentData: str, leave_last: bool = True, skip_if_all_match: bool = False
) -> FragmentRemover:
    """
    creates a FragmentRemover from a string containing parameter data

    C++ signature :RDKit::MolStandardize::FragmentRemover* FragmentRemoverFromData(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=True [,bool=False]])
    """
    ...

def GetV1TautomerEnumerator(self) -> TautomerEnumerator:
    """
    return a TautomerEnumerator using v1 of the enumeration rules

    C++ signature :RDKit::MolStandardize::TautomerEnumerator* GetV1TautomerEnumerator()
    """
    ...

def IsotopeParent(
    self, mol: Mol, params: AtomPairsParameters = None, skipStandardize: bool = False
) -> Mol:
    """
    removes all isotopes specifications from the given molecule

    C++ signature :RDKit::ROMol* IsotopeParent(RDKit::ROMol const* [,boost::python::api::object=None [,bool=False]])
    """
    ...

def Normalize(self, mol: Mol, params: AtomPairsParameters = None) -> Mol:
    """
    Applies a series of standard transformations to correct functional groups and recombine charges

    C++ signature :RDKit::ROMol* Normalize(RDKit::ROMol const* [,boost::python::api::object=None])
    """
    ...

def NormalizerFromData(self, paramData: str, params: CleanupParameters) -> Normalizer:
    """
    creates a Normalizer from a string containing normalization SMARTS

    C++ signature :RDKit::MolStandardize::Normalizer* NormalizerFromData(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,RDKit::MolStandardize::CleanupParameters)
    """
    ...

def NormalizerFromParams(self, params: CleanupParameters) -> Normalizer:
    """
    creates a Normalizer from CleanupParameters

    C++ signature :RDKit::MolStandardize::Normalizer* NormalizerFromParams(RDKit::MolStandardize::CleanupParameters)
    """
    ...

def Reionize(self, mol: Mol, params: AtomPairsParameters = None) -> Mol:
    """
    Ensures the strongest acid groups are charged first

    C++ signature :RDKit::ROMol* Reionize(RDKit::ROMol const* [,boost::python::api::object=None])
    """
    ...

def ReionizerFromData(
    self, paramData: str, chargeCorrections: AtomPairsParameters = []
) -> Reionizer:
    """
    creates a reionizer from a string containing parameter data and a list of charge corrections

    C++ signature :RDKit::MolStandardize::Reionizer* ReionizerFromData(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,boost::python::api::object=[]])
    """
    ...

def RemoveFragments(self, mol: Mol, params: AtomPairsParameters = None) -> Mol:
    """
    Removes fragments from the molecule

    C++ signature :RDKit::ROMol* RemoveFragments(RDKit::ROMol const* [,boost::python::api::object=None])
    """
    ...

def StandardizeSmiles(self, smiles: str) -> str:
    """
    Convenience function for standardizing a SMILES

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > StandardizeSmiles(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
    ...

def StereoParent(
    self, mol: Mol, params: AtomPairsParameters = None, skipStandardize: bool = False
) -> Mol:
    """
    calls removeStereochemistry() on the given molecule

    C++ signature :RDKit::ROMol* StereoParent(RDKit::ROMol const* [,boost::python::api::object=None [,bool=False]])
    """
    ...

def SuperParent(
    self, mol: Mol, params: AtomPairsParameters = None, skipStandardize: bool = False
) -> Mol:
    """
    Returns the super parent. The super parent is the fragment, charge, isotope, stereo, and tautomer parent of the molecule.

    C++ signature :RDKit::ROMol* SuperParent(RDKit::ROMol const* [,boost::python::api::object=None [,bool=False]])
    """
    ...

def TautomerParent(
    self, mol: Mol, params: AtomPairsParameters = None, skipStandardize: bool = False
) -> Mol:
    """
    Returns the tautomer parent of a given molecule. The fragment parent is the standardized canonical tautomer of the molecule

    C++ signature :RDKit::ROMol* TautomerParent(RDKit::ROMol const* [,boost::python::api::object=None [,bool=False]])
    """
    ...

def UpdateParamsFromJSON(self, arg1: CleanupParameters, arg2: str) -> None:
    """
    updates the cleanup parameters from the provided JSON string

    C++ signature :void UpdateParamsFromJSON(RDKit::MolStandardize::CleanupParameters {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
    ...

def ValidateSmiles(self, mol: str) -> list:
    """
    C++ signature :boost::python::list ValidateSmiles(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
    ...
