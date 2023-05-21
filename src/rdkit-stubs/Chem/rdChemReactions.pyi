"""
rdkit.Chem.rdChemReactions module¶
Module containing classes and functions for working with chemical reactions.
"""
from typing import Any, ClassVar, overload

import Boost.Python
from rdkit.Chem import rdmolops
from rdkit.Chem.rdchem import Mol, SubstructMatchParameters
from rdkit.Chem.rdMolDescriptors import AtomPairsParameters
from rdkit.Chem.rdmolops import AdjustQueryParameters
from rdkit.DataStructs.cDataStructs import ExplicitBitVect, UIntSparseIntVect
from rdkit.rdBase import _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE

SANITIZE_ADJUST_REACTANTS: SanitizeFlags
SANITIZE_ALL: SanitizeFlags
SANITIZE_ATOM_MAPS: SanitizeFlags
SANITIZE_MERGEHS: SanitizeFlags
SANITIZE_NONE: SanitizeFlags
SANITIZE_RGROUP_NAMES: SanitizeFlags

class CartesianProductStrategy(EnumerationStrategyBase):
    """
    CartesianProductStrategy produces a standard walk through all possible
    reagent combinations:
    (0,0,0), (1,0,0), (2,0,0) …

    C++ signature :void __init__(_object*)"""

    ...
    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __copy__(cls, RDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class ChemicalReaction(Boost.Python.instance):
    """
    A class for storing and applying chemical reactions.

    Sample Usage:>>> from rdkit import Chem
    >>> from rdkit.Chem import rdChemReactions
    >>> rxn = rdChemReactions.ReactionFromSmarts('[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]')
    >>> reacts = (Chem.MolFromSmiles('C(=O)O'),Chem.MolFromSmiles('CNC'))
    >>> products = rxn.RunReactants(reacts)
    >>> len(products)
    1
    >>> len(products[0])
    1
    >>> Chem.MolToSmiles(products[0][0])
    'CN(C)C=O'

    Constructor, takes no arguments

    C++ signature :void __init__(_object*)

    __init__( (object)arg1, (str)arg2) -> None :

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)

    __init__( (object)arg1, (ChemicalReaction)arg2) -> None :

    C++ signature :void __init__(_object*,RDKit::ChemicalReaction)"""

    __getstate_manages_dict__: ClassVar[bool] = ...
    __instance_size__: ClassVar[int] = ...
    __safe_for_unpickling__: ClassVar[bool] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def AddAgentTemplate(self, arg1: ChemicalReaction, arg2: Mol) -> int:
        """
        adds a agent (a Molecule)

        C++ signature :unsigned int AddAgentTemplate(RDKit::ChemicalReaction {lvalue},boost::shared_ptr<RDKit::ROMol>)
        """
        ...
    def AddProductTemplate(self, arg1: ChemicalReaction, arg2: Mol) -> int:
        """
        adds a product (a Molecule)

        C++ signature :unsigned int AddProductTemplate(RDKit::ChemicalReaction {lvalue},boost::shared_ptr<RDKit::ROMol>)
        """
        ...
    def AddReactantTemplate(self, arg1: ChemicalReaction, arg2: Mol) -> int:
        """
        adds a reactant (a Molecule) to the reaction

        C++ signature :unsigned int AddReactantTemplate(RDKit::ChemicalReaction {lvalue},boost::shared_ptr<RDKit::ROMol>)
        """
        ...
    def AddRecursiveQueriesToReaction(
        self,
        reaction: ChemicalReaction,
        queries: dict = {},
        propName: str = "molFileValue",
        getLabels: bool = False,
    ) -> object:
        """
        adds recursive queries and returns reactant labels

        C++ signature :boost::python::api::object AddRecursiveQueriesToReaction(RDKit::ChemicalReaction {lvalue} [,boost::python::dict={} [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’molFileValue’ [,bool=False]]])
        """
        ...
    def ClearComputedProps(self, arg1: ChemicalReaction) -> None:
        """
        Removes all computed properties from the reaction.

        C++ signature :void ClearComputedProps(RDKit::ChemicalReaction)"""
        ...
    def ClearProp(self, arg1: ChemicalReaction, arg2: str) -> None:
        """
        Removes a property from the reaction.

        ARGUMENTS:
        key: the name of the property to clear (a string).

        C++ signature :void ClearProp(RDKit::ChemicalReaction,char const*)"""
        ...
    def GetAgentTemplate(self: ChemicalReaction, which: int) -> Mol:
        """
        returns one of our agent templates

        C++ signature :RDKit::ROMol* GetAgentTemplate(RDKit::ChemicalReaction const*,unsigned int)
        """
        ...
    def GetAgents(self, arg1: ChemicalReaction) -> MOL_SPTR_VECT:
        """
        get the agent templates

        C++ signature :std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > > GetAgents(RDKit::ChemicalReaction {lvalue})
        """
        ...
    def GetBoolProp(self, arg1: ChemicalReaction, arg2: str) -> bool:
        """
        Returns the Bool value of the property if possible.

        ARGUMENTS:
        key: the name of the property to return (a string).

        RETURNS: a bool

        NOTE:
        If the property has not been set, a KeyError exception will be raised.

        C++ signature :bool GetBoolProp(RDKit::ChemicalReaction*,char const*)"""
        ...
    def GetDoubleProp(self, arg1: ChemicalReaction, arg2: str) -> float:
        """
        Returns the double value of the property if possible.

        ARGUMENTS:
        key: the name of the property to return (a string).

        RETURNS: a double

        NOTE:
        If the property has not been set, a KeyError exception will be raised.

        C++ signature :double GetDoubleProp(RDKit::ChemicalReaction*,char const*)"""
        ...
    def GetIntProp(self, arg1: ChemicalReaction, arg2: str) -> int:
        """
        Returns the integer value of the property if possible.

        ARGUMENTS:
        key: the name of the property to return (a string).

        RETURNS: an integer

        NOTE:
        If the property has not been set, a KeyError exception will be raised.

        C++ signature :int GetIntProp(RDKit::ChemicalReaction*,char const*)"""
        ...
    def GetNumAgentTemplates(self, arg1: ChemicalReaction) -> int:
        """
        returns the number of agents this reaction expects

        C++ signature :unsigned int GetNumAgentTemplates(RDKit::ChemicalReaction {lvalue})
        """
        ...
    def GetNumProductTemplates(self, arg1: ChemicalReaction) -> int:
        """
        returns the number of products this reaction generates

        C++ signature :unsigned int GetNumProductTemplates(RDKit::ChemicalReaction {lvalue})
        """
        ...
    def GetNumReactantTemplates(self, arg1: ChemicalReaction) -> int:
        """
        returns the number of reactants this reaction expects

        C++ signature :unsigned int GetNumReactantTemplates(RDKit::ChemicalReaction {lvalue})
        """
        ...
    def GetProductTemplate(self: ChemicalReaction, which: int) -> Mol:
        """
        returns one of our product templates

        C++ signature :RDKit::ROMol* GetProductTemplate(RDKit::ChemicalReaction const*,unsigned int)
        """
        ...
    def GetProducts(self, arg1: ChemicalReaction) -> MOL_SPTR_VECT:
        """
        get the product templates

        C++ signature :std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > > GetProducts(RDKit::ChemicalReaction {lvalue})
        """
        ...
    def GetProp(self, arg1: ChemicalReaction, arg2: str) -> str:
        """
        Returns the value of the property.

        ARGUMENTS:
        key: the name of the property to return (a string).

        RETURNS: a string

        NOTE:
        If the property has not been set, a KeyError exception will be raised.

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetProp(RDKit::ChemicalReaction*,char const*)
        """
        ...
    def GetPropNames(
        self: ChemicalReaction,
        includePrivate: bool = False,
        includeComputed: bool = False,
    ) -> _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE:
        """
        Returns a tuple with all property names for this reaction.

        ARGUMENTS:

        includePrivate: (optional) toggles inclusion of private properties in the result set.Defaults to 0.

        includeComputed: (optional) toggles inclusion of computed properties in the result set.Defaults to 0.

        RETURNS: a tuple of strings

        C++ signature :std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > GetPropNames(RDKit::ChemicalReaction {lvalue} [,bool=False [,bool=False]])
        """
        ...
    def GetPropsAsDict(
        self: ChemicalReaction,
        includePrivate: bool = False,
        includeComputed: bool = False,
    ) -> dict:
        """
        Returns a dictionary populated with the reaction’s properties.n.b. Some properties are not able to be converted to python types.

        ARGUMENTS:

        includePrivate: (optional) toggles inclusion of private properties in the result set.Defaults to False.

        includeComputed: (optional) toggles inclusion of computed properties in the result set.Defaults to False.

        RETURNS: a dictionary

        C++ signature :boost::python::dict GetPropsAsDict(RDKit::ChemicalReaction [,bool=False [,bool=False]])
        """
        ...
    def GetReactantTemplate(self: ChemicalReaction, which: int) -> Mol:
        """
        returns one of our reactant templates

        C++ signature :RDKit::ROMol* GetReactantTemplate(RDKit::ChemicalReaction const*,unsigned int)
        """
        ...
    def GetReactants(self, arg1: ChemicalReaction) -> MOL_SPTR_VECT:
        """
        get the reactant templates

        C++ signature :std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > > GetReactants(RDKit::ChemicalReaction {lvalue})
        """
        ...
    def GetReactingAtoms(
        self: ChemicalReaction, mappedAtomsOnly: bool = False
    ) -> object:
        """
        returns a sequence of sequences with the atoms that change in the reaction

        C++ signature :boost::python::api::object GetReactingAtoms(RDKit::ChemicalReaction [,bool=False])
        """
        ...
    def GetSubstructParams(self, arg1: ChemicalReaction) -> SubstructMatchParameters:
        """
        get the parameter object controlling the substructure matching

        C++ signature :RDKit::SubstructMatchParameters* GetSubstructParams(RDKit::ChemicalReaction {lvalue})
        """
        ...
    def GetUnsignedProp(self, arg1: ChemicalReaction, arg2: str) -> int:
        """
        Returns the unsigned int value of the property if possible.

        ARGUMENTS:
        key: the name of the property to return (a string).

        RETURNS: an unsigned integer

        NOTE:
        If the property has not been set, a KeyError exception will be raised.

        C++ signature :unsigned int GetUnsignedProp(RDKit::ChemicalReaction*,char const*)
        """
        ...
    def HasProp(self, arg1: ChemicalReaction, arg2: str) -> int:
        """
        Queries a molecule to see if a particular property has been assigned.

        ARGUMENTS:
        key: the name of the property to check for (a string).

        C++ signature :int HasProp(RDKit::ChemicalReaction,char const*)"""
        ...
    def Initialize(self: ChemicalReaction, silent: bool = False) -> None:
        """
        initializes the reaction so that it can be used

        C++ signature :void Initialize(RDKit::ChemicalReaction {lvalue} [,bool=False])
        """
        ...
    def IsInitialized(self, arg1: ChemicalReaction) -> bool:
        """
        checks if the reaction is ready for use

        C++ signature :bool IsInitialized(RDKit::ChemicalReaction {lvalue})"""
        ...
    def IsMoleculeAgent(self, arg1: ChemicalReaction, arg2: Mol) -> bool:
        """
        returns whether or not the molecule has a substructure match to one of the agents.

        C++ signature :bool IsMoleculeAgent(RDKit::ChemicalReaction,RDKit::ROMol)"""
        ...
    def IsMoleculeProduct(self, arg1: ChemicalReaction, arg2: Mol) -> bool:
        """
        returns whether or not the molecule has a substructure match to one of the products.

        C++ signature :bool IsMoleculeProduct(RDKit::ChemicalReaction,RDKit::ROMol)"""
        ...
    def IsMoleculeReactant(self, arg1: ChemicalReaction, arg2: Mol) -> bool:
        """
        returns whether or not the molecule has a substructure match to one of the reactants.

        C++ signature :bool IsMoleculeReactant(RDKit::ChemicalReaction,RDKit::ROMol)"""
        ...
    def RemoveAgentTemplates(
        self: ChemicalReaction, targetList: AtomPairsParameters = None
    ) -> None:
        """
        Removes agents from reaction. If targetList is provide the agents will be transferred to that list.

        C++ signature :void RemoveAgentTemplates(RDKit::ChemicalReaction {lvalue} [,boost::python::api::object=None])
        """
        ...
    def RemoveUnmappedProductTemplates(
        self: ChemicalReaction,
        thresholdUnmappedAtoms: float = 0.2,
        moveToAgentTemplates: bool = True,
        targetList: AtomPairsParameters = None,
    ) -> None:
        """
        Removes molecules with an atom mapping ratio below thresholdUnmappedAtoms from product templates to the agent templates or to a given targetList

        C++ signature :void RemoveUnmappedProductTemplates(RDKit::ChemicalReaction* [,double=0.2 [,bool=True [,boost::python::api::object=None]]])
        """
        ...
    def RemoveUnmappedReactantTemplates(
        self: ChemicalReaction,
        thresholdUnmappedAtoms: float = 0.2,
        moveToAgentTemplates: bool = True,
        targetList: AtomPairsParameters = None,
    ) -> None:
        """
        Removes molecules with an atom mapping ratio below thresholdUnmappedAtoms from reactant templates to the agent templates or to a given targetList

        C++ signature :void RemoveUnmappedReactantTemplates(RDKit::ChemicalReaction* [,double=0.2 [,bool=True [,boost::python::api::object=None]]])
        """
        ...
    def RunReactant(
        self, arg1: ChemicalReaction, arg2: AtomPairsParameters, arg3: int
    ) -> object:
        """
        apply the reaction to a single reactant

        C++ signature :_object* RunReactant(RDKit::ChemicalReaction*,boost::python::api::object,unsigned int)
        """
        ...
    def RunReactantInPlace(self: ChemicalReaction, reactant: Mol) -> bool:
        """
        apply the reaction to a single reactant in place. The reactant itself is modified. This can only be used for single reactant - single product reactions.

        C++ signature :bool RunReactantInPlace(RDKit::ChemicalReaction*,RDKit::ROMol*)
        """
        ...
    @overload
    def RunReactants(
        self: ChemicalReaction, reactants: tuple, maxProducts: int = 1000
    ) -> object:
        """
        apply the reaction to a sequence of reactant molecules and return the products as a tuple of tuples.  If maxProducts is not zero, stop the reaction when maxProducts have been generated [default=1000]

        C++ signature :_object* RunReactants(RDKit::ChemicalReaction*,boost::python::list [,unsigned int=1000])
        """
        ...
    @overload
    def RunReactants(
        self: ChemicalReaction, reactants: list, maxProducts: int = 1000
    ) -> object: ...
    def SetBoolProp(
        self: ChemicalReaction, key: str, val: bool, computed: bool = False
    ) -> None:
        """
        Sets a boolean valued molecular property

        ARGUMENTS:
        key: the name of the property to be set (a string).
        value: the property value as a bool.

        computed: (optional) marks the property as being computed.Defaults to False.

        C++ signature :void SetBoolProp(RDKit::ChemicalReaction,char const*,bool [,bool=False])
        """
        ...
    def SetDoubleProp(
        self: ChemicalReaction, key: str, val: float, computed: bool = False
    ) -> None:
        """
        Sets a double valued molecular property

        ARGUMENTS:
        key: the name of the property to be set (a string).
        value: the property value as a double.

        computed: (optional) marks the property as being computed.Defaults to 0.

        C++ signature :void SetDoubleProp(RDKit::ChemicalReaction,char const*,double [,bool=False])
        """
        ...
    def SetIntProp(
        self: ChemicalReaction, key: str, val: int, computed: bool = False
    ) -> None:
        """
        Sets an integer valued molecular property

        ARGUMENTS:
        key: the name of the property to be set (an unsigned number).
        value: the property value as an integer.

        computed: (optional) marks the property as being computed.Defaults to False.

        C++ signature :void SetIntProp(RDKit::ChemicalReaction,char const*,int [,bool=False])
        """
        ...
    def SetProp(
        self: ChemicalReaction, key: str, val: str, computed: bool = False
    ) -> None:
        """
        Sets a molecular property

        ARGUMENTS:
        key: the name of the property to be set (a string).
        value: the property value (a string).

        computed: (optional) marks the property as being computed.Defaults to False.

        C++ signature :void SetProp(RDKit::ChemicalReaction,char const*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])
        """
        ...
    def SetUnsignedProp(
        self: ChemicalReaction, key: str, val: int, computed: bool = False
    ) -> None:
        """
        Sets an unsigned integer valued molecular property

        ARGUMENTS:
        key: the name of the property to be set (a string).
        value: the property value as an unsigned integer.

        computed: (optional) marks the property as being computed.Defaults to False.

        C++ signature :void SetUnsignedProp(RDKit::ChemicalReaction,char const*,unsigned int [,bool=False])
        """
        ...
    @overload
    def ToBinary(self: ChemicalReaction) -> object:
        """
        Returns a binary string representation of the reaction.

        C++ signature :boost::python::api::object ToBinary(RDKit::ChemicalReaction,unsigned int)
        """
        ...
    @overload
    def ToBinary(self: ChemicalReaction, propertyFlags: int) -> object: ...
    def Validate(self: ChemicalReaction, silent: bool = False) -> tuple:
        """
        checks the reaction for potential problems, returns (numWarnings,numErrors)

        C++ signature :boost::python::tuple Validate(RDKit::ChemicalReaction const* [,bool=False])
        """
        ...
    @classmethod
    def _getImplicitPropertiesFlag(cls, RDKit) -> Any: ...
    @classmethod
    def _setImplicitPropertiesFlag(cls, RDKit, bool) -> Any: ...
    @classmethod
    def __getinitargs__(cls, RDKit) -> Any: ...
    @classmethod
    def __getstate__(cls, boost) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @classmethod
    def __setstate__(cls, state) -> Any: ...

class EnumerateLibrary(EnumerateLibraryBase):
    """
    This class allows easy enumeration of reactions.  Simply provide a reaction
    and a set of reagents and you are off the races.
    Note that this functionality should be considered beta and that the API may
    change in a future release.
    EnumerateLibrary follows the python enumerator protocol, for example:
    library = EnumerateLibrary(rxn, bbs)
    for products in library:

    … do something with the product

    It is useful to sanitize reactions before hand:
    SanitizeRxn(rxn)
    library = EnumerateLibrary(rxn, bbs)
    If ChemDraw style reaction semantics are prefereed, you can apply
    the ChemDraw parameters:
    SanitizeRxn(rxn, params=GetChemDrawRxnAdjustParams())
    For one, this enforces only matching RGroups and assumes all atoms
    have fully satisfied valences.
    Each product has the same output as applying a set of reagents to
    the libraries reaction.
    This can be a bit confusing as each product can have multiple molecules
    generated.  The returned data structure is as follows:

    [ [products1], [products2],… ]

    Where products1 are the molecule products for the reactions first product
    template and products2 are the molecule products for the second product
    template.  Since each reactant can match more than once, there may be
    multiple product molecules for each template.

    for products in library:
    for results_for_product_template in products:
    for mol in results_for_product_template:Chem.MolToSmiles(mol) # finally have a molecule!

    For sufficiently large libraries, using this iteration strategy is not
    recommended as the library may contain more products than atoms in the
    universe.  To help with this, you can supply an enumeration strategy.
    The default strategy is a CartesianProductStrategy which enumerates
    everything.  RandomSampleStrategy randomly samples the products but
    this strategy never terminates, however, python supplies itertools:
    import itertools
    library = EnumerateLibrary(rxn, bbs, rdChemReactions.RandomSampleStrategy())
    for result in itertools.islice(library, 1000):

    # do something with the first 1000 samples

    for result in itertools.islice(library, 1000):# do something with the next 1000 samples

    Libraries are also serializable, including their current state:
    s = library.Serialize()
    library2 = EnumerateLibrary()
    library2.InitFromString(s)
    for result in itertools.islice(libary2, 1000):

    # do something with the next 1000 samples

    __init__( (object)arg1) -> None :

    C++ signature :void __init__(_object*)

    __init__( (object)arg1, (ChemicalReaction)rxn, (list)reagents [, (EnumerationParams)params]) -> None :

    C++ signature :void __init__(_object*,RDKit::ChemicalReaction,boost::python::list [,RDKit::EnumerationParams])

    __init__( (object)arg1, (ChemicalReaction)rxn, (tuple)reagents [, (EnumerationParams)params]) -> None :

    C++ signature :void __init__(_object*,RDKit::ChemicalReaction,boost::python::tuple [,RDKit::EnumerationParams])

    __init__( (object)arg1, (ChemicalReaction)rxn, (list)reagents, (EnumerationStrategyBase)enumerator [, (EnumerationParams)params]) -> None :

    C++ signature :void __init__(_object*,RDKit::ChemicalReaction,boost::python::list,RDKit::EnumerationStrategyBase [,RDKit::EnumerationParams])

    __init__( (object)arg1, (ChemicalReaction)rxn, (tuple)reagents, (EnumerationStrategyBase)enumerator [, (EnumerationParams)params]) -> None :

    C++ signature :void __init__(_object*,RDKit::ChemicalReaction,boost::python::tuple,RDKit::EnumerationStrategyBase [,RDKit::EnumerationParams])
    """

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def GetReagents(self, arg1: EnumerateLibrary) -> VectMolVect:
        """
        Return the reagents used in this library.

        C++ signature :std::vector<std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > >, std::allocator<std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > > > > GetReagents(RDKit::EnumerateLibraryWrap {lvalue})
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class EnumerateLibraryBase(Boost.Python.instance):
    """
    Raises an exception
    This class cannot be instantiated from Python"""

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def GetEnumerator(self, arg1: EnumerateLibraryBase) -> EnumerationStrategyBase:
        """
        Returns the enumation strategy for the current library

        C++ signature :RDKit::EnumerationStrategyBase GetEnumerator(RDKit::EnumerateLibraryBase {lvalue})
        """
        ...
    def GetPosition(self, arg1: EnumerateLibraryBase) -> VectSizeT:
        """
        Returns the current enumeration position into the reagent vectors

        C++ signature :std::vector<unsigned long, std::allocator<unsigned long> > GetPosition(RDKit::EnumerateLibraryBase {lvalue})
        """
        ...
    def GetReaction(self, arg1: EnumerateLibraryBase) -> ChemicalReaction:
        """
        Returns the chemical reaction for this library

        C++ signature :RDKit::ChemicalReaction GetReaction(RDKit::EnumerateLibraryBase {lvalue})
        """
        ...
    def GetState(self, arg1: EnumerateLibraryBase) -> str:
        """
        Returns the current enumeration state (position) of the library.
        This position can be used to restart the library from a known position

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetState(RDKit::EnumerateLibraryBase {lvalue})
        """
        ...
    def InitFromString(self, arg1: EnumerateLibraryBase, data: str) -> None:
        """
        Inititialize the library from a binary string

        C++ signature :void InitFromString(RDKit::EnumerateLibraryBase {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    def ResetState(self, arg1: EnumerateLibraryBase) -> None:
        """
        Returns the current enumeration state (position) of the library to the start.

        C++ signature :void ResetState(RDKit::EnumerateLibraryBase {lvalue})"""
        ...
    def Serialize(self, arg1: EnumerateLibraryBase) -> object:
        """
        Serialize the library to a binary string.
        Note that the position in the library is serialized as well.  Care should
        be taken when serializing.  See GetState/SetState for position manipulation.

        C++ signature :boost::python::api::object Serialize(RDKit::EnumerateLibraryBase)
        """
        ...
    def SetState(self, arg1: EnumerateLibraryBase, state: str) -> None:
        """
        Sets the enumeration state (position) of the library.

        C++ signature :void SetState(RDKit::EnumerateLibraryBase {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    def next(self, arg1: EnumerateLibraryBase) -> object:
        """
        Return the next molecule from the enumeration.

        C++ signature :_object* next(RDKit::EnumerateLibraryBase*)"""
        ...
    def nextSmiles(self, arg1: EnumerateLibraryBase) -> VectorOfStringVectors:
        """
        Return the next smiles string from the enumeration.

        C++ signature :std::vector<std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >, std::allocator<std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > > > nextSmiles(RDKit::EnumerateLibraryBase {lvalue})
        """
        ...
    @classmethod
    def __bool__(cls, RDKit) -> Any: ...
    @classmethod
    def __iter__(cls, boost) -> Any: ...
    @classmethod
    def __next__(cls, RDKit) -> Any: ...
    @classmethod
    def __nonzero__(cls, RDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class EnumerationParams(Boost.Python.instance):
    """
    Controls some aspects of how the enumeration is performed.
    Options:

    reagentMaxMatchCount [ default Infinite ]This specifies how many times the reactant template can match a reagent.

    sanePartialProducts [default false]
    If true, forces all products of the reagent plus the product templatespass chemical sanitization.  Note that if the product template itself
    does not pass sanitization, then none of the products will.

    __init__( (object)arg1) -> None :

    C++ signature :void __init__(_object*)

    property reagentMaxMatchCount¶

    property sanePartialProducts¶"""

    ...
    __instance_size__: ClassVar[int] = ...
    reagentMaxMatchCount: Any
    sanePartialProducts: Any

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class EnumerationStrategyBase(Boost.Python.instance):
    """
    Raises an exception
    This class cannot be instantiated from Python"""

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def GetNumPermutations(self, arg1: EnumerationStrategyBase) -> int:
        """
        Returns the total number of results for this enumeration strategy.
        Note that some strategies are effectively infinite.

        C++ signature :unsigned long GetNumPermutations(RDKit::EnumerationStrategyBase {lvalue})
        """
        ...
    def GetPosition(self, arg1: EnumerationStrategyBase) -> VectSizeT:
        """
        Return the current indices into the arrays of reagents

        C++ signature :std::vector<unsigned long, std::allocator<unsigned long> > GetPosition(RDKit::EnumerationStrategyBase {lvalue})
        """
        ...
    def Initialize(
        self, arg1: EnumerationStrategyBase, arg2: ChemicalReaction, arg3: list
    ) -> None:
        """
        C++ signature :void Initialize(RDKit::EnumerationStrategyBase {lvalue},RDKit::ChemicalReaction {lvalue},boost::python::list)
        """
        ...
    def Skip(self, arg1: EnumerationStrategyBase, skipCount: int) -> bool:
        """
        Skip the next Nth results. note: this may be an expensive operation
        depending on the enumeration strategy used. It is recommended to use
        the enumerator state to advance to a known position

        C++ signature :bool Skip(RDKit::EnumerationStrategyBase {lvalue},unsigned long)
        """
        ...
    def Type(self, arg1: EnumerationStrategyBase) -> str:
        """
        Returns the enumeration strategy type as a string.

        C++ signature :char const* Type(RDKit::EnumerationStrategyBase {lvalue})"""
        ...
    @overload
    def next(self, arg1: EnumerationStrategyBase) -> VectSizeT:
        """


        C++ signature :void next(RDKit::EnumerationStrategyBase* {lvalue})"""
        ...
    @overload
    def next(self, arg1: EnumerationStrategyBase) -> None: ...
    @classmethod
    def __bool__(cls, RDKit) -> Any: ...
    @overload
    @classmethod
    def __copy__(cls, RDKit) -> Any: ...
    @overload
    @classmethod
    def __copy__(cls, RDKit) -> Any: ...
    @overload
    @classmethod
    def __next__(cls, RDKit) -> Any: ...
    @overload
    @classmethod
    def __next__(cls, RDKit) -> Any: ...
    @classmethod
    def __nonzero__(cls, RDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class EvenSamplePairsStrategy(EnumerationStrategyBase):
    """
    Randomly sample Pairs evenly from a collection of building blocks
    This is a good strategy for choosing a relatively small selection
    of building blocks from a larger set.  As the amount of work needed
    to retrieve the next evenly sample building block grows with the
    number of samples, this method performs progressively worse as the
    number of samples gets larger.
    See EnumerationStrategyBase for more details.

    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def Stats(self, arg1: EvenSamplePairsStrategy) -> str:
        """
        Return the statistics log of the pairs used in the current enumeration.

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > Stats(RDKit::EvenSamplePairsStrategy {lvalue})
        """
        ...
    @classmethod
    def __copy__(cls, RDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class FingerprintType(Boost.Python.enum):
    AtomPairFP: FingerprintType = ...
    MorganFP: FingerprintType = ...
    PatternFP: FingerprintType = ...
    RDKitFP: FingerprintType = ...
    TopologicalTorsion: FingerprintType = ...
    names: dict[str, FingerprintType] = ...
    values: dict[int, FingerprintType] = ...
    __slots__: ClassVar[tuple] = ...

class MOL_SPTR_VECT(Boost.Python.instance):
    """
    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def append(self, arg1: MOL_SPTR_VECT, arg2: AtomPairsParameters) -> None:
        """
        C++ signature :void append(std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > > {lvalue},boost::python::api::object)
        """
        ...
    def extend(self, arg1: MOL_SPTR_VECT, arg2: AtomPairsParameters) -> None:
        """
        C++ signature :void extend(std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > > {lvalue},boost::python::api::object)
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

class RandomSampleAllBBsStrategy(EnumerationStrategyBase):
    """
    RandomSampleAllBBsStrategy randomly samples from the reagent sets
    with the constraint that all building blocks are samples as early as possible.
    Note that this strategy never halts and can produce duplicates.

    C++ signature :void __init__(_object*)"""

    ...
    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __copy__(cls, RDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class RandomSampleStrategy(EnumerationStrategyBase):
    """
    RandomSampleStrategy simply randomly samples from the reagent sets.
    Note that this strategy never halts and can produce duplicates.

    C++ signature :void __init__(_object*)"""

    ...
    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __copy__(cls, RDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class ReactionFingerprintParams(Boost.Python.instance):
    """
    A class for storing parameters to manipulate the calculation of fingerprints of chemical reactions.
    Constructor, takes no arguments

    C++ signature :void __init__(_object*)

    __init__( (object)arg1, (bool)arg2, (float)arg3, (int)arg4, (int)arg5, (int)arg6, (FingerprintType)arg7) -> None :

    C++ signature :void __init__(_object*,bool,double,unsigned int,int,unsigned int,RDKit::FingerprintType)

    property agentWeight¶

    property bitRatioAgents¶

    property fpSize¶

    property fpType¶

    property includeAgents¶

    property nonAgentWeight¶"""

    ...
    __instance_size__: ClassVar[int] = ...
    agentWeight: Any
    bitRatioAgents: Any
    fpSize: Any
    fpType: Any
    includeAgents: Any
    nonAgentWeight: Any

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class SanitizeFlags(Boost.Python.enum):
    SANITIZE_ADJUST_REACTANTS: SanitizeFlags = ...
    SANITIZE_ALL: SanitizeFlags = ...
    SANITIZE_ATOM_MAPS: SanitizeFlags = ...
    SANITIZE_MERGEHS: SanitizeFlags = ...
    SANITIZE_NONE: SanitizeFlags = ...
    SANITIZE_RGROUP_NAMES: SanitizeFlags = ...
    names: dict[str, SanitizeFlags] = ...
    values: dict[int, SanitizeFlags] = ...
    __slots__: ClassVar[tuple] = ...

class VectMolVect(Boost.Python.instance):
    """
    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def append(self, arg1: VectMolVect, arg2: AtomPairsParameters) -> None:
        """
        C++ signature :void append(std::vector<std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > >, std::allocator<std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > > > > {lvalue},boost::python::api::object)
        """
        ...
    def extend(self, arg1: VectMolVect, arg2: AtomPairsParameters) -> None:
        """
        C++ signature :void extend(std::vector<std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > >, std::allocator<std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > > > > {lvalue},boost::python::api::object)
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

class VectSizeT(Boost.Python.instance):
    """
    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def append(self, arg1: VectSizeT, arg2: AtomPairsParameters) -> None:
        """
        C++ signature :void append(std::vector<unsigned long, std::allocator<unsigned long> > {lvalue},boost::python::api::object)
        """
        ...
    def extend(self, arg1: VectSizeT, arg2: AtomPairsParameters) -> None:
        """
        C++ signature :void extend(std::vector<unsigned long, std::allocator<unsigned long> > {lvalue},boost::python::api::object)
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

class VectorOfStringVectors(Boost.Python.instance):
    """
    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def append(self, arg1: VectorOfStringVectors, arg2: AtomPairsParameters) -> None:
        """
        C++ signature :void append(std::vector<std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >, std::allocator<std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > > > {lvalue},boost::python::api::object)
        """
        ...
    def extend(self, arg1: VectorOfStringVectors, arg2: AtomPairsParameters) -> None:
        """
        C++ signature :void extend(std::vector<std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >, std::allocator<std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > > > {lvalue},boost::python::api::object)
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

def Compute2DCoordsForReaction(
    self,
    reaction: ChemicalReaction,
    spacing: float = 2.0,
    updateProps: bool = True,
    canonOrient: bool = True,
    nFlipsPerSample: int = 0,
    nSample: int = 0,
    sampleSeed: int = 0,
    permuteDeg4Nodes: bool = False,
    bondLength: float = -1.0,
) -> None:
    """
    Compute 2D coordinates for a reaction.
    ARGUMENTS:
    reaction - the reaction of interest
    spacing - the amount of space left between components of the reaction
    canonOrient - orient the reactants and products in a canonical way

    updateProps - if set, properties such as conjugation andhybridization will be calculated for the reactant and product
    templates before generating coordinates. This should result in
    better depictions, but can lead to errors in some cases.

    nFlipsPerSample - number of rotatable bonds that areflipped at random at a time.

    nSample - Number of random samplings of rotatable bonds.
    sampleSeed - seed for the random sampling process.

    permuteDeg4Nodes - allow permutation of bonds at a degree 4node during the sampling process

    bondLength - change the default bond length for depiction

    C++ signature :void Compute2DCoordsForReaction(RDKit::ChemicalReaction {lvalue} [,double=2.0 [,bool=True [,bool=True [,unsigned int=0 [,unsigned int=0 [,int=0 [,bool=False [,double=-1.0]]]]]]]])
    """
    ...

def CreateDifferenceFingerprintForReaction(
    self,
    reaction: ChemicalReaction,
    ReactionFingerPrintParams: ReactionFingerprintParams = ...,
) -> UIntSparseIntVect:
    """
    construct a difference fingerprint for a ChemicalReaction by subtracting the reactant fingerprint from the product fingerprint

    C++ signature :RDKit::SparseIntVect<unsigned int>* CreateDifferenceFingerprintForReaction(RDKit::ChemicalReaction [,RDKit::ReactionFingerprintParams=<rdkit.Chem.rdChemReactions.ReactionFingerprintParams object at 0x7fd481545940>])
    """
    ...

def CreateStructuralFingerprintForReaction(
    self,
    reaction: ChemicalReaction,
    ReactionFingerPrintParams: ReactionFingerprintParams = ...,
) -> ExplicitBitVect:
    """
    construct a structural fingerprint for a ChemicalReaction by concatenating the reactant fingerprint and the product fingerprint

    C++ signature :ExplicitBitVect* CreateStructuralFingerprintForReaction(RDKit::ChemicalReaction [,RDKit::ReactionFingerprintParams=<rdkit.Chem.rdChemReactions.ReactionFingerprintParams object at 0x7fd481545a40>])
    """
    ...

def EnumerateLibraryCanSerialize(self) -> bool:
    """
    Returns True if the EnumerateLibrary is serializable (requires boost serialization

    C++ signature :bool EnumerateLibraryCanSerialize()"""
    ...

def GetChemDrawRxnAdjustParams(self) -> AdjustQueryParameters:
    """
    (deprecated, see MatchOnlyAtRgroupsAdjustParams)Returns the chemdraw style adjustment parameters for reactant templates

    C++ signature :RDKit::MolOps::AdjustQueryParameters GetChemDrawRxnAdjustParams()"""
    ...

def GetDefaultAdjustParams(self) -> AdjustQueryParameters:
    """
    Returns the default adjustment parameters for reactant templates

    C++ signature :RDKit::MolOps::AdjustQueryParameters GetDefaultAdjustParams()"""
    ...

def HasAgentTemplateSubstructMatch(
    self, reaction: ChemicalReaction, queryReaction: ChemicalReaction
) -> bool:
    """
    tests if the agents of a queryReaction are the same as those of a reaction

    C++ signature :bool HasAgentTemplateSubstructMatch(RDKit::ChemicalReaction,RDKit::ChemicalReaction)
    """
    ...

def HasProductTemplateSubstructMatch(
    self, reaction: ChemicalReaction, queryReaction: ChemicalReaction
) -> bool:
    """
    tests if the products of a queryReaction are substructures of the products of a reaction

    C++ signature :bool HasProductTemplateSubstructMatch(RDKit::ChemicalReaction,RDKit::ChemicalReaction)
    """
    ...

def HasReactantTemplateSubstructMatch(
    self, reaction: ChemicalReaction, queryReaction: ChemicalReaction
) -> bool:
    """
    tests if the reactants of a queryReaction are substructures of the reactants of a reaction

    C++ signature :bool HasReactantTemplateSubstructMatch(RDKit::ChemicalReaction,RDKit::ChemicalReaction)
    """
    ...

def HasReactionAtomMapping(self, arg1: ChemicalReaction) -> bool:
    """
    tests if a reaction obtains any atom mapping

    C++ signature :bool HasReactionAtomMapping(RDKit::ChemicalReaction)"""
    ...

def HasReactionSubstructMatch(
    self,
    reaction: ChemicalReaction,
    queryReaction: ChemicalReaction,
    includeAgents: bool = False,
) -> bool:
    """
    tests if the queryReaction is a substructure of a reaction

    C++ signature :bool HasReactionSubstructMatch(RDKit::ChemicalReaction,RDKit::ChemicalReaction [,bool=False])
    """
    ...

def IsReactionTemplateMoleculeAgent(self, molecule: Mol, agentThreshold: float) -> bool:
    """
    tests if a molecule can be classified as an agent depending on the ratio of mapped atoms and a give threshold

    C++ signature :bool IsReactionTemplateMoleculeAgent(RDKit::ROMol,double)"""
    ...

def MatchOnlyAtRgroupsAdjustParams(self) -> AdjustQueryParameters:
    """
    Only match at the specified rgroup locations in the reactant templates

    C++ signature :RDKit::MolOps::AdjustQueryParameters MatchOnlyAtRgroupsAdjustParams()
    """
    ...

def PreprocessReaction(
    self, reaction: ChemicalReaction, queries: dict = {}, propName: str = "molFileValue"
) -> object:
    """
    A function for preprocessing reactions with more specific queries.
    Queries are indicated by labels on atoms (molFileAlias property by default)
    When these labels are found, more specific queries are placed on the atoms.
    By default, the available quieries come from

    FilterCatalog.GetFlattenedFunctionalGroupHierarchy(True)n

    Sample Usage:>>> from rdkit import Chem, RDConfig
    >>> from rdkit.Chem import MolFromSmiles, AllChem
    >>> from rdkit.Chem.rdChemReactions import PreprocessReaction
    >>> import os
    >>> testFile = os.path.join(RDConfig.RDCodeDir,'Chem','SimpleEnum','test_data','boronic1.rxn')
    >>> rxn = AllChem.ReactionFromRxnFile(testFile)
    >>> rxn.Initialize()
    >>> nWarn,nError,nReacts,nProds,reactantLabels = PreprocessReaction(rxn)
    >>> nWarn
    0
    >>> nError
    0
    >>> nReacts
    2
    >>> nProds
    1
    >>> reactantLabels
    (((0, 'halogen.bromine.aromatic'),), ((1, 'boronicacid'),))

    If there are functional group labels in the input reaction (via atoms with molFileValue properties),
    the corresponding atoms will have queries added to them so that they only match such things. We can
    see this here:
    >>> rxn = AllChem.ReactionFromRxnFile(testFile)
    >>> rxn.Initialize()
    >>> r1 = rxn.GetReactantTemplate(0)
    >>> m1 = Chem.MolFromSmiles('CCBr')
    >>> m2 = Chem.MolFromSmiles('c1ccccc1Br')

    These both match because the reaction file itself just has R1-Br:>>> m1.HasSubstructMatch(r1)
    True
    >>> m2.HasSubstructMatch(r1)
    True

    After preprocessing, we only match the aromatic Br:>>> d = PreprocessReaction(rxn)
    >>> m1.HasSubstructMatch(r1)
    False
    >>> m2.HasSubstructMatch(r1)
    True

    We also support or queries in the values field (separated by commas):>>> testFile = os.path.join(RDConfig.RDCodeDir,'Chem','SimpleEnum','test_data','azide_reaction.rxn')
    >>> rxn = AllChem.ReactionFromRxnFile(testFile)
    >>> rxn.Initialize()
    >>> reactantLabels = PreprocessReaction(rxn)[-1]
    >>> reactantLabels
    (((1, 'azide'),), ((1, 'carboxylicacid,acidchloride'),))
    >>> m1 = Chem.MolFromSmiles('CC(=O)O')
    >>> m2 = Chem.MolFromSmiles('CC(=O)Cl')
    >>> m3 = Chem.MolFromSmiles('CC(=O)N')
    >>> r2 = rxn.GetReactantTemplate(1)
    >>> m1.HasSubstructMatch(r2)
    True
    >>> m2.HasSubstructMatch(r2)
    True
    >>> m3.HasSubstructMatch(r2)
    False

    unrecognized final group types are returned as None:>>> testFile = os.path.join(RDConfig.RDCodeDir,'Chem','SimpleEnum','test_data','bad_value1.rxn')
    >>> rxn = AllChem.ReactionFromRxnFile(testFile)
    >>> rxn.Initialize()
    >>> nWarn,nError,nReacts,nProds,reactantLabels = PreprocessReaction(rxn)
    Traceback (most recent call last):
      ...
    KeyError: 'boromicacid'

    One unrecognized group type in a comma-separated list makes the whole thing fail:>>> testFile = os.path.join(RDConfig.RDCodeDir,'Chem','SimpleEnum','test_data','bad_value2.rxn')
    >>> rxn = AllChem.ReactionFromRxnFile(testFile)
    >>> rxn.Initialize()
    >>> nWarn,nError,nReacts,nProds,reactantLabels = PreprocessReaction(rxn)
    Traceback (most recent call last):
      ...
    KeyError: 'carboxylicacid,acidchlroide'
    >>> testFile = os.path.join(RDConfig.RDCodeDir,'Chem','SimpleEnum','test_data','bad_value3.rxn')
    >>> rxn = AllChem.ReactionFromRxnFile(testFile)
    >>> rxn.Initialize()
    >>> nWarn,nError,nReacts,nProds,reactantLabels = PreprocessReaction(rxn)
    Traceback (most recent call last):
      ...
    KeyError: 'carboxyliccaid,acidchloride'
    >>> rxn = rdChemReactions.ChemicalReaction()
    >>> rxn.Initialize()
    >>> nWarn,nError,nReacts,nProds,reactantLabels = PreprocessReaction(rxn)
    >>> reactantLabels
    ()
    >>> reactantLabels == ()
    True

    C++ signature :boost::python::api::object PreprocessReaction(RDKit::ChemicalReaction {lvalue} [,boost::python::dict={} [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’molFileValue’]])
    """
    ...

def ReactionFromMolecule(self, arg1: Mol) -> ChemicalReaction:
    """
    construct a ChemicalReaction from an molecule if the RXN role property of the molecule is set

    C++ signature :RDKit::ChemicalReaction* ReactionFromMolecule(RDKit::ROMol)"""
    ...

def ReactionFromPNGFile(self, arg1: str) -> ChemicalReaction:
    """
    construct a ChemicalReaction from metadata in a PNG file

    C++ signature :RDKit::ChemicalReaction* ReactionFromPNGFile(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
    ...

def ReactionFromPNGString(self, arg1: str) -> ChemicalReaction:
    """
    construct a ChemicalReaction from an string with PNG data

    C++ signature :RDKit::ChemicalReaction* ReactionFromPNGString(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
    ...

def ReactionFromRxnBlock(
    self,
    rxnblock: str,
    sanitize: bool = False,
    removeHs: bool = False,
    strictParsing: bool = True,
) -> ChemicalReaction:
    """
    construct a ChemicalReaction from a string in MDL rxn format

    C++ signature :RDKit::ChemicalReaction* ReactionFromRxnBlock(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False [,bool=False [,bool=True]]])
    """
    ...

def ReactionFromRxnFile(
    self,
    filename: str,
    sanitize: bool = False,
    removeHs: bool = False,
    strictParsing: bool = True,
) -> ChemicalReaction:
    """
    construct a ChemicalReaction from an MDL rxn file

    C++ signature :RDKit::ChemicalReaction* ReactionFromRxnFile(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False [,bool=False [,bool=True]]])
    """
    ...

def ReactionFromSmarts(
    self, SMARTS: str, replacements: dict = {}, useSmiles: bool = False
) -> ChemicalReaction:
    """
    construct a ChemicalReaction from a reaction SMARTS string.
    see the documentation for rdkit.Chem.MolFromSmiles for an explanation
    of the replacements argument.

    C++ signature :RDKit::ChemicalReaction* ReactionFromSmarts(char const* [,boost::python::dict={} [,bool=False]])
    """
    ...

def ReactionMetadataToPNGFile(
    self,
    mol: ChemicalReaction,
    filename: AtomPairsParameters,
    includePkl: bool = True,
    includeSmiles: bool = True,
    includeSmarts: bool = False,
    includeMol: bool = False,
) -> object:
    """
    Reads the contents of a PNG file and adds metadata about a reaction to it. The modified file contents are returned.

    C++ signature :boost::python::api::object ReactionMetadataToPNGFile(RDKit::ChemicalReaction,boost::python::api::object [,bool=True [,bool=True [,bool=False [,bool=False]]]])
    """
    ...

def ReactionMetadataToPNGString(
    self,
    mol: ChemicalReaction,
    pngdata: AtomPairsParameters,
    includePkl: bool = True,
    includeSmiles: bool = True,
    includeSmarts: bool = False,
    includeRxn: bool = False,
) -> object:
    """
    Adds metadata about a reaction to the PNG string passed in.The modified string is returned.

    C++ signature :boost::python::api::object ReactionMetadataToPNGString(RDKit::ChemicalReaction,boost::python::api::object [,bool=True [,bool=True [,bool=False [,bool=False]]]])
    """
    ...

def ReactionToMolecule(self, reaction: ChemicalReaction) -> Mol:
    """
    construct a molecule for a ChemicalReaction with RXN role property set

    C++ signature :RDKit::ROMol* ReactionToMolecule(RDKit::ChemicalReaction)"""
    ...

def ReactionToRxnBlock(
    self,
    reaction: ChemicalReaction,
    separateAgents: bool = False,
    forceV3000: bool = False,
) -> str:
    """
    construct a string in MDL rxn format for a ChemicalReaction

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > ReactionToRxnBlock(RDKit::ChemicalReaction [,bool=False [,bool=False]])
    """
    ...

def ReactionToSmarts(self, reaction: ChemicalReaction) -> str:
    """
    construct a reaction SMARTS string for a ChemicalReaction

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > ReactionToSmarts(RDKit::ChemicalReaction)
    """
    ...

def ReactionToSmiles(self, reaction: ChemicalReaction, canonical: bool = True) -> str:
    """
    construct a reaction SMILES string for a ChemicalReaction

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > ReactionToSmiles(RDKit::ChemicalReaction [,bool=True])
    """
    ...

def ReactionToV3KRxnBlock(
    self, reaction: ChemicalReaction, separateAgents: bool = False
) -> str:
    """
    construct a string in MDL v3000 rxn format for a ChemicalReaction

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > ReactionToV3KRxnBlock(RDKit::ChemicalReaction [,bool=False])
    """
    ...

def ReduceProductToSideChains(self, product: Mol, addDummyAtoms: bool = True) -> Mol:
    """
    reduce the product of a reaction to the side chains added by the reaction.              The output is a molecule with attached wildcards indicating where the product was attached.              The dummy atom has the same reaction-map number as the product atom (if available).

    C++ signature :RDKit::ROMol* ReduceProductToSideChains(boost::shared_ptr<RDKit::ROMol> [,bool=True])
    """
    ...

def RemoveMappingNumbersFromReactions(self, reaction: ChemicalReaction) -> None:
    """
    Removes the mapping numbers from the molecules of a reaction

    C++ signature :void RemoveMappingNumbersFromReactions(RDKit::ChemicalReaction)"""
    ...

def SanitizeRxn(
    self,
    rxn: ChemicalReaction,
    sanitizeOps: int = 4294967295,
    params: rdmolops.AdjustQueryParameters = ...,
    catchErrors: bool = False,
) -> SanitizeFlags:
    """
    Does some sanitization of the reactant and product templates of a reaction.

    The reaction is modified in place.
    If sanitization fails, an exception will be thrown unless catchErrors is set

    ARGUMENTS:

    rxn: the reaction to be modified
    sanitizeOps: (optional) reaction sanitization operations to be carried out
    these should be constructed by or’ing together the
    operations in rdkit.Chem.rdChemReactions.SanitizeFlags
    optional adjustment parameters for changing the meaning of the substructure
    matching done in the templates.  The default is
    rdkit.Chem.rdChemReactions.DefaultRxnAdjustParams which aromatizes
    kekule structures if possible.
    catchErrors: (optional) if provided, instead of raising an exception
    when sanitization fails (the default behavior), the
    first operation that failed (as defined in rdkit.Chem.rdChemReactions.SanitizeFlags)
    is returned. Zero is returned on success.

    The operations carried out by default are:
    fixRGroups(): sets R group labels on mapped dummy atoms when possible
    fixAtomMaps(): attempts to set atom maps on unmapped R groups
    adjustTemplate(): calls adjustQueryProperties() on all reactant templates
    fixHs(): merges explicit Hs in the reactant templates that don’t map to heavy atoms

    C++ signature :RDKit::RxnOps::SanitizeRxnFlags SanitizeRxn(RDKit::ChemicalReaction {lvalue} [,unsigned long=4294967295 [,RDKit::MolOps::AdjustQueryParameters=<rdkit.Chem.rdmolops.AdjustQueryParameters object at 0x7fd4815545e0> [,bool=False]]])
    """
    ...

def UpdateProductsStereochemistry(self, reaction: ChemicalReaction) -> None:
    """
    Caution: This is an expert-user function which will change a property (molInversionFlag) of your products.          This function is called by default using the RXN or SMARTS parser for reactions and should really only be called if reactions have been constructed some other way.          The function updates the stereochemistry of the product by considering 4 different cases: inversion, retention, removal, and introduction

    C++ signature :void UpdateProductsStereochemistry(RDKit::ChemicalReaction*)"""
    ...
