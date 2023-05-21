"""
rdkit.Chem.rdfiltercatalog module¶
"""
from typing import Any, ClassVar, overload

import Boost.Python
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdfiltercatalog import FilterCatalogs
from rdkit.Chem.rdMolDescriptors import AtomPairsParameters
from rdkit.rdBase import _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE

class ExclusionList(FilterMatcherBase):
    """
    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def AddPattern(self, arg1: ExclusionList, arg2: FilterMatcherBase) -> None:
        """
        Add a FilterMatcherBase that should not appear in a molecule

        C++ signature :void AddPattern(RDKit::ExclusionList {lvalue},RDKit::FilterMatcherBase)
        """
        ...
    def SetExclusionPatterns(
        self, arg1: ExclusionList, arg2: AtomPairsParameters
    ) -> None:
        """
        Set a list of FilterMatcherBases that should not appear in a molecule

        C++ signature :void SetExclusionPatterns(RDKit::ExclusionList {lvalue},boost::python::api::object)
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class FilterCatalog(Boost.Python.instance):
    """
    C++ signature :void __init__(_object*)

    __init__( (object)arg1, (str)arg2) -> None :

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)

    __init__( (object)arg1, (FilterCatalogParams)arg2) -> None :

    C++ signature :void __init__(_object*,RDKit::FilterCatalogParams)

    __init__( (object)arg1, (FilterCatalogs)arg2) -> None :

    C++ signature :void __init__(_object*,RDKit::FilterCatalogParams::FilterCatalogs)"""

    __getstate_manages_dict__: ClassVar[bool] = ...
    __instance_size__: ClassVar[int] = ...
    __safe_for_unpickling__: ClassVar[bool] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def AddEntry(
        self, entry: FilterCatalog, updateFPLength: FilterCatalogEntry = False
    ) -> None:
        """
        Add a FilterCatalogEntry to the catalog

        C++ signature :void AddEntry(RDKit::FilterCatalog {lvalue} [,RDKit::FilterCatalogEntry*=False])
        """
        ...
    def GetEntry(self, arg1: FilterCatalog, idx: int) -> FilterCatalogEntry:
        """
        Return the FilterCatalogEntry at the specified index

        C++ signature :boost::shared_ptr<RDKit::FilterCatalogEntry const> GetEntry(RDKit::FilterCatalog {lvalue},unsigned int)
        """
        ...
    def GetEntryWithIdx(self, arg1: FilterCatalog, idx: int) -> FilterCatalogEntry:
        """
        Return the FilterCatalogEntry at the specified index

        C++ signature :boost::shared_ptr<RDKit::FilterCatalogEntry const> GetEntryWithIdx(RDKit::FilterCatalog {lvalue},unsigned int)
        """
        ...
    def GetFilterMatches(self, arg1: FilterCatalog, mol: Mol) -> VectFilterMatch:
        """
        Return every matching filter from all catalog entries that match mol

        C++ signature :std::vector<RDKit::FilterMatch, std::allocator<RDKit::FilterMatch> > GetFilterMatches(RDKit::FilterCatalog {lvalue},RDKit::ROMol)
        """
        ...
    def GetFirstMatch(self, arg1: FilterCatalog, mol: Mol) -> FilterCatalogEntry:
        """
        Return the first catalog entry that matches mol

        C++ signature :boost::shared_ptr<RDKit::FilterCatalogEntry const> GetFirstMatch(RDKit::FilterCatalog {lvalue},RDKit::ROMol)
        """
        ...
    def GetMatches(self, arg1: FilterCatalog, mol: Mol) -> FilterCatalogEntryList:
        """
        Return all catalog entries that match mol

        C++ signature :std::vector<boost::shared_ptr<RDKit::FilterCatalogEntry const>, std::allocator<boost::shared_ptr<RDKit::FilterCatalogEntry const> > > GetMatches(RDKit::FilterCatalog {lvalue},RDKit::ROMol)
        """
        ...
    def GetNumEntries(self, arg1: FilterCatalog) -> int:
        """
        Returns the number of entries in the catalog

        C++ signature :unsigned int GetNumEntries(RDKit::FilterCatalog {lvalue})"""
        ...
    def HasMatch(self, arg1: FilterCatalog, mol: Mol) -> bool:
        """
        Returns True if the catalog has an entry that matches mol

        C++ signature :bool HasMatch(RDKit::FilterCatalog {lvalue},RDKit::ROMol)"""
        ...
    def RemoveEntry(self, arg1: FilterCatalog, arg2: AtomPairsParameters) -> bool:
        """
        Remove the given entry from the catalog

        C++ signature :bool RemoveEntry(RDKit::FilterCatalog {lvalue},boost::python::api::object)
        """
        ...
    def Serialize(self, arg1: FilterCatalog) -> object:
        """
        C++ signature :boost::python::api::object Serialize(RDKit::FilterCatalog)"""
        ...
    @classmethod
    def __getinitargs__(cls, RDKit) -> Any: ...
    @classmethod
    def __getstate__(cls, boost) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @classmethod
    def __setstate__(cls, state) -> Any: ...

class FilterCatalogEntry(Boost.Python.instance):
    """
    A filter catalog entry is an entry in a filter catalog.
    Each filter is named and is used to flag a molecule usually for some
    undesirable property.
    For example, a PAINS (Pan Assay INterference) catalog entry be appear as
    follows:
    >>> from rdkit.Chem.FilterCatalog import *
    >>> params = FilterCatalogParams()
    >>> params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
    True
    >>> catalog = FilterCatalog(params)
    >>> mol = Chem.MolFromSmiles('O=C(Cn1cnc2c1c(=O)n(C)c(=O)n2C)N/N=C/c1c(O)ccc2c1cccc2')
    >>> entry = catalog.GetFirstMatch(mol)
    >>> print (entry.GetProp('Scope'))
    PAINS filters (family A)
    >>> print (entry.GetDescription())
    hzone_phenol_A(479)

    __init__( (object)arg1) -> None :

    C++ signature :void __init__(_object*)

    __init__( (object)arg1, (str)arg2, (FilterMatcherBase)arg3) -> None :

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,RDKit::FilterMatcherBase {lvalue})
    """

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def ClearProp(self, arg1: FilterCatalogEntry, arg2: str) -> None:
        """
        C++ signature :void ClearProp(RDKit::FilterCatalogEntry {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    def GetDescription(self, arg1: FilterCatalogEntry) -> str:
        """
        Get the description of the catalog entry

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetDescription(RDKit::FilterCatalogEntry {lvalue})
        """
        ...
    def GetFilterMatches(self, arg1: FilterCatalogEntry, mol: Mol) -> VectFilterMatch:
        """
        Retrieve the list of filters that match the molecule

        C++ signature :std::vector<RDKit::FilterMatch, std::allocator<RDKit::FilterMatch> > GetFilterMatches(RDKit::FilterCatalogEntry {lvalue},RDKit::ROMol)
        """
        ...
    def GetProp(self, arg1: FilterCatalogEntry, arg2: str) -> str:
        """
        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetProp(RDKit::FilterCatalogEntry {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    def GetPropList(
        self, arg1: FilterCatalogEntry
    ) -> _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE:
        """
        C++ signature :std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > GetPropList(RDKit::FilterCatalogEntry {lvalue})
        """
        ...
    def HasFilterMatch(self, arg1: FilterCatalogEntry, mol: Mol) -> bool:
        """
        Returns True if the catalog entry contains filters that match the molecule

        C++ signature :bool HasFilterMatch(RDKit::FilterCatalogEntry {lvalue},RDKit::ROMol)
        """
        ...
    def IsValid(self, arg1: FilterCatalogEntry) -> bool:
        """
        C++ signature :bool IsValid(RDKit::FilterCatalogEntry {lvalue})"""
        ...
    def Serialize(self, arg1: FilterCatalogEntry) -> object:
        """
        C++ signature :boost::python::api::object Serialize(RDKit::FilterCatalogEntry)
        """
        ...
    def SetDescription(self, arg1: FilterCatalogEntry, description: str) -> None:
        """
        Set the description of the catalog entry

        C++ signature :void SetDescription(RDKit::FilterCatalogEntry {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    def SetProp(self, arg1: FilterCatalogEntry, arg2: str, arg3: str) -> None:
        """
        C++ signature :void SetProp(RDKit::FilterCatalogEntry {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class FilterCatalogEntryList(Boost.Python.instance):
    """
    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def append(self, arg1: FilterCatalogEntryList, arg2: AtomPairsParameters) -> None:
        """
        C++ signature :void append(std::vector<boost::shared_ptr<RDKit::FilterCatalogEntry const>, std::allocator<boost::shared_ptr<RDKit::FilterCatalogEntry const> > > {lvalue},boost::python::api::object)
        """
        ...
    def extend(self, arg1: FilterCatalogEntryList, arg2: AtomPairsParameters) -> None:
        """
        C++ signature :void extend(std::vector<boost::shared_ptr<RDKit::FilterCatalogEntry const>, std::allocator<boost::shared_ptr<RDKit::FilterCatalogEntry const> > > {lvalue},boost::python::api::object)
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

class FilterCatalogListOfEntryList(Boost.Python.instance):
    """
    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def append(
        self, arg1: FilterCatalogListOfEntryList, arg2: AtomPairsParameters
    ) -> None:
        """
        C++ signature :void append(std::vector<std::vector<boost::shared_ptr<RDKit::FilterCatalogEntry const>, std::allocator<boost::shared_ptr<RDKit::FilterCatalogEntry const> > >, std::allocator<std::vector<boost::shared_ptr<RDKit::FilterCatalogEntry const>, std::allocator<boost::shared_ptr<RDKit::FilterCatalogEntry const> > > > > {lvalue},boost::python::api::object)
        """
        ...
    def extend(
        self, arg1: FilterCatalogListOfEntryList, arg2: AtomPairsParameters
    ) -> None:
        """
        C++ signature :void extend(std::vector<std::vector<boost::shared_ptr<RDKit::FilterCatalogEntry const>, std::allocator<boost::shared_ptr<RDKit::FilterCatalogEntry const> > >, std::allocator<std::vector<boost::shared_ptr<RDKit::FilterCatalogEntry const>, std::allocator<boost::shared_ptr<RDKit::FilterCatalogEntry const> > > > > {lvalue},boost::python::api::object)
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

class FilterCatalogParams(Boost.Python.instance):
    """
    C++ signature :void __init__(_object*)

    __init__( (object)arg1, (FilterCatalogs)arg2) -> None :Construct from a FilterCatalogs identifier (i.e. FilterCatalogParams.PAINS)

    C++ signature :void __init__(_object*,RDKit::FilterCatalogParams::FilterCatalogs)"""

    class FilterCatalogs(Boost.Python.enum):
        ALL: FilterCatalogs = ...
        BRENK: FilterCatalogs = ...
        CHEMBL: FilterCatalogs = ...
        CHEMBL_BMS: FilterCatalogs = ...
        CHEMBL_Dundee: FilterCatalogs = ...
        CHEMBL_Glaxo: FilterCatalogs = ...
        CHEMBL_Inpharmatica: FilterCatalogs = ...
        CHEMBL_LINT: FilterCatalogs = ...
        CHEMBL_MLSMR: FilterCatalogs = ...
        CHEMBL_SureChEMBL: FilterCatalogs = ...
        NIH: FilterCatalogs = ...
        PAINS: FilterCatalogs = ...
        PAINS_A: FilterCatalogs = ...
        PAINS_B: FilterCatalogs = ...
        PAINS_C: FilterCatalogs = ...
        ZINC: FilterCatalogs = ...
        names: dict[str, FilterCatalogs] = ...
        values: dict[int, FilterCatalogs] = ...
        __slots__: ClassVar[tuple] = ...
    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def AddCatalog(
        self, arg1: FilterCatalogParams, arg2: FilterCatalogParams.FilterCatalogs
    ) -> bool:
        """
        C++ signature :bool AddCatalog(RDKit::FilterCatalogParams {lvalue},RDKit::FilterCatalogParams::FilterCatalogs)
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class FilterHierarchyMatcher(FilterMatcherBase):
    """
    Hierarchical Filter
    basic constructors: FilterHierarchyMatcher( matcher )
    where can be any FilterMatcherBase (SmartsMatcher, etc)

    FilterHierarchyMatcher’s have children and can form matchingtrees.  then GetFilterMatches is called, the most specific (
    i.e. lowest node in a branch) is returned.

    n.b. A FilterHierarchicalMatcher of functional groups is returnedby calling GetFunctionalGroupHierarchy()

    >>> from rdkit.Chem import MolFromSmiles
    >>> from rdkit.Chem.FilterCatalog import *
    >>> functionalGroups = GetFunctionalGroupHierarchy()
    >>> [match.filterMatch.GetName()
    ...     for match in functionalGroups.GetFilterMatches(
    ...         MolFromSmiles('c1ccccc1Cl'))]
    ['Halogen.Aromatic', 'Halogen.NotFluorine.Aromatic']

    C++ signature :void __init__(_object*)

    __init__( (object)arg1, (FilterMatcherBase)arg2) -> None :Construct from a filtermatcher

    C++ signature :void __init__(_object*,RDKit::FilterMatcherBase)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def AddChild(
        self, arg1: FilterHierarchyMatcher, arg2: FilterHierarchyMatcher
    ) -> FilterHierarchyMatcher:
        """
        Add a child node to this hierarchy.

        C++ signature :boost::shared_ptr<RDKit::FilterHierarchyMatcher> AddChild(RDKit::FilterHierarchyMatcher {lvalue},RDKit::FilterHierarchyMatcher)
        """
        ...
    def SetPattern(self, arg1: FilterHierarchyMatcher, arg2: FilterMatcherBase) -> None:
        """
        Set the filtermatcher pattern for this node.  An empty node is considered a root node and passes along the matches to the children.

        C++ signature :void SetPattern(RDKit::FilterHierarchyMatcher {lvalue},RDKit::FilterMatcherBase)
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class FilterMatch(Boost.Python.instance):
    """
    Object that holds the result of running FilterMatcherBase::GetMatches

    filterMatch holds the FilterMatchBase that triggered the match
    atomParis holds the [ (query_atom_idx, target_atom_idx) ] pairs for the matches.

    Note that some matches may not have atom pairs (especially matches that use FilterMatchOps.Not

    C++ signature :void __init__(_object*,boost::shared_ptr<RDKit::FilterMatcherBase>,std::vector<std::pair<int, int>, std::allocator<std::pair<int, int> > >)

    property atomPairs¶

    property filterMatch¶"""

    ...
    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @property
    def atomPairs(self) -> Any: ...
    @property
    def filterMatch(self) -> Any: ...

class FilterMatcherBase(Boost.Python.instance):
    """
    Base class for matching molecules to filters.

    A FilterMatcherBase supplies the following API
    - IsValid() returns True if the matcher is valid for use, False otherwise
    - HasMatch(mol) returns True if the molecule matches the filter
    - GetMatches(mol) -> [FilterMatch, FilterMatch] returns all the FilterMatch data

    that matches the molecule

    print( FilterMatcherBase ) will print user-friendly information about the filter
    Note that a FilterMatcherBase can be combined from may FilterMatcherBases
    This is why GetMatches can return multiple FilterMatcherBases.
    >>> from rdkit.Chem.FilterCatalog import *
    >>> carbon_matcher = SmartsMatcher(‘Carbon’, ‘[#6]’, 0, 1)
    >>> oxygen_matcher = SmartsMatcher(‘Oxygen’, ‘[#8]’, 0, 1)
    >>> co_matcher = FilterMatchOps.Or(carbon_matcher, oxygen_matcher)
    >>> mol = Chem.MolFromSmiles(‘C’)
    >>> matches = co_matcher.GetMatches(mol)
    >>> len(matches)
    1
    >>> print(matches[0].filterMatch)
    Carbon
    Raises an exception
    This class cannot be instantiated from Python"""

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def GetMatches(self, arg1: FilterMatcherBase, mol: Mol) -> VectFilterMatch:
        """
        Returns the list of matching subfilters mol matches any filter

        C++ signature :std::vector<RDKit::FilterMatch, std::allocator<RDKit::FilterMatch> > GetMatches(RDKit::FilterMatcherBase {lvalue},RDKit::ROMol)
        """
        ...
    def GetName(self, arg1: FilterMatcherBase) -> str:
        """
        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetName(RDKit::FilterMatcherBase {lvalue})
        """
        ...
    def HasMatch(self, arg1: FilterMatcherBase, mol: Mol) -> bool:
        """
        Returns True if mol matches the filter

        C++ signature :bool HasMatch(RDKit::FilterMatcherBase {lvalue},RDKit::ROMol)"""
        ...
    def IsValid(self, arg1: FilterMatcherBase) -> bool:
        """
        Return True if the filter matcher is valid, False otherwise

        C++ signature :bool IsValid(RDKit::FilterMatcherBase {lvalue})"""
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class IntPair(Boost.Python.instance):
    """
    C++ signature :void __init__(_object*)

    __init__( (object)arg1, (int)arg2, (int)arg3) -> None :

    C++ signature :void __init__(_object*,int,int)

    property query¶

    property target¶"""

    ...
    __instance_size__: ClassVar[int] = ...
    query: Any
    target: Any

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __getitem__(cls, index) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class MatchTypeVect(Boost.Python.instance):
    """
    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def append(self, arg1: MatchTypeVect, arg2: AtomPairsParameters) -> None:
        """
        C++ signature :void append(std::vector<std::pair<int, int>, std::allocator<std::pair<int, int> > > {lvalue},boost::python::api::object)
        """
        ...
    def extend(self, arg1: MatchTypeVect, arg2: AtomPairsParameters) -> None:
        """
        C++ signature :void extend(std::vector<std::pair<int, int>, std::allocator<std::pair<int, int> > > {lvalue},boost::python::api::object)
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

class MolList(Boost.Python.instance):
    """
    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def append(self, arg1: MolList, arg2: AtomPairsParameters) -> None:
        """
        C++ signature :void append(std::vector<RDKit::ROMol*, std::allocator<RDKit::ROMol*> > {lvalue},boost::python::api::object)
        """
        ...
    def extend(self, arg1: MolList, arg2: AtomPairsParameters) -> None:
        """
        C++ signature :void extend(std::vector<RDKit::ROMol*, std::allocator<RDKit::ROMol*> > {lvalue},boost::python::api::object)
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

class PythonFilterMatcher(FilterMatcherBase):
    """
    C++ signature :void __init__(_object*,_object*)"""

    ...
    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class SmartsMatcher(FilterMatcherBase):
    """
    Smarts Matcher Filter
    basic constructors:
    SmartsMatcher( name, smarts_pattern, minCount=1, maxCount=UINT_MAX )
    SmartsMatcher( name, molecule, minCount=1, maxCount=UINT_MAX )

    note: If the supplied smarts pattern is not valid, the IsValid() function willreturn False

    >>> from rdkit.Chem.FilterCatalog import *
    >>> minCount, maxCount = 1,2
    >>> carbon_matcher = SmartsMatcher('Carbon', '[#6]', minCount, maxCount)
    >>> print (carbon_matcher.HasMatch(Chem.MolFromSmiles('CC')))
    True
    >>> print (carbon_matcher.HasMatch(Chem.MolFromSmiles('CCC')))
    False
    >>> carbon_matcher.SetMinCount(2)
    >>> print (carbon_matcher.HasMatch(Chem.MolFromSmiles('C')))
    False
    >>> carbon_matcher.SetMaxCount(3)
    >>> print (carbon_matcher.HasMatch(Chem.MolFromSmiles('CCC')))
    True

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)

    __init__( (object)arg1, (Mol)arg2) -> None :Construct from a molecule

    C++ signature :void __init__(_object*,RDKit::ROMol)

    __init__( (object)arg1, (str)name, (Mol)mol [, (int)minCount=1 [, (int)maxCount=4294967295]]) -> None :Construct from a name, molecule, minimum and maximum count

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,RDKit::ROMol [,unsigned int=1 [,unsigned int=4294967295]])

    __init__( (object)arg1, (str)name, (str)smarts [, (int)minCount=1 [, (int)maxCount=4294967295]]) -> None :Construct from a name,smarts pattern, minimum and maximum count

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,unsigned int=1 [,unsigned int=4294967295]])
    """

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def GetMaxCount(self, arg1: SmartsMatcher) -> int:
        """
        Get the maximum times pattern can appear for the filter to match

        C++ signature :unsigned int GetMaxCount(RDKit::SmartsMatcher {lvalue})"""
        ...
    def GetMinCount(self, arg1: SmartsMatcher) -> int:
        """
        Get the minimum times pattern must appear for the filter to match

        C++ signature :unsigned int GetMinCount(RDKit::SmartsMatcher {lvalue})"""
        ...
    def GetPattern(self, arg1: SmartsMatcher) -> Mol:
        """
        C++ signature :boost::shared_ptr<RDKit::ROMol> GetPattern(RDKit::SmartsMatcher {lvalue})
        """
        ...
    def IsValid(self, arg1: SmartsMatcher) -> bool:
        """
        Returns True if the SmartsMatcher is valid

        C++ signature :bool IsValid(RDKit::SmartsMatcher {lvalue})"""
        ...
    def SetMaxCount(self, arg1: SmartsMatcher, count: int) -> None:
        """
        Set the maximum times pattern can appear for the filter to match

        C++ signature :void SetMaxCount(RDKit::SmartsMatcher {lvalue},unsigned int)"""
        ...
    def SetMinCount(self, arg1: SmartsMatcher, count: int) -> None:
        """
        Set the minimum times pattern must appear to match

        C++ signature :void SetMinCount(RDKit::SmartsMatcher {lvalue},unsigned int)"""
        ...
    @overload
    def SetPattern(self, arg1: SmartsMatcher, arg2: Mol) -> None:
        """
        Set the smarts pattern for the Smarts Matcher (warning: MinimumCount is not reset)

        C++ signature :void SetPattern(RDKit::SmartsMatcher {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    @overload
    def SetPattern(self, arg1: SmartsMatcher, arg2: str) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class VectFilterMatch(Boost.Python.instance):
    """
    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def append(self, arg1: VectFilterMatch, arg2: AtomPairsParameters) -> None:
        """
        C++ signature :void append(std::vector<RDKit::FilterMatch, std::allocator<RDKit::FilterMatch> > {lvalue},boost::python::api::object)
        """
        ...
    def extend(self, arg1: VectFilterMatch, arg2: AtomPairsParameters) -> None:
        """
        C++ signature :void extend(std::vector<RDKit::FilterMatch, std::allocator<RDKit::FilterMatch> > {lvalue},boost::python::api::object)
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

def FilterCatalogCanSerialize(self) -> bool:
    """
    Returns True if the FilterCatalog is serializable (requires boost serialization

    C++ signature :bool FilterCatalogCanSerialize()"""
    ...

def GetFlattenedFunctionalGroupHierarchy(self, normalized: bool = False) -> dict:
    """
    Returns the flattened functional group hierarchy as a dictionary  of name:ROMOL_SPTR substructure items

    C++ signature :boost::python::dict GetFlattenedFunctionalGroupHierarchy([ bool=False])
    """
    ...

def GetFunctionalGroupHierarchy(self) -> FilterCatalog:
    """
    Returns the functional group hierarchy filter catalog

    C++ signature :RDKit::FilterCatalog GetFunctionalGroupHierarchy()"""
    ...

def RunFilterCatalog(
    self,
    filterCatalog: FilterCatalog,
    smiles: _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE,
    numThreads: int = 1,
) -> FilterCatalogListOfEntryList:
    """
    Run the filter catalog on the input list of smiles strings.
    Use numThreads=0 to use all available processors. Returns a vector of vectors.  For each input smiles, a vector of FilterCatalogEntry objects are returned for each matched filter.  If a molecule matches no filter, the vector will be empty. If a smiles string can’t be parsed, a ‘Bad smiles’ entry is returned.

    C++ signature :std::vector<std::vector<boost::shared_ptr<RDKit::FilterCatalogEntry const>, std::allocator<boost::shared_ptr<RDKit::FilterCatalogEntry const> > >, std::allocator<std::vector<boost::shared_ptr<RDKit::FilterCatalogEntry const>, std::allocator<boost::shared_ptr<RDKit::FilterCatalogEntry const> > > > > RunFilterCatalog(RDKit::FilterCatalog,std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > [,int=1])
    """
    ...
