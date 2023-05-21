"""
rdkit.Chem.FilterCatalog module¶
"""
from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit.Chem import rdfiltercatalog
from rdkit.Chem.rdfiltercatalog import *
from rdkit.Chem.rdfiltercatalog import FilterMatcherBase

class FilterMatcher(rdfiltercatalog.PythonFilterMatcher):
    """
    FilterMatcher - This class allows creation of Python based
    filters.  Subclass this class to create a Filter useable
    in a FilterCatalogEntry
    Simple Example:
    from rdkit.Chem import rdMolDescriptors
    class MWFilter(FilterMatcher):

    def __init__(self, minMw, maxMw):FilterMatcher.__init__(self, “MW violation”)
    self.minMw = minMw
    self.maxMw = maxMw

    def IsValid(self):return True

    def HasMatch(self, mol):mw = rdMolDescriptors.CalcExactMolWt(mol)
    return not self.minMw <= mw <= self.maxMw

    C++ signature :void __init__(_object*,_object*)"""

    def GetMatches(self, mol, matchVect) -> bool:
        """
        (By default, this calls HasMatch and does not modify matchVect)
        matchVect is a vector of FilterMatch’s which hold the matching
        filter and the matched query_atom, mol_atom pairs if applicable.
        To append to this vector:
        v = MatchTypeVect()
        v.append(IntPair( query_atom_idx, mol_atom_idx ) )
        match = FilterMatch(self, v)
        matchVect.append( match )"""
        ...
    def GetName(self, arg1: FilterMatcherBase) -> str:
        """
        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetName(RDKit::FilterMatcherBase {lvalue})
        """
        ...
    name: Incomplete

    def __init__(self, name: str = ...) -> None: ...
    def HasMatch(self, mol):
        """
        Return True if the filter matches the molecule"""
        ...
    def GetMatches(self, mol, matchVect) -> bool:
        """
        (By default, this calls HasMatch and does not modify matchVect)
        matchVect is a vector of FilterMatch’s which hold the matching
        filter and the matched query_atom, mol_atom pairs if applicable.
        To append to this vector:
        v = MatchTypeVect()
        v.append(IntPair( query_atom_idx, mol_atom_idx ) )
        match = FilterMatch(self, v)
        matchVect.append( match )"""
        ...
    def IsValid(self, mol):
        """
        Must override this function"""
        ...
    def GetName(self, arg1: FilterMatcherBase) -> str:
        """
        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetName(RDKit::FilterMatcherBase {lvalue})
        """
        ...
