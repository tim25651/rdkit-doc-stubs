"""
rdkit.Avalon.pyAvalonTools module¶
Module containing functionality from the Avalon toolkit.

The functions currently exposed are:
GetCanonSmiles()   : return the canonical smiles for a molecule

GetAvalonFP()return the Avalon fingerprint for a molecule asan RDKit ExplicitBitVector

GetAvalonCountFP()return the Avalon fingerprint for a molecule asan RDKit SparseIntVector

Generate2DCoords()use the Avalon coordinate generator to createa set of 2D coordinates for a molecule

Each function can be called with either an RDKit molecule or some
molecule data as text (e.g. a SMILES or an MDL mol block).
See the individual docstrings for more information.
"""
from typing import Any, ClassVar, overload

import Boost.Python

avalonSSSBits: int
avalonSimilarityBits: int

class StruChkFlag(Boost.Python.enum):
    alias_conversion_failed: StruChkFlag = ...
    atom_check_failed: StruChkFlag = ...
    atom_clash: StruChkFlag = ...
    bad_molecule: StruChkFlag = ...
    dubious_stereo_removed: StruChkFlag = ...
    either_warning: StruChkFlag = ...
    fragments_found: StruChkFlag = ...
    names: dict[str, StruChkFlag] = ...
    recharged: StruChkFlag = ...
    size_check_failed: StruChkFlag = ...
    stereo_error: StruChkFlag = ...
    stereo_forced_bad: StruChkFlag = ...
    stereo_transformed: StruChkFlag = ...
    template_transformed: StruChkFlag = ...
    transformed: StruChkFlag = ...
    values: dict[int, StruChkFlag] = ...
    __slots__: ClassVar[tuple] = ...

class StruChkResult(Boost.Python.enum):
    bad_set: StruChkResult = ...
    names: dict[str, StruChkResult] = ...
    success: StruChkResult = ...
    transformed_set: StruChkResult = ...
    values: dict[int, StruChkResult] = ...
    __slots__: ClassVar[tuple] = ...

@overload
def CheckMolecule(self, molstring: str, isSmiles: bool) -> tuple:
    """
    check a molecule passed in as an RDKit molecule.
    The first member of the return tuple contains the bit-encoded corrections made to the molecule.
    If possible, the molecule (corrected when appropriate) is returned as the second member of
    the return tuple. Otherwise, None is returned.

    C++ signature :boost::python::tuple CheckMolecule(RDKit::ROMol {lvalue})"""
    ...

@overload
def CheckMolecule(self, mol: object) -> tuple: ...
def CheckMoleculeString(self, molstring: str, isSmiles: bool) -> tuple:
    """
    check a molecule passed in as a string and returns the result as a string.
    If the isSmiles argument is true, the string should represent the SMILES encoding
    of the molecule, otherwise it should be encoded as an MDL molfile.
    The first member of the return tuple contains the bit-encoded corrections made to the molecule.
    If possible, a corrected CTAB for the molecule is returned as the second member of
    the return tuple.

    C++ signature :boost::python::tuple CheckMoleculeString(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,bool)
    """
    ...

def CloseCheckMolFiles(self) -> None:
    """
    close open files used by molecule-checking functions.

    C++ signature :void CloseCheckMolFiles()"""
    ...

@overload
def Generate2DCoords(self, mol: object, clearConfs: bool = True) -> int:
    """
    returns an MDL mol block with 2D coordinates for some molecule data.
    If the isSmiles argument is true, the data is assumed to be SMILES, otherwise
    MDL mol data is assumed.

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > Generate2DCoords(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,bool)
    """
    ...

@overload
def Generate2DCoords(self, molData: str, isSmiles: bool) -> str: ...
@overload
def GetAvalonCountFP(
    self, mol: object, nBits: int = 512, isQuery: bool = False, bitFlags: int = 15761407
) -> object:
    """
    returns the Avalon count fingerprint for some molecule data.
    If the isSmiles argument is true, the data is assumed to be SMILES, otherwise
    MDL mol data is assumed.

    C++ signature :RDKit::SparseIntVect<unsigned int>* GetAvalonCountFP(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,bool [,unsigned int=512 [,bool=False [,unsigned int=15761407]]])
    """
    ...

@overload
def GetAvalonCountFP(
    self,
    molData: str,
    isSmiles: bool,
    nBits: int = 512,
    isQuery: bool = False,
    bitFlags: int = 15761407,
) -> object: ...
@overload
def GetAvalonFP(
    self,
    mol: object,
    nBits: int = 512,
    isQuery: bool = False,
    resetVect: bool = False,
    bitFlags: int = 15761407,
) -> object:
    """
    returns the Avalon fingerprint for some molecule data.
    If the isSmiles argument is true, the data is assumed to be SMILES, otherwise
    MDL mol data is assumed.

    C++ signature :ExplicitBitVect* GetAvalonFP(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,bool [,unsigned int=512 [,bool=False [,bool=False [,unsigned int=15761407]]]])
    """
    ...

@overload
def GetAvalonFP(
    self,
    molData: str,
    isSmiles: bool,
    nBits: int = 512,
    isQuery: bool = False,
    resetVect: bool = False,
    bitFlags: int = 15761407,
) -> object: ...
@overload
def GetAvalonFPAsWords(
    self,
    mol: object,
    nBits: int = 512,
    isQuery: bool = False,
    resetVect: bool = False,
    bitFlags: int = 15761407,
) -> list:
    """
    returns the Avalon fingerprint for some molecule data as a list of ints.
    If the isSmiles argument is true, the data is assumed to be SMILES, otherwise
    MDL mol data is assumed.

    C++ signature :boost::python::list GetAvalonFPAsWords(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,bool [,unsigned int=512 [,bool=False [,bool=False [,unsigned int=15761407]]]])
    """
    ...

@overload
def GetAvalonFPAsWords(
    self,
    molData: str,
    isSmiles: bool,
    nBits: int = 512,
    isQuery: bool = False,
    resetVect: bool = False,
    bitFlags: int = 15761407,
) -> list: ...
@overload
def GetCanonSmiles(self, mol: object, flags: int = -1) -> str:
    """
    Returns canonical smiles for some molecule data.
    If the isSmiles argument is true, the data is assumed to be SMILES, otherwise
    MDL mol data is assumed.

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetCanonSmiles(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,bool [,int=-1])
    """
    ...

@overload
def GetCanonSmiles(self, molData: str, isSmiles: bool, flags: int = -1) -> str: ...
def GetCheckMolLog(self) -> str:
    """
    Returns the Struchk log for the last molecules processed.

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetCheckMolLog()
    """
    ...

def InitializeCheckMol(self, options: str = "") -> int:
    """
    initializes the structure checker.
    The argument should contain option lines separated by embedded newlines.An empty string will be used if the argument is omitted.An non-zero error code is returned in case of failure.

    C++ signature :int InitializeCheckMol([ std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’’])
    """
    ...
