"""
rdkit.DataStructs.cDataStructs module¶
Module containing an assortment of functionality for basic data structures.

At the moment the data structures defined are:
Bit Vector classes (for storing signatures, fingerprints and the like:

ExplicitBitVect: class for relatively small (10s of thousands of bits) ordense bit vectors.

SparseBitVect:   class for large, sparse bit vectors

DiscreteValueVect:   class for storing vectors of integers
SparseIntVect:       class for storing sparse vectors of integers
"""
from typing import Any, ClassVar, overload

import Boost.Python
from rdkit.Chem.rdMolDescriptors import AtomPairsParameters
from rdkit.rdBase import _vectd, _vecti

EIGHTBITVALUE: DiscreteValueType
FOURBITVALUE: DiscreteValueType
ONEBITVALUE: DiscreteValueType
SIXTEENBITVALUE: DiscreteValueType
TWOBITVALUE: DiscreteValueType

class DiscreteValueType(Boost.Python.enum):
    EIGHTBITVALUE: DiscreteValueType = ...
    FOURBITVALUE: DiscreteValueType = ...
    ONEBITVALUE: DiscreteValueType = ...
    SIXTEENBITVALUE: DiscreteValueType = ...
    TWOBITVALUE: DiscreteValueType = ...
    names: dict[str, DiscreteValueType] = ...
    values: dict[int, DiscreteValueType] = ...
    __slots__: ClassVar[tuple] = ...

class DiscreteValueVect(Boost.Python.instance):
    """
    A container class for storing unsigned integer
    values within a particular range.
    The length of the vector and type of its elements (determines the maximum value
    that can be stored) are both set at construction time.
    As you would expect, _DiscreteValueVects_ support a set of binary operations
    so you can do things like:

    dvv3 = dvv1 & dvv2  the result contains the smallest value in each entry
    dvv3 = dvv1 | dvv2  the result contains the largest value in each entry
    dvv1 += dvv2     values are truncated when necessary
    dvv3 = dvv1 + dvv2    values are truncated when necessary
    dvv1 -= dvv3    would-be negative values are set to zero
    dvv3 = dvv1 - dvv2    would-be negative values are set to zero

    Elements can be set and read using indexing (i.e. bv[i] = 4 or val=bv[i])
    Constructor

    C++ signature :void __init__(_object*,RDKit::DiscreteValueVect::DiscreteValueType,unsigned int)

    __init__( (object)arg1, (str)arg2) -> None :

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """

    __getstate_manages_dict__: ClassVar[bool] = ...
    __instance_size__: ClassVar[int] = ...
    __safe_for_unpickling__: ClassVar[bool] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def GetTotalVal(self, arg1: DiscreteValueVect) -> int:
        """
        Get the sum of the values in the vector, basically L1 norm

        C++ signature :unsigned int GetTotalVal(RDKit::DiscreteValueVect {lvalue})"""
        ...
    def GetValueType(self, arg1: DiscreteValueVect) -> DiscreteValueType:
        """
        Get the type of value stored in the vector

        C++ signature :RDKit::DiscreteValueVect::DiscreteValueType GetValueType(RDKit::DiscreteValueVect {lvalue})
        """
        ...
    @classmethod
    def __add__(cls, other) -> Any: ...
    @classmethod
    def __and__(cls, other) -> Any: ...
    @classmethod
    def __getinitargs__(cls, RDKit) -> Any: ...
    @classmethod
    def __getitem__(cls, RDKit, unsignedint) -> Any: ...
    @classmethod
    def __getstate__(cls, boost) -> Any: ...
    @classmethod
    def __iadd__(cls, boost, RDKit) -> Any: ...
    @classmethod
    def __isub__(cls, boost, RDKit) -> Any: ...
    @classmethod
    def __len__(cls, RDKit) -> Any: ...
    @classmethod
    def __or__(cls, other) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @classmethod
    def __setitem__(cls, index, object) -> Any: ...
    @classmethod
    def __setstate__(cls, state) -> Any: ...
    @classmethod
    def __sub__(cls, other) -> Any: ...

class ExplicitBitVect(Boost.Python.instance):
    """
    A class to store explicit bit vectors.
    This class is most useful for situations where the size of the vector
    is relatively small (tens of thousands or smaller).
    For larger vectors, use the _SparseBitVect_ class instead.
    As you would expect, _ExplicitBitVects_ support a set of binary operations
    so you can do things like:

    bv3 = bv1 & bv2  (bitwise and)
    bv3 = bv1 | bv2  (bitwise or)
    bv3 = bv1 ^ bv2  (bitwise xor)
    bv3 = ~bv1       (bitwise negation)

    Bits can be set and read using either the Set/UnsetBit() and GetBit() methods
    or by indexing (i.e. bv[i] = 1 or if bv[i]).

    C++ signature :void __init__(_object*,unsigned int)

    __init__( (object)arg1, (str)arg2) -> None :

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)

    __init__( (object)arg1, (int)size, (bool)bitsSet) -> None :

    C++ signature :void __init__(_object*,unsigned int,bool)"""

    __getstate_manages_dict__: ClassVar[bool] = ...
    __instance_size__: ClassVar[int] = ...
    __safe_for_unpickling__: ClassVar[bool] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def FromBase64(self, arg1: ExplicitBitVect, arg2: str) -> None:
        """
        Initializes the vector from a base64 encoded binary string.

        C++ signature :void FromBase64(ExplicitBitVect {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    def GetBit(self, arg1: ExplicitBitVect, arg2: int) -> bool:
        """
        Returns the value of a bit.

        C++ signature :bool GetBit(ExplicitBitVect {lvalue},unsigned int)"""
        ...
    def GetNumBits(self, arg1: ExplicitBitVect) -> int:
        """
        Returns the number of bits in the vector (the vector’s size).

        C++ signature :unsigned int GetNumBits(ExplicitBitVect {lvalue})"""
        ...
    def GetNumOffBits(self, arg1: ExplicitBitVect) -> int:
        """
        Returns the number of off bits.

        C++ signature :unsigned int GetNumOffBits(ExplicitBitVect {lvalue})"""
        ...
    def GetNumOnBits(self, arg1: ExplicitBitVect) -> int:
        """
        Returns the number of on bits.

        C++ signature :unsigned int GetNumOnBits(ExplicitBitVect {lvalue})"""
        ...
    def GetOnBits(self, arg1: ExplicitBitVect) -> _vecti:
        """
        Returns a tuple containing IDs of the on bits.

        C++ signature :std::vector<int, std::allocator<int> > GetOnBits(ExplicitBitVect)
        """
        ...
    def SetBit(self, arg1: ExplicitBitVect, arg2: int) -> bool:
        """
        Turns on a particular bit.  Returns the original state of the bit.

        C++ signature :bool SetBit(ExplicitBitVect {lvalue},unsigned int)"""
        ...
    def SetBitsFromList(self, arg1: ExplicitBitVect, arg2: AtomPairsParameters) -> None:
        """
        Turns on a set of bits.  The argument should be a tuple or list of bit ids.

        C++ signature :void SetBitsFromList(ExplicitBitVect*,boost::python::api::object)
        """
        ...
    def ToBase64(self, arg1: ExplicitBitVect) -> str:
        """
        Converts the vector to a base64 string (the base64 encoded version of the results of ToString()).

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > ToBase64(ExplicitBitVect {lvalue})
        """
        ...
    def ToBinary(self, arg1: ExplicitBitVect) -> object:
        """
        Returns an internal binary representation of the vector.

        C++ signature :boost::python::api::object ToBinary(ExplicitBitVect)"""
        ...
    def ToBitString(self):
        """
        BitVectToText( (SparseBitVect)arg1) -> str :

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > BitVectToText(SparseBitVect)

        BitVectToText( (ExplicitBitVect)arg1) -> str :Returns a string of zeros and ones representing the bit vector.

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > BitVectToText(ExplicitBitVect)
        """
        ...
    def ToList(self, arg1: ExplicitBitVect) -> list:
        """
        Return the Bitvector as a python list (faster than list(vect))

        C++ signature :boost::python::list ToList(ExplicitBitVect)"""
        ...
    def UnSetBit(self, arg1: ExplicitBitVect, arg2: int) -> bool:
        """
        Turns off a particular bit.  Returns the original state of the bit.

        C++ signature :bool UnSetBit(ExplicitBitVect {lvalue},unsigned int)"""
        ...
    def UnSetBitsFromList(
        self, arg1: ExplicitBitVect, arg2: AtomPairsParameters
    ) -> None:
        """
        Turns off a set of bits.  The argument should be a tuple or list of bit ids.

        C++ signature :void UnSetBitsFromList(ExplicitBitVect*,boost::python::api::object)
        """
        ...
    @classmethod
    def __add__(cls, other) -> Any: ...
    @classmethod
    def __and__(cls, other) -> Any: ...
    @classmethod
    def __eq__(cls, other) -> Any: ...
    @classmethod
    def __getinitargs__(cls, ExplicitBitVect) -> Any: ...
    @classmethod
    def __getitem__(cls, ExplicitBitVect, int) -> Any: ...
    @classmethod
    def __getstate__(cls, boost) -> Any: ...
    @classmethod
    def __iadd__(cls, boost, ExplicitBitVect) -> Any: ...
    @classmethod
    def __invert__(cls) -> Any: ...
    @classmethod
    def __len__(cls) -> Any: ...
    @classmethod
    def __ne__(cls, other) -> Any: ...
    @classmethod
    def __or__(cls, other) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @classmethod
    def __setitem__(cls, index, object) -> Any: ...
    @classmethod
    def __setstate__(cls, state) -> Any: ...
    @classmethod
    def __xor__(cls, other) -> Any: ...

class FPBReader(Boost.Python.instance):
    """
    A class for reading and searching FPB files from Andrew Dalke’s chemfp.
    Note that this functionality is still experimental and the API may
    change in future releases.
    docstring

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])
    """

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def GetBytes(self, arg1: FPBReader, arg2: int) -> object:
        """
        returns a particular fingerprint as bytes

        C++ signature :boost::python::api::object GetBytes(RDKit::FPBReader const*,unsigned int)
        """
        ...
    def GetContainingNeighbors(self, arg1: FPBReader, bv: str) -> tuple:
        """
        returns indices of neighbors that contain this fingerprint (where all bits from this fingerprint are also set)

        C++ signature :boost::python::tuple GetContainingNeighbors(RDKit::FPBReader const*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    def GetFP(self, arg1: FPBReader, arg2: int) -> ExplicitBitVect:
        """
        returns a particular fingerprint as an ExplicitBitVect

        C++ signature :boost::shared_ptr<ExplicitBitVect> GetFP(RDKit::FPBReader {lvalue},unsigned int)
        """
        ...
    def GetId(self, arg1: FPBReader, arg2: int) -> str:
        """
        returns the id of a particular fingerprint

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetId(RDKit::FPBReader {lvalue},unsigned int)
        """
        ...
    def GetNumBits(self, arg1: FPBReader) -> int:
        """
        returns the number of bits in a fingerprint

        C++ signature :unsigned int GetNumBits(RDKit::FPBReader {lvalue})"""
        ...
    def GetTanimoto(self, arg1: FPBReader, arg2: int, arg3: str) -> float:
        """
        return the tanimoto similarity of a particular fingerprint to the bytes provided

        C++ signature :double GetTanimoto(RDKit::FPBReader const*,unsigned int,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    def GetTanimotoNeighbors(
        self, arg1: FPBReader, bv: str, threshold: float = 0.7
    ) -> tuple:
        """
        returns tanimoto similarities to and indices of all neighbors above the specified threshold

        C++ signature :boost::python::tuple GetTanimotoNeighbors(RDKit::FPBReader const*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,double=0.7])
        """
        ...
    def GetTversky(
        self, arg1: FPBReader, arg2: int, arg3: str, arg4: float, arg5: float
    ) -> float:
        """
        return the Tverksy similarity of a particular fingerprint to the bytes provided

        C++ signature :double GetTversky(RDKit::FPBReader const*,unsigned int,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,double,double)
        """
        ...
    def GetTverskyNeighbors(
        self, arg1: FPBReader, bv: str, ca: float, cb: float, threshold: float = 0.7
    ) -> tuple:
        """
        returns Tversky similarities to and indices of all neighbors above the specified threshold

        C++ signature :boost::python::tuple GetTverskyNeighbors(RDKit::FPBReader const*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,double,double [,double=0.7])
        """
        ...
    def Init(self, arg1: FPBReader) -> None:
        """
        Read the fingerprints from the file. This can take a while.

        C++ signature :void Init(RDKit::FPBReader {lvalue})"""
        ...
    @classmethod
    def __getitem__(cls, RDKit, unsignedint) -> Any: ...
    @classmethod
    def __len__(cls, RDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class IntSparseIntVect(Boost.Python.instance):
    """
    A container class for storing integer
    values within a particular range.
    The length of the vector is set at construction time.
    As you would expect, _SparseIntVects_ support a set of binary operations
    so you can do things like:

    Arithmetic:
    siv1 += siv2
    siv3 = siv1 + siv2
    siv1 -= siv3
    siv3 = siv1 - siv2
    “Fuzzy” binary operations:
    siv3 = siv1 & siv2  the result contains the smallest value in each entry
    siv3 = siv1 | siv2  the result contains the largest value in each entry

    Elements can be set and read using indexing (i.e. siv[i] = 4 or val=siv[i])
    Constructor

    C++ signature :void __init__(_object*,int)

    __init__( (object)arg1, (str)arg2) -> None :

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """

    __getstate_manages_dict__: ClassVar[bool] = ...
    __instance_size__: ClassVar[int] = ...
    __safe_for_unpickling__: ClassVar[bool] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def GetLength(self, arg1: IntSparseIntVect) -> int:
        """
        Returns the length of the vector

        C++ signature :int GetLength(RDKit::SparseIntVect<int> {lvalue})"""
        ...
    def GetNonzeroElements(self, arg1: IntSparseIntVect) -> dict:
        """
        returns a dictionary of the nonzero elements

        C++ signature :boost::python::dict GetNonzeroElements(RDKit::SparseIntVect<int> {lvalue})
        """
        ...
    def GetTotalVal(self, arg1: IntSparseIntVect, useAbs: bool = False) -> int:
        """
        Get the sum of the values in the vector, basically L1 norm

        C++ signature :int GetTotalVal(RDKit::SparseIntVect<int> {lvalue} [,bool=False])
        """
        ...
    def ToBinary(self, arg1: IntSparseIntVect) -> object:
        """
        returns a binary (pickle) representation of the vector

        C++ signature :boost::python::api::object ToBinary(RDKit::SparseIntVect<int>)"""
        ...
    def ToList(self, arg1: IntSparseIntVect) -> list:
        """
        Return the SparseIntVect as a python list

        C++ signature :boost::python::list ToList(RDKit::SparseIntVect<int> {lvalue})"""
        ...
    def UpdateFromSequence(
        self, arg1: IntSparseIntVect, arg2: AtomPairsParameters
    ) -> None:
        """
        update the vector based on the values in the list or tuple

        C++ signature :void UpdateFromSequence(RDKit::SparseIntVect<int> {lvalue},boost::python::api::object {lvalue})
        """
        ...
    @classmethod
    def __add__(cls, other) -> Any: ...
    @classmethod
    def __and__(cls, other) -> Any: ...
    @classmethod
    def __eq__(cls, other) -> Any: ...
    @classmethod
    def __getinitargs__(cls, RDKit) -> Any: ...
    @classmethod
    def __getitem__(cls, RDKit, int) -> Any: ...
    @classmethod
    def __getstate__(cls, boost) -> Any: ...
    @overload
    @classmethod
    def __iadd__(cls, boost, RDKit) -> Any: ...
    @overload
    @classmethod
    def __iadd__(cls, boost, int) -> Any: ...
    @classmethod
    def __idiv__(cls, boost, int) -> Any: ...
    @classmethod
    def __imul__(cls, boost, int) -> Any: ...
    @overload
    @classmethod
    def __isub__(cls, boost, RDKit) -> Any: ...
    @overload
    @classmethod
    def __isub__(cls, boost, int) -> Any: ...
    @classmethod
    def __ne__(cls, other) -> Any: ...
    @classmethod
    def __or__(cls, other) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @classmethod
    def __setitem__(cls, index, object) -> Any: ...
    @classmethod
    def __setstate__(cls, state) -> Any: ...
    @classmethod
    def __sub__(cls, other) -> Any: ...

class LongSparseIntVect(Boost.Python.instance):
    """
    A container class for storing integer
    values within a particular range.
    The length of the vector is set at construction time.
    As you would expect, _SparseIntVects_ support a set of binary operations
    so you can do things like:

    Arithmetic:
    siv1 += siv2
    siv3 = siv1 + siv2
    siv1 -= siv3
    siv3 = siv1 - siv2
    “Fuzzy” binary operations:
    siv3 = siv1 & siv2  the result contains the smallest value in each entry
    siv3 = siv1 | siv2  the result contains the largest value in each entry

    Elements can be set and read using indexing (i.e. siv[i] = 4 or val=siv[i])
    Constructor

    C++ signature :void __init__(_object*,long)

    __init__( (object)arg1, (str)arg2) -> None :

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """

    __getstate_manages_dict__: ClassVar[bool] = ...
    __instance_size__: ClassVar[int] = ...
    __safe_for_unpickling__: ClassVar[bool] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def GetLength(self, arg1: LongSparseIntVect) -> int:
        """
        Returns the length of the vector

        C++ signature :long GetLength(RDKit::SparseIntVect<long> {lvalue})"""
        ...
    def GetNonzeroElements(self, arg1: LongSparseIntVect) -> dict:
        """
        returns a dictionary of the nonzero elements

        C++ signature :boost::python::dict GetNonzeroElements(RDKit::SparseIntVect<long> {lvalue})
        """
        ...
    def GetTotalVal(self, arg1: LongSparseIntVect, useAbs: bool = False) -> int:
        """
        Get the sum of the values in the vector, basically L1 norm

        C++ signature :int GetTotalVal(RDKit::SparseIntVect<long> {lvalue} [,bool=False])
        """
        ...
    def ToBinary(self, arg1: LongSparseIntVect) -> object:
        """
        returns a binary (pickle) representation of the vector

        C++ signature :boost::python::api::object ToBinary(RDKit::SparseIntVect<long>)
        """
        ...
    def ToList(self, arg1: LongSparseIntVect) -> list:
        """
        Return the SparseIntVect as a python list

        C++ signature :boost::python::list ToList(RDKit::SparseIntVect<long> {lvalue})
        """
        ...
    def UpdateFromSequence(
        self, arg1: LongSparseIntVect, arg2: AtomPairsParameters
    ) -> None:
        """
        update the vector based on the values in the list or tuple

        C++ signature :void UpdateFromSequence(RDKit::SparseIntVect<long> {lvalue},boost::python::api::object {lvalue})
        """
        ...
    @classmethod
    def __add__(cls, other) -> Any: ...
    @classmethod
    def __and__(cls, other) -> Any: ...
    @classmethod
    def __eq__(cls, other) -> Any: ...
    @classmethod
    def __getinitargs__(cls, RDKit) -> Any: ...
    @classmethod
    def __getitem__(cls, RDKit, long) -> Any: ...
    @classmethod
    def __getstate__(cls, boost) -> Any: ...
    @overload
    @classmethod
    def __iadd__(cls, boost, RDKit) -> Any: ...
    @overload
    @classmethod
    def __iadd__(cls, boost, int) -> Any: ...
    @classmethod
    def __idiv__(cls, boost, int) -> Any: ...
    @classmethod
    def __imul__(cls, boost, int) -> Any: ...
    @overload
    @classmethod
    def __isub__(cls, boost, RDKit) -> Any: ...
    @overload
    @classmethod
    def __isub__(cls, boost, int) -> Any: ...
    @classmethod
    def __ne__(cls, other) -> Any: ...
    @classmethod
    def __or__(cls, other) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @classmethod
    def __setitem__(cls, RDKit, long, int) -> Any: ...
    @classmethod
    def __setstate__(cls, state) -> Any: ...
    @classmethod
    def __sub__(cls, other) -> Any: ...

class MultiFPBReader(Boost.Python.instance):
    """
    A class for reading and searching multiple FPB files from Andrew Dalke’s chemfp.
    Note that this functionality is still experimental and the API may
    change in future releases.
    docstring

    C++ signature :void __init__(_object* [,bool=False])"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def AddReader(self, arg1: MultiFPBReader, arg2: FPBReader) -> int:
        """
        adds an FPBReader to our set of readers

        C++ signature :unsigned int AddReader(RDKit::MultiFPBReader {lvalue},RDKit::FPBReader*)
        """
        ...
    def GetContainingNeighbors(
        self, arg1: MultiFPBReader, bv: str, numThreads: int = 1
    ) -> tuple:
        """
        returns indices of neighbors that contain this fingerprint (where all bits from this fingerprint are also set)

        C++ signature :boost::python::tuple GetContainingNeighbors(RDKit::MultiFPBReader const*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,unsigned int=1])
        """
        ...
    def GetNumBits(self, arg1: MultiFPBReader) -> int:
        """
        returns the number of bits in a fingerprint

        C++ signature :unsigned int GetNumBits(RDKit::MultiFPBReader {lvalue})"""
        ...
    def GetReader(self, arg1: MultiFPBReader, arg2: int) -> FPBReader:
        """
        returns one of our readers

        C++ signature :RDKit::FPBReader* GetReader(RDKit::MultiFPBReader {lvalue},unsigned int)
        """
        ...
    def GetTanimotoNeighbors(
        self, arg1: MultiFPBReader, bv: str, threshold: float = 0.7, numThreads: int = 1
    ) -> tuple:
        """
        returns tanimoto similarities to and indices of all neighbors above the specified threshold

        C++ signature :boost::python::tuple GetTanimotoNeighbors(RDKit::MultiFPBReader const*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,double=0.7 [,unsigned int=1]])
        """
        ...
    def GetTverskyNeighbors(
        self,
        arg1: MultiFPBReader,
        bv: str,
        ca: float,
        cb: float,
        threshold: float = 0.7,
        numThreads: int = 1,
    ) -> tuple:
        """
        returns Tversky similarities to and indices of all neighbors above the specified threshold

        C++ signature :boost::python::tuple GetTverskyNeighbors(RDKit::MultiFPBReader const*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,double,double [,double=0.7 [,unsigned int=1]])
        """
        ...
    def Init(self, arg1: MultiFPBReader) -> None:
        """
        Call Init() on each of our children. This can take a while.

        C++ signature :void Init(RDKit::MultiFPBReader {lvalue})"""
        ...
    @classmethod
    def __len__(cls, RDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class SparseBitVect(Boost.Python.instance):
    """
    A class to store sparse bit vectors.
    This class is most useful for situations where the size of the vector
    is large and relatively few bits are set
    For smaller or denser vectors, the _ExplicitBitVect_ class is much faster.
    As you would expect, _SparseBitVects_ support a set of binary operations
    so you can do things like:

    bv3 = bv1 & bv2  (bitwise and)
    bv3 = bv1 | bv2  (bitwise or)
    bv3 = bv1 ^ bv2  (bitwise xor)
    bv3 = ~bv1       (bitwise negation) NOTE: this operation is likely

    to be VERY slow and inefficient.

    Bits can be set and read using either the Set/UnsetBit() and GetBit() methods
    or by indexing (i.e. bv[i] = 1 or if bv[i]).

    C++ signature :void __init__(_object*,unsigned int)

    __init__( (object)arg1, (str)arg2) -> None :

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """

    __getstate_manages_dict__: ClassVar[bool] = ...
    __instance_size__: ClassVar[int] = ...
    __safe_for_unpickling__: ClassVar[bool] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def FromBase64(self, arg1: SparseBitVect, arg2: str) -> None:
        """
        Initializes the vector from a base64 encoded binary string.

        C++ signature :void FromBase64(SparseBitVect {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    def GetBit(self, arg1: SparseBitVect, arg2: int) -> bool:
        """
        Returns the value of a bit.

        C++ signature :bool GetBit(SparseBitVect {lvalue},unsigned int)"""
        ...
    def GetNumBits(self, arg1: SparseBitVect) -> int:
        """
        Returns the number of bits in the vector (the vector’s size).

        C++ signature :unsigned int GetNumBits(SparseBitVect {lvalue})"""
        ...
    def GetNumOffBits(self, arg1: SparseBitVect) -> int:
        """
        Returns the number of off bits.

        C++ signature :unsigned int GetNumOffBits(SparseBitVect {lvalue})"""
        ...
    def GetNumOnBits(self, arg1: SparseBitVect) -> int:
        """
        Returns the number of on bits.

        C++ signature :unsigned int GetNumOnBits(SparseBitVect {lvalue})"""
        ...
    def GetOnBits(self, arg1: SparseBitVect) -> _vecti:
        """
        Returns a tuple containing IDs of the on bits.

        C++ signature :std::vector<int, std::allocator<int> > GetOnBits(SparseBitVect)
        """
        ...
    def SetBit(self, arg1: SparseBitVect, arg2: int) -> bool:
        """
        Turns on a particular bit.  Returns the original state of the bit.

        C++ signature :bool SetBit(SparseBitVect {lvalue},unsigned int)"""
        ...
    def SetBitsFromList(self, arg1: SparseBitVect, arg2: AtomPairsParameters) -> None:
        """
        Turns on a set of bits.  The argument should be a tuple or list of bit ids.

        C++ signature :void SetBitsFromList(SparseBitVect*,boost::python::api::object)
        """
        ...
    def ToBase64(self, arg1: SparseBitVect) -> str:
        """
        Converts the vector to a base64 string (the base64 encoded version of the results of ToString()).

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > ToBase64(SparseBitVect {lvalue})
        """
        ...
    def ToBinary(self, arg1: SparseBitVect) -> object:
        """
        Returns an internal binary representation of the vector.

        C++ signature :boost::python::api::object ToBinary(SparseBitVect)"""
        ...
    def ToBitString(self):
        """
        BitVectToText( (SparseBitVect)arg1) -> str :

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > BitVectToText(SparseBitVect)

        BitVectToText( (ExplicitBitVect)arg1) -> str :Returns a string of zeros and ones representing the bit vector.

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > BitVectToText(ExplicitBitVect)
        """
        ...
    def ToList(self, arg1: SparseBitVect) -> list:
        """
        Return the BitVector as a python list

        C++ signature :boost::python::list ToList(SparseBitVect)"""
        ...
    def UnSetBit(self, arg1: SparseBitVect, arg2: int) -> bool:
        """
        Turns off a particular bit.  Returns the original state of the bit.

        C++ signature :bool UnSetBit(SparseBitVect {lvalue},unsigned int)"""
        ...
    def UnSetBitsFromList(self, arg1: SparseBitVect, arg2: AtomPairsParameters) -> None:
        """
        Turns off a set of bits.  The argument should be a tuple or list of bit ids.

        C++ signature :void UnSetBitsFromList(SparseBitVect*,boost::python::api::object)
        """
        ...
    @classmethod
    def __and__(cls, other) -> Any: ...
    @classmethod
    def __eq__(cls, other) -> Any: ...
    @classmethod
    def __getinitargs__(cls, SparseBitVect) -> Any: ...
    @classmethod
    def __getitem__(cls, SparseBitVect, int) -> Any: ...
    @classmethod
    def __getstate__(cls, boost) -> Any: ...
    @classmethod
    def __invert__(cls) -> Any: ...
    @classmethod
    def __len__(cls) -> Any: ...
    @classmethod
    def __ne__(cls, other) -> Any: ...
    @classmethod
    def __or__(cls, other) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @classmethod
    def __setitem__(cls, index, object) -> Any: ...
    @classmethod
    def __setstate__(cls, state) -> Any: ...
    @classmethod
    def __xor__(cls, other) -> Any: ...

class UIntSparseIntVect(Boost.Python.instance):
    """
    A container class for storing integer
    values within a particular range.
    The length of the vector is set at construction time.
    As you would expect, _SparseIntVects_ support a set of binary operations
    so you can do things like:

    Arithmetic:
    siv1 += siv2
    siv3 = siv1 + siv2
    siv1 -= siv3
    siv3 = siv1 - siv2
    “Fuzzy” binary operations:
    siv3 = siv1 & siv2  the result contains the smallest value in each entry
    siv3 = siv1 | siv2  the result contains the largest value in each entry

    Elements can be set and read using indexing (i.e. siv[i] = 4 or val=siv[i])
    Constructor

    C++ signature :void __init__(_object*,unsigned int)

    __init__( (object)arg1, (str)arg2) -> None :

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """

    __getstate_manages_dict__: ClassVar[bool] = ...
    __instance_size__: ClassVar[int] = ...
    __safe_for_unpickling__: ClassVar[bool] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def GetLength(self, arg1: UIntSparseIntVect) -> int:
        """
        Returns the length of the vector

        C++ signature :unsigned int GetLength(RDKit::SparseIntVect<unsigned int> {lvalue})
        """
        ...
    def GetNonzeroElements(self, arg1: UIntSparseIntVect) -> dict:
        """
        returns a dictionary of the nonzero elements

        C++ signature :boost::python::dict GetNonzeroElements(RDKit::SparseIntVect<unsigned int> {lvalue})
        """
        ...
    def GetTotalVal(self, arg1: UIntSparseIntVect, useAbs: bool = False) -> int:
        """
        Get the sum of the values in the vector, basically L1 norm

        C++ signature :int GetTotalVal(RDKit::SparseIntVect<unsigned int> {lvalue} [,bool=False])
        """
        ...
    def ToBinary(self, arg1: UIntSparseIntVect) -> object:
        """
        returns a binary (pickle) representation of the vector

        C++ signature :boost::python::api::object ToBinary(RDKit::SparseIntVect<unsigned int>)
        """
        ...
    def ToList(self, arg1: UIntSparseIntVect) -> list:
        """
        Return the SparseIntVect as a python list

        C++ signature :boost::python::list ToList(RDKit::SparseIntVect<unsigned int> {lvalue})
        """
        ...
    def UpdateFromSequence(
        self, arg1: UIntSparseIntVect, arg2: AtomPairsParameters
    ) -> None:
        """
        update the vector based on the values in the list or tuple

        C++ signature :void UpdateFromSequence(RDKit::SparseIntVect<unsigned int> {lvalue},boost::python::api::object {lvalue})
        """
        ...
    @classmethod
    def __add__(cls, other) -> Any: ...
    @classmethod
    def __and__(cls, other) -> Any: ...
    @classmethod
    def __eq__(cls, other) -> Any: ...
    @classmethod
    def __getinitargs__(cls, RDKit) -> Any: ...
    @classmethod
    def __getitem__(cls, RDKit, unsignedint) -> Any: ...
    @classmethod
    def __getstate__(cls, boost) -> Any: ...
    @overload
    @classmethod
    def __iadd__(cls, boost, RDKit) -> Any: ...
    @overload
    @classmethod
    def __iadd__(cls, boost, int) -> Any: ...
    @classmethod
    def __idiv__(cls, boost, int) -> Any: ...
    @classmethod
    def __imul__(cls, boost, int) -> Any: ...
    @overload
    @classmethod
    def __isub__(cls, boost, RDKit) -> Any: ...
    @overload
    @classmethod
    def __isub__(cls, boost, int) -> Any: ...
    @classmethod
    def __ne__(cls, other) -> Any: ...
    @classmethod
    def __or__(cls, other) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @classmethod
    def __setitem__(cls, RDKit, unsignedint, int) -> Any: ...
    @classmethod
    def __setstate__(cls, state) -> Any: ...
    @classmethod
    def __sub__(cls, other) -> Any: ...

class ULongSparseIntVect(Boost.Python.instance):
    """
    A container class for storing integer
    values within a particular range.
    The length of the vector is set at construction time.
    As you would expect, _SparseIntVects_ support a set of binary operations
    so you can do things like:

    Arithmetic:
    siv1 += siv2
    siv3 = siv1 + siv2
    siv1 -= siv3
    siv3 = siv1 - siv2
    “Fuzzy” binary operations:
    siv3 = siv1 & siv2  the result contains the smallest value in each entry
    siv3 = siv1 | siv2  the result contains the largest value in each entry

    Elements can be set and read using indexing (i.e. siv[i] = 4 or val=siv[i])
    Constructor

    C++ signature :void __init__(_object*,unsigned long)

    __init__( (object)arg1, (str)arg2) -> None :

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """

    __getstate_manages_dict__: ClassVar[bool] = ...
    __instance_size__: ClassVar[int] = ...
    __safe_for_unpickling__: ClassVar[bool] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def GetLength(self, arg1: ULongSparseIntVect) -> int:
        """
        Returns the length of the vector

        C++ signature :unsigned long GetLength(RDKit::SparseIntVect<unsigned long> {lvalue})
        """
        ...
    def GetNonzeroElements(self, arg1: ULongSparseIntVect) -> dict:
        """
        returns a dictionary of the nonzero elements

        C++ signature :boost::python::dict GetNonzeroElements(RDKit::SparseIntVect<unsigned long> {lvalue})
        """
        ...
    def GetTotalVal(self, arg1: ULongSparseIntVect, useAbs: bool = False) -> int:
        """
        Get the sum of the values in the vector, basically L1 norm

        C++ signature :int GetTotalVal(RDKit::SparseIntVect<unsigned long> {lvalue} [,bool=False])
        """
        ...
    def ToBinary(self, arg1: ULongSparseIntVect) -> object:
        """
        returns a binary (pickle) representation of the vector

        C++ signature :boost::python::api::object ToBinary(RDKit::SparseIntVect<unsigned long>)
        """
        ...
    def ToList(self, arg1: ULongSparseIntVect) -> list:
        """
        Return the SparseIntVect as a python list

        C++ signature :boost::python::list ToList(RDKit::SparseIntVect<unsigned long> {lvalue})
        """
        ...
    def UpdateFromSequence(
        self, arg1: ULongSparseIntVect, arg2: AtomPairsParameters
    ) -> None:
        """
        update the vector based on the values in the list or tuple

        C++ signature :void UpdateFromSequence(RDKit::SparseIntVect<unsigned long> {lvalue},boost::python::api::object {lvalue})
        """
        ...
    @classmethod
    def __add__(cls, other) -> Any: ...
    @classmethod
    def __and__(cls, other) -> Any: ...
    @classmethod
    def __eq__(cls, other) -> Any: ...
    @classmethod
    def __getinitargs__(cls, RDKit) -> Any: ...
    @classmethod
    def __getitem__(cls, RDKit, unsignedlong) -> Any: ...
    @classmethod
    def __getstate__(cls, boost) -> Any: ...
    @overload
    @classmethod
    def __iadd__(cls, boost, RDKit) -> Any: ...
    @overload
    @classmethod
    def __iadd__(cls, boost, int) -> Any: ...
    @classmethod
    def __idiv__(cls, boost, int) -> Any: ...
    @classmethod
    def __imul__(cls, boost, int) -> Any: ...
    @overload
    @classmethod
    def __isub__(cls, boost, RDKit) -> Any: ...
    @overload
    @classmethod
    def __isub__(cls, boost, int) -> Any: ...
    @classmethod
    def __ne__(cls, other) -> Any: ...
    @classmethod
    def __or__(cls, other) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @classmethod
    def __setitem__(cls, RDKit, unsignedlong, int) -> Any: ...
    @classmethod
    def __setstate__(cls, state) -> Any: ...
    @classmethod
    def __sub__(cls, other) -> Any: ...

@overload
def AllBitSimilarity(self, v1: SparseBitVect, v2: SparseBitVect) -> float:
    """
    (B(bv1) - B(bv1^bv2)) / B(bv1)

    C++ signature :double AllBitSimilarity(ExplicitBitVect,ExplicitBitVect)"""
    ...

@overload
def AllBitSimilarity(self, v1: ExplicitBitVect, v2: ExplicitBitVect) -> float: ...
@overload
def AllProbeBitsMatch(self, arg1: SparseBitVect, arg2: SparseBitVect) -> bool:
    """


    C++ signature :bool AllProbeBitsMatch(ExplicitBitVect,ExplicitBitVect)

    AllProbeBitsMatch( (SparseBitVect)arg1, (str)arg2) -> bool :

    C++ signature :bool AllProbeBitsMatch(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)

    AllProbeBitsMatch( (ExplicitBitVect)arg1, (str)arg2) -> bool :
    Returns True if all bits in the first argument match all bits in the vector defined by the pickle in the second argument.

    C++ signature :bool AllProbeBitsMatch(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
    ...

@overload
def AllProbeBitsMatch(self, arg1: ExplicitBitVect, arg2: ExplicitBitVect) -> bool:
    """


    C++ signature :bool AllProbeBitsMatch(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)

    AllProbeBitsMatch( (ExplicitBitVect)arg1, (str)arg2) -> bool :
    Returns True if all bits in the first argument match all bits in the vector defined by the pickle in the second argument.

    C++ signature :bool AllProbeBitsMatch(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
    ...

@overload
def AllProbeBitsMatch(self, arg1: SparseBitVect, arg2: str) -> bool:
    """

    Returns True if all bits in the first argument match all bits in the vector defined by the pickle in the second argument.

    C++ signature :bool AllProbeBitsMatch(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
    ...

@overload
def AllProbeBitsMatch(self, arg1: ExplicitBitVect, arg2: str) -> bool: ...
@overload
def AsymmetricSimilarity(
    self, bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0
) -> float:
    """
    B(bv1&bv2) / min(B(bv1),B(bv2))

    C++ signature :double AsymmetricSimilarity(ExplicitBitVect,ExplicitBitVect [,bool=0])

    AsymmetricSimilarity( (SparseBitVect)bv1, (str)pkl [, (bool)returnDistance=0]) -> float :

    C++ signature :double AsymmetricSimilarity(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

    AsymmetricSimilarity( (ExplicitBitVect)bv1, (str)pkl [, (bool)returnDistance=0]) -> float :B(bv1&bv2) / min(B(bv1),B(bv2))

    C++ signature :double AsymmetricSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])
    """
    ...

@overload
def AsymmetricSimilarity(
    self, bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = 0
) -> float:
    """


    C++ signature :double AsymmetricSimilarity(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

    AsymmetricSimilarity( (ExplicitBitVect)bv1, (str)pkl [, (bool)returnDistance=0]) -> float :B(bv1&bv2) / min(B(bv1),B(bv2))

    C++ signature :double AsymmetricSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])
    """
    ...

@overload
def AsymmetricSimilarity(
    self, bv1: SparseBitVect, pkl: str, returnDistance: bool = 0
) -> float:
    """
    B(bv1&bv2) / min(B(bv1),B(bv2))

    C++ signature :double AsymmetricSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])
    """
    ...

@overload
def AsymmetricSimilarity(
    self, bv1: ExplicitBitVect, pkl: str, returnDistance: bool = 0
) -> float: ...
def AsymmetricSimilarityNeighbors(
    self, bvqueries: AtomPairsParameters, bvList: AtomPairsParameters
) -> list:
    """
    B(bv1&bv2) / min(B(bv1),B(bv2))

    C++ signature :boost::python::list AsymmetricSimilarityNeighbors(boost::python::api::object,boost::python::api::object)
    """
    ...

def AsymmetricSimilarityNeighbors_sparse(
    self, bvqueries: AtomPairsParameters, bvList: AtomPairsParameters
) -> list:
    """
    B(bv1&bv2) / min(B(bv1),B(bv2))

    C++ signature :boost::python::list AsymmetricSimilarityNeighbors_sparse(boost::python::api::object,boost::python::api::object)
    """
    ...

@overload
def BitVectToBinaryText(self, arg1: SparseBitVect) -> object:
    """
    Returns a binary string (byte array) representing the bit vector.

    C++ signature :boost::python::api::object BitVectToBinaryText(ExplicitBitVect)"""
    ...

@overload
def BitVectToBinaryText(self, arg1: ExplicitBitVect) -> object: ...
@overload
def BitVectToFPSText(self, arg1: SparseBitVect) -> str:
    """
    Returns an FPS string representing the bit vector.

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > BitVectToFPSText(ExplicitBitVect)
    """
    ...

@overload
def BitVectToFPSText(self, arg1: ExplicitBitVect) -> str: ...
@overload
def BitVectToText(self, arg1: SparseBitVect) -> str:
    """
    Returns a string of zeros and ones representing the bit vector.

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > BitVectToText(ExplicitBitVect)
    """
    ...

@overload
def BitVectToText(self, arg1: ExplicitBitVect) -> str: ...
@overload
def BraunBlanquetSimilarity(
    self, bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0
) -> float:
    """
    B(bv1&bv2) / max(B(bv1),B(bv2))

    C++ signature :double BraunBlanquetSimilarity(ExplicitBitVect,ExplicitBitVect [,bool=0])

    BraunBlanquetSimilarity( (SparseBitVect)bv1, (str)pkl [, (bool)returnDistance=0]) -> float :

    C++ signature :double BraunBlanquetSimilarity(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

    BraunBlanquetSimilarity( (ExplicitBitVect)bv1, (str)pkl [, (bool)returnDistance=0]) -> float :B(bv1&bv2) / max(B(bv1),B(bv2))

    C++ signature :double BraunBlanquetSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])
    """
    ...

@overload
def BraunBlanquetSimilarity(
    self, bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = 0
) -> float:
    """


    C++ signature :double BraunBlanquetSimilarity(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

    BraunBlanquetSimilarity( (ExplicitBitVect)bv1, (str)pkl [, (bool)returnDistance=0]) -> float :B(bv1&bv2) / max(B(bv1),B(bv2))

    C++ signature :double BraunBlanquetSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])
    """
    ...

@overload
def BraunBlanquetSimilarity(
    self, bv1: SparseBitVect, pkl: str, returnDistance: bool = 0
) -> float:
    """
    B(bv1&bv2) / max(B(bv1),B(bv2))

    C++ signature :double BraunBlanquetSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])
    """
    ...

@overload
def BraunBlanquetSimilarity(
    self, bv1: ExplicitBitVect, pkl: str, returnDistance: bool = 0
) -> float: ...
def BraunBlanquetSimilarityNeighbors(
    self, bvqueries: AtomPairsParameters, bvList: AtomPairsParameters
) -> list:
    """
    B(bv1&bv2) / max(B(bv1),B(bv2))

    C++ signature :boost::python::list BraunBlanquetSimilarityNeighbors(boost::python::api::object,boost::python::api::object)
    """
    ...

def BraunBlanquetSimilarityNeighbors_sparse(
    self, bvqueries: AtomPairsParameters, bvList: AtomPairsParameters
) -> list:
    """
    B(bv1&bv2) / max(B(bv1),B(bv2))

    C++ signature :boost::python::list BraunBlanquetSimilarityNeighbors_sparse(boost::python::api::object,boost::python::api::object)
    """
    ...

@overload
def BulkAllBitSimilarity(
    self, v1: ExplicitBitVect, v2: AtomPairsParameters, returnDistance: bool = 0
) -> list:
    """
    (B(bv1) - B(bv1^bv2)) / B(bv1)

    C++ signature :boost::python::list BulkAllBitSimilarity(ExplicitBitVect const*,boost::python::api::object [,bool=0])
    """
    ...

@overload
def BulkAllBitSimilarity(
    self, v1: ExplicitBitVect, v2: AtomPairsParameters, returnDistance: bool = 0
) -> list: ...
@overload
def BulkAsymmetricSimilarity(
    self, bv1: SparseBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0
) -> list:
    """
    B(bv1&bv2) / min(B(bv1),B(bv2))

    C++ signature :boost::python::list BulkAsymmetricSimilarity(ExplicitBitVect const*,boost::python::api::object [,bool=0])
    """
    ...

@overload
def BulkAsymmetricSimilarity(
    self, bv1: ExplicitBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0
) -> list: ...
@overload
def BulkBraunBlanquetSimilarity(
    self, bv1: SparseBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0
) -> list:
    """
    B(bv1&bv2) / max(B(bv1),B(bv2))

    C++ signature :boost::python::list BulkBraunBlanquetSimilarity(ExplicitBitVect const*,boost::python::api::object [,bool=0])
    """
    ...

@overload
def BulkBraunBlanquetSimilarity(
    self, bv1: ExplicitBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0
) -> list: ...
@overload
def BulkCosineSimilarity(
    self, bv1: SparseBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0
) -> list:
    """
    B(bv1&bv2) / sqrt(B(bv1) * B(bv2))

    C++ signature :boost::python::list BulkCosineSimilarity(ExplicitBitVect const*,boost::python::api::object [,bool=0])
    """
    ...

@overload
def BulkCosineSimilarity(
    self, bv1: ExplicitBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0
) -> list: ...
@overload
def BulkDiceSimilarity(
    self, bv1: SparseBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0
) -> list:
    """
    2*B(bv1&bv2) / (B(bv1) + B(bv2))

    C++ signature :boost::python::list BulkDiceSimilarity(ExplicitBitVect const*,boost::python::api::object [,bool=0])

    BulkDiceSimilarity( (IntSparseIntVect)v1, (list)v2 [, (bool)returnDistance=False]) -> list :return the Dice similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkDiceSimilarity(RDKit::SparseIntVect<int>,boost::python::list [,bool=False])

    BulkDiceSimilarity( (LongSparseIntVect)v1, (list)v2 [, (bool)returnDistance=False]) -> list :return the Dice similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkDiceSimilarity(RDKit::SparseIntVect<long>,boost::python::list [,bool=False])

    BulkDiceSimilarity( (UIntSparseIntVect)v1, (list)v2 [, (bool)returnDistance=False]) -> list :return the Dice similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkDiceSimilarity(RDKit::SparseIntVect<unsigned int>,boost::python::list [,bool=False])

    BulkDiceSimilarity( (ULongSparseIntVect)v1, (list)v2 [, (bool)returnDistance=False]) -> list :return the Dice similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkDiceSimilarity(RDKit::SparseIntVect<unsigned long>,boost::python::list [,bool=False])
    """
    ...

@overload
def BulkDiceSimilarity(
    self, bv1: ExplicitBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0
) -> list:
    """
    return the Dice similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkDiceSimilarity(RDKit::SparseIntVect<int>,boost::python::list [,bool=False])

    BulkDiceSimilarity( (LongSparseIntVect)v1, (list)v2 [, (bool)returnDistance=False]) -> list :return the Dice similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkDiceSimilarity(RDKit::SparseIntVect<long>,boost::python::list [,bool=False])

    BulkDiceSimilarity( (UIntSparseIntVect)v1, (list)v2 [, (bool)returnDistance=False]) -> list :return the Dice similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkDiceSimilarity(RDKit::SparseIntVect<unsigned int>,boost::python::list [,bool=False])

    BulkDiceSimilarity( (ULongSparseIntVect)v1, (list)v2 [, (bool)returnDistance=False]) -> list :return the Dice similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkDiceSimilarity(RDKit::SparseIntVect<unsigned long>,boost::python::list [,bool=False])
    """
    ...

@overload
def BulkDiceSimilarity(
    self, v1: IntSparseIntVect, v2: list, returnDistance: bool = False
) -> list:
    """
    return the Dice similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkDiceSimilarity(RDKit::SparseIntVect<long>,boost::python::list [,bool=False])

    BulkDiceSimilarity( (UIntSparseIntVect)v1, (list)v2 [, (bool)returnDistance=False]) -> list :return the Dice similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkDiceSimilarity(RDKit::SparseIntVect<unsigned int>,boost::python::list [,bool=False])

    BulkDiceSimilarity( (ULongSparseIntVect)v1, (list)v2 [, (bool)returnDistance=False]) -> list :return the Dice similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkDiceSimilarity(RDKit::SparseIntVect<unsigned long>,boost::python::list [,bool=False])
    """
    ...

@overload
def BulkDiceSimilarity(
    self, v1: LongSparseIntVect, v2: list, returnDistance: bool = False
) -> list:
    """
    return the Dice similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkDiceSimilarity(RDKit::SparseIntVect<unsigned int>,boost::python::list [,bool=False])

    BulkDiceSimilarity( (ULongSparseIntVect)v1, (list)v2 [, (bool)returnDistance=False]) -> list :return the Dice similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkDiceSimilarity(RDKit::SparseIntVect<unsigned long>,boost::python::list [,bool=False])
    """
    ...

@overload
def BulkDiceSimilarity(
    self, v1: UIntSparseIntVect, v2: list, returnDistance: bool = False
) -> list:
    """
    return the Dice similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkDiceSimilarity(RDKit::SparseIntVect<unsigned long>,boost::python::list [,bool=False])
    """
    ...

@overload
def BulkDiceSimilarity(
    self, v1: ULongSparseIntVect, v2: list, returnDistance: bool = False
) -> list: ...
@overload
def BulkKulczynskiSimilarity(
    self, bv1: SparseBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0
) -> list:
    """
    B(bv1&bv2)*(B(bv1) + B(bv2)) / (2 * B(bv1) * B(bv2))

    C++ signature :boost::python::list BulkKulczynskiSimilarity(ExplicitBitVect const*,boost::python::api::object [,bool=0])
    """
    ...

@overload
def BulkKulczynskiSimilarity(
    self, bv1: ExplicitBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0
) -> list: ...
@overload
def BulkMcConnaugheySimilarity(
    self, bv1: SparseBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0
) -> list:
    """
    (B(bv1&bv2) * (B(bv1)+B(bv2)) - B(bv1)*B(bv2)) / (B(bv1) * B(bv2))

    C++ signature :boost::python::list BulkMcConnaugheySimilarity(ExplicitBitVect const*,boost::python::api::object [,bool=0])
    """
    ...

@overload
def BulkMcConnaugheySimilarity(
    self, bv1: ExplicitBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0
) -> list: ...
@overload
def BulkOnBitSimilarity(
    self, v1: ExplicitBitVect, v2: AtomPairsParameters, returnDistance: bool = 0
) -> list:
    """
    B(bv1&bv2) / B(bv1|bv2)

    C++ signature :boost::python::list BulkOnBitSimilarity(ExplicitBitVect const*,boost::python::api::object [,bool=0])
    """
    ...

@overload
def BulkOnBitSimilarity(
    self, v1: ExplicitBitVect, v2: AtomPairsParameters, returnDistance: bool = 0
) -> list: ...
@overload
def BulkRogotGoldbergSimilarity(
    self, bv1: SparseBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0
) -> list:
    """
    B(bv1&bv2) / B(bv1)

    C++ signature :boost::python::list BulkRogotGoldbergSimilarity(ExplicitBitVect const*,boost::python::api::object [,bool=0])
    """
    ...

@overload
def BulkRogotGoldbergSimilarity(
    self, bv1: ExplicitBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0
) -> list: ...
@overload
def BulkRusselSimilarity(
    self, bv1: SparseBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0
) -> list:
    """
    B(bv1&bv2) / B(bv1)

    C++ signature :boost::python::list BulkRusselSimilarity(ExplicitBitVect const*,boost::python::api::object [,bool=0])
    """
    ...

@overload
def BulkRusselSimilarity(
    self, bv1: ExplicitBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0
) -> list: ...
@overload
def BulkSokalSimilarity(
    self, bv1: SparseBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0
) -> list:
    """
    B(bv1&bv2) / (2*B(bv1) + 2*B(bv2) - 3*B(bv1&bv2))

    C++ signature :boost::python::list BulkSokalSimilarity(ExplicitBitVect const*,boost::python::api::object [,bool=0])
    """
    ...

@overload
def BulkSokalSimilarity(
    self, bv1: ExplicitBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0
) -> list: ...
@overload
def BulkTanimotoSimilarity(
    self, bv1: SparseBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0
) -> list:
    """
    B(bv1&bv2) / (B(bv1) + B(bv2) - B(bv1&bv2))

    C++ signature :boost::python::list BulkTanimotoSimilarity(ExplicitBitVect const*,boost::python::api::object [,bool=0])

    BulkTanimotoSimilarity( (IntSparseIntVect)v1, (list)v2 [, (bool)returnDistance=False]) -> list :return the Tanimoto similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkTanimotoSimilarity(RDKit::SparseIntVect<int>,boost::python::list [,bool=False])

    BulkTanimotoSimilarity( (LongSparseIntVect)v1, (list)v2 [, (bool)returnDistance=False]) -> list :return the Tanimoto similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkTanimotoSimilarity(RDKit::SparseIntVect<long>,boost::python::list [,bool=False])

    BulkTanimotoSimilarity( (UIntSparseIntVect)v1, (list)v2 [, (bool)returnDistance=False]) -> list :return the Tanimoto similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkTanimotoSimilarity(RDKit::SparseIntVect<unsigned int>,boost::python::list [,bool=False])

    BulkTanimotoSimilarity( (ULongSparseIntVect)v1, (list)v2 [, (bool)returnDistance=False]) -> list :return the Tanimoto similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkTanimotoSimilarity(RDKit::SparseIntVect<unsigned long>,boost::python::list [,bool=False])
    """
    ...

@overload
def BulkTanimotoSimilarity(
    self, bv1: ExplicitBitVect, bvList: AtomPairsParameters, returnDistance: bool = 0
) -> list:
    """
    return the Tanimoto similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkTanimotoSimilarity(RDKit::SparseIntVect<int>,boost::python::list [,bool=False])

    BulkTanimotoSimilarity( (LongSparseIntVect)v1, (list)v2 [, (bool)returnDistance=False]) -> list :return the Tanimoto similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkTanimotoSimilarity(RDKit::SparseIntVect<long>,boost::python::list [,bool=False])

    BulkTanimotoSimilarity( (UIntSparseIntVect)v1, (list)v2 [, (bool)returnDistance=False]) -> list :return the Tanimoto similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkTanimotoSimilarity(RDKit::SparseIntVect<unsigned int>,boost::python::list [,bool=False])

    BulkTanimotoSimilarity( (ULongSparseIntVect)v1, (list)v2 [, (bool)returnDistance=False]) -> list :return the Tanimoto similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkTanimotoSimilarity(RDKit::SparseIntVect<unsigned long>,boost::python::list [,bool=False])
    """
    ...

@overload
def BulkTanimotoSimilarity(
    self, v1: IntSparseIntVect, v2: list, returnDistance: bool = False
) -> list:
    """
    return the Tanimoto similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkTanimotoSimilarity(RDKit::SparseIntVect<long>,boost::python::list [,bool=False])

    BulkTanimotoSimilarity( (UIntSparseIntVect)v1, (list)v2 [, (bool)returnDistance=False]) -> list :return the Tanimoto similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkTanimotoSimilarity(RDKit::SparseIntVect<unsigned int>,boost::python::list [,bool=False])

    BulkTanimotoSimilarity( (ULongSparseIntVect)v1, (list)v2 [, (bool)returnDistance=False]) -> list :return the Tanimoto similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkTanimotoSimilarity(RDKit::SparseIntVect<unsigned long>,boost::python::list [,bool=False])
    """
    ...

@overload
def BulkTanimotoSimilarity(
    self, v1: LongSparseIntVect, v2: list, returnDistance: bool = False
) -> list:
    """
    return the Tanimoto similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkTanimotoSimilarity(RDKit::SparseIntVect<unsigned int>,boost::python::list [,bool=False])

    BulkTanimotoSimilarity( (ULongSparseIntVect)v1, (list)v2 [, (bool)returnDistance=False]) -> list :return the Tanimoto similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkTanimotoSimilarity(RDKit::SparseIntVect<unsigned long>,boost::python::list [,bool=False])
    """
    ...

@overload
def BulkTanimotoSimilarity(
    self, v1: UIntSparseIntVect, v2: list, returnDistance: bool = False
) -> list:
    """
    return the Tanimoto similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkTanimotoSimilarity(RDKit::SparseIntVect<unsigned long>,boost::python::list [,bool=False])
    """
    ...

@overload
def BulkTanimotoSimilarity(
    self, v1: ULongSparseIntVect, v2: list, returnDistance: bool = False
) -> list: ...
@overload
def BulkTverskySimilarity(
    self,
    bv1: SparseBitVect,
    bvList: AtomPairsParameters,
    a: float,
    b: float,
    returnDistance: bool = 0,
) -> list:
    """
    B(bv1&bv2) / (a*B(bv1)+b*B(bv2)+(1-a-b)*B(bv1&bv2)

    C++ signature :boost::python::list BulkTverskySimilarity(ExplicitBitVect const*,boost::python::api::object,double,double [,bool=0])

    BulkTverskySimilarity( (IntSparseIntVect)v1, (list)v2, (float)a, (float)b [, (bool)returnDistance=False]) -> list :return the Tversky similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkTverskySimilarity(RDKit::SparseIntVect<int>,boost::python::list,double,double [,bool=False])

    BulkTverskySimilarity( (LongSparseIntVect)v1, (list)v2, (float)a, (float)b [, (bool)returnDistance=False]) -> list :return the Tversky similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkTverskySimilarity(RDKit::SparseIntVect<long>,boost::python::list,double,double [,bool=False])

    BulkTverskySimilarity( (UIntSparseIntVect)v1, (list)v2, (float)a, (float)b [, (bool)returnDistance=False]) -> list :return the Tversky similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkTverskySimilarity(RDKit::SparseIntVect<unsigned int>,boost::python::list,double,double [,bool=False])

    BulkTverskySimilarity( (ULongSparseIntVect)v1, (list)v2, (float)a, (float)b [, (bool)returnDistance=False]) -> list :return the Tversky similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkTverskySimilarity(RDKit::SparseIntVect<unsigned long>,boost::python::list,double,double [,bool=False])
    """
    ...

@overload
def BulkTverskySimilarity(
    self,
    bv1: ExplicitBitVect,
    bvList: AtomPairsParameters,
    a: float,
    b: float,
    returnDistance: bool = 0,
) -> list:
    """
    return the Tversky similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkTverskySimilarity(RDKit::SparseIntVect<int>,boost::python::list,double,double [,bool=False])

    BulkTverskySimilarity( (LongSparseIntVect)v1, (list)v2, (float)a, (float)b [, (bool)returnDistance=False]) -> list :return the Tversky similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkTverskySimilarity(RDKit::SparseIntVect<long>,boost::python::list,double,double [,bool=False])

    BulkTverskySimilarity( (UIntSparseIntVect)v1, (list)v2, (float)a, (float)b [, (bool)returnDistance=False]) -> list :return the Tversky similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkTverskySimilarity(RDKit::SparseIntVect<unsigned int>,boost::python::list,double,double [,bool=False])

    BulkTverskySimilarity( (ULongSparseIntVect)v1, (list)v2, (float)a, (float)b [, (bool)returnDistance=False]) -> list :return the Tversky similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkTverskySimilarity(RDKit::SparseIntVect<unsigned long>,boost::python::list,double,double [,bool=False])
    """
    ...

@overload
def BulkTverskySimilarity(
    self,
    v1: IntSparseIntVect,
    v2: list,
    a: float,
    b: float,
    returnDistance: bool = False,
) -> list:
    """
    return the Tversky similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkTverskySimilarity(RDKit::SparseIntVect<long>,boost::python::list,double,double [,bool=False])

    BulkTverskySimilarity( (UIntSparseIntVect)v1, (list)v2, (float)a, (float)b [, (bool)returnDistance=False]) -> list :return the Tversky similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkTverskySimilarity(RDKit::SparseIntVect<unsigned int>,boost::python::list,double,double [,bool=False])

    BulkTverskySimilarity( (ULongSparseIntVect)v1, (list)v2, (float)a, (float)b [, (bool)returnDistance=False]) -> list :return the Tversky similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkTverskySimilarity(RDKit::SparseIntVect<unsigned long>,boost::python::list,double,double [,bool=False])
    """
    ...

@overload
def BulkTverskySimilarity(
    self,
    v1: LongSparseIntVect,
    v2: list,
    a: float,
    b: float,
    returnDistance: bool = False,
) -> list:
    """
    return the Tversky similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkTverskySimilarity(RDKit::SparseIntVect<unsigned int>,boost::python::list,double,double [,bool=False])

    BulkTverskySimilarity( (ULongSparseIntVect)v1, (list)v2, (float)a, (float)b [, (bool)returnDistance=False]) -> list :return the Tversky similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkTverskySimilarity(RDKit::SparseIntVect<unsigned long>,boost::python::list,double,double [,bool=False])
    """
    ...

@overload
def BulkTverskySimilarity(
    self,
    v1: UIntSparseIntVect,
    v2: list,
    a: float,
    b: float,
    returnDistance: bool = False,
) -> list:
    """
    return the Tversky similarities between one vector and a sequence of others

    C++ signature :boost::python::list BulkTverskySimilarity(RDKit::SparseIntVect<unsigned long>,boost::python::list,double,double [,bool=False])
    """
    ...

@overload
def BulkTverskySimilarity(
    self,
    v1: ULongSparseIntVect,
    v2: list,
    a: float,
    b: float,
    returnDistance: bool = False,
) -> list: ...
def ComputeL1Norm(self, arg1: DiscreteValueVect, arg2: DiscreteValueVect) -> int:
    """
    Compute the distance between two discrete vector values

    C++ signature :unsigned int ComputeL1Norm(RDKit::DiscreteValueVect,RDKit::DiscreteValueVect)
    """
    ...

def ConvertToExplicit(self, arg1: SparseBitVect) -> ExplicitBitVect:
    """
    Converts a SparseBitVector to an ExplicitBitVector and returns the ExplicitBitVector

    C++ signature :ExplicitBitVect* ConvertToExplicit(SparseBitVect const*)"""
    ...

@overload
def ConvertToNumpyArray(
    self, bv: ExplicitBitVect, destArray: AtomPairsParameters
) -> None:
    """


    C++ signature :void ConvertToNumpyArray(RDKit::DiscreteValueVect,boost::python::api::object)

    ConvertToNumpyArray( (IntSparseIntVect)bv, (AtomPairsParameters)destArray) -> None :

    C++ signature :void ConvertToNumpyArray(RDKit::SparseIntVect<int>,boost::python::api::object)

    ConvertToNumpyArray( (LongSparseIntVect)bv, (AtomPairsParameters)destArray) -> None :

    C++ signature :void ConvertToNumpyArray(RDKit::SparseIntVect<long>,boost::python::api::object)

    ConvertToNumpyArray( (UIntSparseIntVect)bv, (AtomPairsParameters)destArray) -> None :

    C++ signature :void ConvertToNumpyArray(RDKit::SparseIntVect<unsigned int>,boost::python::api::object)

    ConvertToNumpyArray( (ULongSparseIntVect)bv, (AtomPairsParameters)destArray) -> None :

    C++ signature :void ConvertToNumpyArray(RDKit::SparseIntVect<unsigned long>,boost::python::api::object)
    """
    ...

@overload
def ConvertToNumpyArray(
    self, bv: DiscreteValueVect, destArray: AtomPairsParameters
) -> None:
    """


    C++ signature :void ConvertToNumpyArray(RDKit::SparseIntVect<int>,boost::python::api::object)

    ConvertToNumpyArray( (LongSparseIntVect)bv, (AtomPairsParameters)destArray) -> None :

    C++ signature :void ConvertToNumpyArray(RDKit::SparseIntVect<long>,boost::python::api::object)

    ConvertToNumpyArray( (UIntSparseIntVect)bv, (AtomPairsParameters)destArray) -> None :

    C++ signature :void ConvertToNumpyArray(RDKit::SparseIntVect<unsigned int>,boost::python::api::object)

    ConvertToNumpyArray( (ULongSparseIntVect)bv, (AtomPairsParameters)destArray) -> None :

    C++ signature :void ConvertToNumpyArray(RDKit::SparseIntVect<unsigned long>,boost::python::api::object)
    """
    ...

@overload
def ConvertToNumpyArray(
    self, bv: IntSparseIntVect, destArray: AtomPairsParameters
) -> None:
    """


    C++ signature :void ConvertToNumpyArray(RDKit::SparseIntVect<long>,boost::python::api::object)

    ConvertToNumpyArray( (UIntSparseIntVect)bv, (AtomPairsParameters)destArray) -> None :

    C++ signature :void ConvertToNumpyArray(RDKit::SparseIntVect<unsigned int>,boost::python::api::object)

    ConvertToNumpyArray( (ULongSparseIntVect)bv, (AtomPairsParameters)destArray) -> None :

    C++ signature :void ConvertToNumpyArray(RDKit::SparseIntVect<unsigned long>,boost::python::api::object)
    """
    ...

@overload
def ConvertToNumpyArray(
    self, bv: LongSparseIntVect, destArray: AtomPairsParameters
) -> None:
    """


    C++ signature :void ConvertToNumpyArray(RDKit::SparseIntVect<unsigned int>,boost::python::api::object)

    ConvertToNumpyArray( (ULongSparseIntVect)bv, (AtomPairsParameters)destArray) -> None :

    C++ signature :void ConvertToNumpyArray(RDKit::SparseIntVect<unsigned long>,boost::python::api::object)
    """
    ...

@overload
def ConvertToNumpyArray(
    self, bv: UIntSparseIntVect, destArray: AtomPairsParameters
) -> None:
    """


    C++ signature :void ConvertToNumpyArray(RDKit::SparseIntVect<unsigned long>,boost::python::api::object)
    """
    ...

@overload
def ConvertToNumpyArray(
    self, bv: ULongSparseIntVect, destArray: AtomPairsParameters
) -> None: ...
@overload
def CosineSimilarity(
    self, bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0
) -> float:
    """
    B(bv1&bv2) / sqrt(B(bv1) * B(bv2))

    C++ signature :double CosineSimilarity(ExplicitBitVect,ExplicitBitVect [,bool=0])

    CosineSimilarity( (SparseBitVect)bv1, (str)pkl [, (bool)returnDistance=0]) -> float :

    C++ signature :double CosineSimilarity(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

    CosineSimilarity( (ExplicitBitVect)bv1, (str)pkl [, (bool)returnDistance=0]) -> float :B(bv1&bv2) / sqrt(B(bv1) * B(bv2))

    C++ signature :double CosineSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])
    """
    ...

@overload
def CosineSimilarity(
    self, bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = 0
) -> float:
    """


    C++ signature :double CosineSimilarity(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

    CosineSimilarity( (ExplicitBitVect)bv1, (str)pkl [, (bool)returnDistance=0]) -> float :B(bv1&bv2) / sqrt(B(bv1) * B(bv2))

    C++ signature :double CosineSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])
    """
    ...

@overload
def CosineSimilarity(
    self, bv1: SparseBitVect, pkl: str, returnDistance: bool = 0
) -> float:
    """
    B(bv1&bv2) / sqrt(B(bv1) * B(bv2))

    C++ signature :double CosineSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])
    """
    ...

@overload
def CosineSimilarity(
    self, bv1: ExplicitBitVect, pkl: str, returnDistance: bool = 0
) -> float: ...
def CosineSimilarityNeighbors(
    self, bvqueries: AtomPairsParameters, bvList: AtomPairsParameters
) -> list:
    """
    B(bv1&bv2) / sqrt(B(bv1) * B(bv2))

    C++ signature :boost::python::list CosineSimilarityNeighbors(boost::python::api::object,boost::python::api::object)
    """
    ...

def CosineSimilarityNeighbors_sparse(
    self, bvqueries: AtomPairsParameters, bvList: AtomPairsParameters
) -> list:
    """
    B(bv1&bv2) / sqrt(B(bv1) * B(bv2))

    C++ signature :boost::python::list CosineSimilarityNeighbors_sparse(boost::python::api::object,boost::python::api::object)
    """
    ...

def CreateFromBinaryText(self, arg1: str) -> ExplicitBitVect:
    """
    Creates an ExplicitBitVect from a binary string (byte array).

    C++ signature :ExplicitBitVect* CreateFromBinaryText(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
    ...

def CreateFromBitString(self, arg1: str) -> ExplicitBitVect:
    """
    Creates an ExplicitBitVect from a bit string (string of 0s and 1s).

    C++ signature :ExplicitBitVect* CreateFromBitString(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
    ...

def CreateFromFPSText(self, arg1: str) -> ExplicitBitVect:
    """
    Creates an ExplicitBitVect from an FPS string.

    C++ signature :ExplicitBitVect* CreateFromFPSText(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
    ...

@overload
def DiceSimilarity(
    self, bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0
) -> float:
    """
    2*B(bv1&bv2) / (B(bv1) + B(bv2))

    C++ signature :double DiceSimilarity(ExplicitBitVect,ExplicitBitVect [,bool=0])

    DiceSimilarity( (SparseBitVect)bv1, (str)pkl [, (bool)returnDistance=0]) -> float :

    C++ signature :double DiceSimilarity(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

    DiceSimilarity( (ExplicitBitVect)bv1, (str)pkl [, (bool)returnDistance=0]) -> float :2*B(bv1&bv2) / (B(bv1) + B(bv2))

    C++ signature :double DiceSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

    DiceSimilarity( (IntSparseIntVect)siv1, (IntSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Dice similarity between two vectors

    C++ signature :double DiceSimilarity(RDKit::SparseIntVect<int>,RDKit::SparseIntVect<int> [,bool=False [,double=0.0]])

    DiceSimilarity( (LongSparseIntVect)siv1, (LongSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Dice similarity between two vectors

    C++ signature :double DiceSimilarity(RDKit::SparseIntVect<long>,RDKit::SparseIntVect<long> [,bool=False [,double=0.0]])

    DiceSimilarity( (UIntSparseIntVect)siv1, (UIntSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Dice similarity between two vectors

    C++ signature :double DiceSimilarity(RDKit::SparseIntVect<unsigned int>,RDKit::SparseIntVect<unsigned int> [,bool=False [,double=0.0]])

    DiceSimilarity( (ULongSparseIntVect)siv1, (ULongSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Dice similarity between two vectors

    C++ signature :double DiceSimilarity(RDKit::SparseIntVect<unsigned long>,RDKit::SparseIntVect<unsigned long> [,bool=False [,double=0.0]])
    """
    ...

@overload
def DiceSimilarity(
    self, bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = 0
) -> float:
    """


    C++ signature :double DiceSimilarity(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

    DiceSimilarity( (ExplicitBitVect)bv1, (str)pkl [, (bool)returnDistance=0]) -> float :2*B(bv1&bv2) / (B(bv1) + B(bv2))

    C++ signature :double DiceSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

    DiceSimilarity( (IntSparseIntVect)siv1, (IntSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Dice similarity between two vectors

    C++ signature :double DiceSimilarity(RDKit::SparseIntVect<int>,RDKit::SparseIntVect<int> [,bool=False [,double=0.0]])

    DiceSimilarity( (LongSparseIntVect)siv1, (LongSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Dice similarity between two vectors

    C++ signature :double DiceSimilarity(RDKit::SparseIntVect<long>,RDKit::SparseIntVect<long> [,bool=False [,double=0.0]])

    DiceSimilarity( (UIntSparseIntVect)siv1, (UIntSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Dice similarity between two vectors

    C++ signature :double DiceSimilarity(RDKit::SparseIntVect<unsigned int>,RDKit::SparseIntVect<unsigned int> [,bool=False [,double=0.0]])

    DiceSimilarity( (ULongSparseIntVect)siv1, (ULongSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Dice similarity between two vectors

    C++ signature :double DiceSimilarity(RDKit::SparseIntVect<unsigned long>,RDKit::SparseIntVect<unsigned long> [,bool=False [,double=0.0]])
    """
    ...

@overload
def DiceSimilarity(
    self, bv1: SparseBitVect, pkl: str, returnDistance: bool = 0
) -> float:
    """
    2*B(bv1&bv2) / (B(bv1) + B(bv2))

    C++ signature :double DiceSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

    DiceSimilarity( (IntSparseIntVect)siv1, (IntSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Dice similarity between two vectors

    C++ signature :double DiceSimilarity(RDKit::SparseIntVect<int>,RDKit::SparseIntVect<int> [,bool=False [,double=0.0]])

    DiceSimilarity( (LongSparseIntVect)siv1, (LongSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Dice similarity between two vectors

    C++ signature :double DiceSimilarity(RDKit::SparseIntVect<long>,RDKit::SparseIntVect<long> [,bool=False [,double=0.0]])

    DiceSimilarity( (UIntSparseIntVect)siv1, (UIntSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Dice similarity between two vectors

    C++ signature :double DiceSimilarity(RDKit::SparseIntVect<unsigned int>,RDKit::SparseIntVect<unsigned int> [,bool=False [,double=0.0]])

    DiceSimilarity( (ULongSparseIntVect)siv1, (ULongSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Dice similarity between two vectors

    C++ signature :double DiceSimilarity(RDKit::SparseIntVect<unsigned long>,RDKit::SparseIntVect<unsigned long> [,bool=False [,double=0.0]])
    """
    ...

@overload
def DiceSimilarity(
    self, bv1: ExplicitBitVect, pkl: str, returnDistance: bool = 0
) -> float:
    """
    return the Dice similarity between two vectors

    C++ signature :double DiceSimilarity(RDKit::SparseIntVect<int>,RDKit::SparseIntVect<int> [,bool=False [,double=0.0]])

    DiceSimilarity( (LongSparseIntVect)siv1, (LongSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Dice similarity between two vectors

    C++ signature :double DiceSimilarity(RDKit::SparseIntVect<long>,RDKit::SparseIntVect<long> [,bool=False [,double=0.0]])

    DiceSimilarity( (UIntSparseIntVect)siv1, (UIntSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Dice similarity between two vectors

    C++ signature :double DiceSimilarity(RDKit::SparseIntVect<unsigned int>,RDKit::SparseIntVect<unsigned int> [,bool=False [,double=0.0]])

    DiceSimilarity( (ULongSparseIntVect)siv1, (ULongSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Dice similarity between two vectors

    C++ signature :double DiceSimilarity(RDKit::SparseIntVect<unsigned long>,RDKit::SparseIntVect<unsigned long> [,bool=False [,double=0.0]])
    """
    ...

@overload
def DiceSimilarity(
    self,
    siv1: IntSparseIntVect,
    siv2: IntSparseIntVect,
    returnDistance: bool = False,
    bounds: float = 0.0,
) -> float:
    """
    return the Dice similarity between two vectors

    C++ signature :double DiceSimilarity(RDKit::SparseIntVect<long>,RDKit::SparseIntVect<long> [,bool=False [,double=0.0]])

    DiceSimilarity( (UIntSparseIntVect)siv1, (UIntSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Dice similarity between two vectors

    C++ signature :double DiceSimilarity(RDKit::SparseIntVect<unsigned int>,RDKit::SparseIntVect<unsigned int> [,bool=False [,double=0.0]])

    DiceSimilarity( (ULongSparseIntVect)siv1, (ULongSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Dice similarity between two vectors

    C++ signature :double DiceSimilarity(RDKit::SparseIntVect<unsigned long>,RDKit::SparseIntVect<unsigned long> [,bool=False [,double=0.0]])
    """
    ...

@overload
def DiceSimilarity(
    self,
    siv1: LongSparseIntVect,
    siv2: LongSparseIntVect,
    returnDistance: bool = False,
    bounds: float = 0.0,
) -> float:
    """
    return the Dice similarity between two vectors

    C++ signature :double DiceSimilarity(RDKit::SparseIntVect<unsigned int>,RDKit::SparseIntVect<unsigned int> [,bool=False [,double=0.0]])

    DiceSimilarity( (ULongSparseIntVect)siv1, (ULongSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Dice similarity between two vectors

    C++ signature :double DiceSimilarity(RDKit::SparseIntVect<unsigned long>,RDKit::SparseIntVect<unsigned long> [,bool=False [,double=0.0]])
    """
    ...

@overload
def DiceSimilarity(
    self,
    siv1: UIntSparseIntVect,
    siv2: UIntSparseIntVect,
    returnDistance: bool = False,
    bounds: float = 0.0,
) -> float:
    """
    return the Dice similarity between two vectors

    C++ signature :double DiceSimilarity(RDKit::SparseIntVect<unsigned long>,RDKit::SparseIntVect<unsigned long> [,bool=False [,double=0.0]])
    """
    ...

@overload
def DiceSimilarity(
    self,
    siv1: ULongSparseIntVect,
    siv2: ULongSparseIntVect,
    returnDistance: bool = False,
    bounds: float = 0.0,
) -> float: ...
def DiceSimilarityNeighbors(
    self, bvqueries: AtomPairsParameters, bvList: AtomPairsParameters
) -> list:
    """
    2*B(bv1&bv2) / (B(bv1) + B(bv2))

    C++ signature :boost::python::list DiceSimilarityNeighbors(boost::python::api::object,boost::python::api::object)
    """
    ...

def DiceSimilarityNeighbors_sparse(
    self, bvqueries: AtomPairsParameters, bvList: AtomPairsParameters
) -> list:
    """
    2*B(bv1&bv2) / (B(bv1) + B(bv2))

    C++ signature :boost::python::list DiceSimilarityNeighbors_sparse(boost::python::api::object,boost::python::api::object)
    """
    ...

@overload
def FoldFingerprint(self, bv: SparseBitVect, foldFactor: int = 2) -> SparseBitVect:
    """
    Folds the fingerprint by the provided amount. The default, foldFactor=2, returns a fingerprint that is half the size of the original.

    C++ signature :ExplicitBitVect* FoldFingerprint(ExplicitBitVect [,unsigned int=2])
    """
    ...

@overload
def FoldFingerprint(
    self, bv: ExplicitBitVect, foldFactor: int = 2
) -> ExplicitBitVect: ...
@overload
def InitFromDaylightString(self, arg1: SparseBitVect, arg2: str) -> None:
    """
    Fill a BitVect using an ASCII (Daylight) encoding of a fingerprint.

    Arguments
    bv: either a _SparseBitVect_ or an _ExplicitBitVect_

    txt: a string with the Daylight encoding (this is the text thatthe Daylight tools put in the FP field of a TDT)

    C++ signature :void InitFromDaylightString(ExplicitBitVect {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
    ...

@overload
def InitFromDaylightString(self, arg1: ExplicitBitVect, arg2: str) -> None: ...
@overload
def KulczynskiSimilarity(
    self, bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0
) -> float:
    """
    B(bv1&bv2)*(B(bv1) + B(bv2)) / (2 * B(bv1) * B(bv2))

    C++ signature :double KulczynskiSimilarity(ExplicitBitVect,ExplicitBitVect [,bool=0])

    KulczynskiSimilarity( (SparseBitVect)bv1, (str)pkl [, (bool)returnDistance=0]) -> float :

    C++ signature :double KulczynskiSimilarity(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

    KulczynskiSimilarity( (ExplicitBitVect)bv1, (str)pkl [, (bool)returnDistance=0]) -> float :B(bv1&bv2)*(B(bv1) + B(bv2)) / (2 * B(bv1) * B(bv2))

    C++ signature :double KulczynskiSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])
    """
    ...

@overload
def KulczynskiSimilarity(
    self, bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = 0
) -> float:
    """


    C++ signature :double KulczynskiSimilarity(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

    KulczynskiSimilarity( (ExplicitBitVect)bv1, (str)pkl [, (bool)returnDistance=0]) -> float :B(bv1&bv2)*(B(bv1) + B(bv2)) / (2 * B(bv1) * B(bv2))

    C++ signature :double KulczynskiSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])
    """
    ...

@overload
def KulczynskiSimilarity(
    self, bv1: SparseBitVect, pkl: str, returnDistance: bool = 0
) -> float:
    """
    B(bv1&bv2)*(B(bv1) + B(bv2)) / (2 * B(bv1) * B(bv2))

    C++ signature :double KulczynskiSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])
    """
    ...

@overload
def KulczynskiSimilarity(
    self, bv1: ExplicitBitVect, pkl: str, returnDistance: bool = 0
) -> float: ...
def KulczynskiSimilarityNeighbors(
    self, bvqueries: AtomPairsParameters, bvList: AtomPairsParameters
) -> list:
    """
    B(bv1&bv2)*(B(bv1) + B(bv2)) / (2 * B(bv1) * B(bv2))

    C++ signature :boost::python::list KulczynskiSimilarityNeighbors(boost::python::api::object,boost::python::api::object)
    """
    ...

def KulczynskiSimilarityNeighbors_sparse(
    self, bvqueries: AtomPairsParameters, bvList: AtomPairsParameters
) -> list:
    """
    B(bv1&bv2)*(B(bv1) + B(bv2)) / (2 * B(bv1) * B(bv2))

    C++ signature :boost::python::list KulczynskiSimilarityNeighbors_sparse(boost::python::api::object,boost::python::api::object)
    """
    ...

@overload
def McConnaugheySimilarity(
    self, bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0
) -> float:
    """
    (B(bv1&bv2) * (B(bv1)+B(bv2)) - B(bv1)*B(bv2)) / (B(bv1) * B(bv2))

    C++ signature :double McConnaugheySimilarity(ExplicitBitVect,ExplicitBitVect [,bool=0])

    McConnaugheySimilarity( (SparseBitVect)bv1, (str)pkl [, (bool)returnDistance=0]) -> float :

    C++ signature :double McConnaugheySimilarity(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

    McConnaugheySimilarity( (ExplicitBitVect)bv1, (str)pkl [, (bool)returnDistance=0]) -> float :(B(bv1&bv2) * (B(bv1)+B(bv2)) - B(bv1)*B(bv2)) / (B(bv1) * B(bv2))

    C++ signature :double McConnaugheySimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])
    """
    ...

@overload
def McConnaugheySimilarity(
    self, bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = 0
) -> float:
    """


    C++ signature :double McConnaugheySimilarity(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

    McConnaugheySimilarity( (ExplicitBitVect)bv1, (str)pkl [, (bool)returnDistance=0]) -> float :(B(bv1&bv2) * (B(bv1)+B(bv2)) - B(bv1)*B(bv2)) / (B(bv1) * B(bv2))

    C++ signature :double McConnaugheySimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])
    """
    ...

@overload
def McConnaugheySimilarity(
    self, bv1: SparseBitVect, pkl: str, returnDistance: bool = 0
) -> float:
    """
    (B(bv1&bv2) * (B(bv1)+B(bv2)) - B(bv1)*B(bv2)) / (B(bv1) * B(bv2))

    C++ signature :double McConnaugheySimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])
    """
    ...

@overload
def McConnaugheySimilarity(
    self, bv1: ExplicitBitVect, pkl: str, returnDistance: bool = 0
) -> float: ...
def McConnaugheySimilarityNeighbors(
    self, bvqueries: AtomPairsParameters, bvList: AtomPairsParameters
) -> list:
    """
    (B(bv1&bv2) * (B(bv1)+B(bv2)) - B(bv1)*B(bv2)) / (B(bv1) * B(bv2))

    C++ signature :boost::python::list McConnaugheySimilarityNeighbors(boost::python::api::object,boost::python::api::object)
    """
    ...

def McConnaugheySimilarityNeighbors_sparse(
    self, bvqueries: AtomPairsParameters, bvList: AtomPairsParameters
) -> list:
    """
    (B(bv1&bv2) * (B(bv1)+B(bv2)) - B(bv1)*B(bv2)) / (B(bv1) * B(bv2))

    C++ signature :boost::python::list McConnaugheySimilarityNeighbors_sparse(boost::python::api::object,boost::python::api::object)
    """
    ...

@overload
def NumBitsInCommon(self, arg1: SparseBitVect, arg2: SparseBitVect) -> int:
    """
    Returns the total number of bits in common between the two bit vectors

    C++ signature :int NumBitsInCommon(ExplicitBitVect,ExplicitBitVect)"""
    ...

@overload
def NumBitsInCommon(self, arg1: ExplicitBitVect, arg2: ExplicitBitVect) -> int: ...
@overload
def OffBitProjSimilarity(self, arg1: SparseBitVect, arg2: SparseBitVect) -> _vectd:
    """


    C++ signature :std::vector<double, std::allocator<double> > OffBitProjSimilarity(ExplicitBitVect,ExplicitBitVect)
    """
    ...

@overload
def OffBitProjSimilarity(
    self, arg1: ExplicitBitVect, arg2: ExplicitBitVect
) -> _vectd: ...
@overload
def OffBitsInCommon(self, arg1: SparseBitVect, arg2: SparseBitVect) -> _vecti:
    """
    Returns the number of off bits in common between the two bit vectors

    C++ signature :std::vector<int, std::allocator<int> > OffBitsInCommon(ExplicitBitVect,ExplicitBitVect)
    """
    ...

@overload
def OffBitsInCommon(self, arg1: ExplicitBitVect, arg2: ExplicitBitVect) -> _vecti: ...
@overload
def OnBitProjSimilarity(self, arg1: SparseBitVect, arg2: SparseBitVect) -> _vectd:
    """
    Returns a 2-tuple: (B(bv1&bv2) / B(bv1), B(bv1&bv2) / B(bv2))

    C++ signature :std::vector<double, std::allocator<double> > OnBitProjSimilarity(ExplicitBitVect,ExplicitBitVect)
    """
    ...

@overload
def OnBitProjSimilarity(
    self, arg1: ExplicitBitVect, arg2: ExplicitBitVect
) -> _vectd: ...
@overload
def OnBitSimilarity(self, v1: SparseBitVect, v2: SparseBitVect) -> float:
    """
    B(bv1&bv2) / B(bv1|bv2)

    C++ signature :double OnBitSimilarity(ExplicitBitVect,ExplicitBitVect)"""
    ...

@overload
def OnBitSimilarity(self, v1: ExplicitBitVect, v2: ExplicitBitVect) -> float: ...
@overload
def OnBitsInCommon(self, arg1: SparseBitVect, arg2: SparseBitVect) -> _vecti:
    """
    Returns the number of on bits in common between the two bit vectors

    C++ signature :std::vector<int, std::allocator<int> > OnBitsInCommon(ExplicitBitVect,ExplicitBitVect)
    """
    ...

@overload
def OnBitsInCommon(self, arg1: ExplicitBitVect, arg2: ExplicitBitVect) -> _vecti: ...
@overload
def RogotGoldbergSimilarity(
    self, bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0
) -> float:
    """
    B(bv1&bv2) / B(bv1)

    C++ signature :double RogotGoldbergSimilarity(ExplicitBitVect,ExplicitBitVect [,bool=0])

    RogotGoldbergSimilarity( (SparseBitVect)bv1, (str)pkl [, (bool)returnDistance=0]) -> float :

    C++ signature :double RogotGoldbergSimilarity(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

    RogotGoldbergSimilarity( (ExplicitBitVect)bv1, (str)pkl [, (bool)returnDistance=0]) -> float :B(bv1&bv2) / B(bv1)

    C++ signature :double RogotGoldbergSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])
    """
    ...

@overload
def RogotGoldbergSimilarity(
    self, bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = 0
) -> float:
    """


    C++ signature :double RogotGoldbergSimilarity(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

    RogotGoldbergSimilarity( (ExplicitBitVect)bv1, (str)pkl [, (bool)returnDistance=0]) -> float :B(bv1&bv2) / B(bv1)

    C++ signature :double RogotGoldbergSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])
    """
    ...

@overload
def RogotGoldbergSimilarity(
    self, bv1: SparseBitVect, pkl: str, returnDistance: bool = 0
) -> float:
    """
    B(bv1&bv2) / B(bv1)

    C++ signature :double RogotGoldbergSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])
    """
    ...

@overload
def RogotGoldbergSimilarity(
    self, bv1: ExplicitBitVect, pkl: str, returnDistance: bool = 0
) -> float: ...
def RogotGoldbergSimilarityNeighbors(
    self, bvqueries: AtomPairsParameters, bvList: AtomPairsParameters
) -> list:
    """
    B(bv1&bv2) / B(bv1)

    C++ signature :boost::python::list RogotGoldbergSimilarityNeighbors(boost::python::api::object,boost::python::api::object)
    """
    ...

def RogotGoldbergSimilarityNeighbors_sparse(
    self, bvqueries: AtomPairsParameters, bvList: AtomPairsParameters
) -> list:
    """
    B(bv1&bv2) / B(bv1)

    C++ signature :boost::python::list RogotGoldbergSimilarityNeighbors_sparse(boost::python::api::object,boost::python::api::object)
    """
    ...

@overload
def RusselSimilarity(
    self, bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0
) -> float:
    """
    B(bv1&bv2) / B(bv1)

    C++ signature :double RusselSimilarity(ExplicitBitVect,ExplicitBitVect [,bool=0])

    RusselSimilarity( (SparseBitVect)bv1, (str)pkl [, (bool)returnDistance=0]) -> float :

    C++ signature :double RusselSimilarity(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

    RusselSimilarity( (ExplicitBitVect)bv1, (str)pkl [, (bool)returnDistance=0]) -> float :B(bv1&bv2) / B(bv1)

    C++ signature :double RusselSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])
    """
    ...

@overload
def RusselSimilarity(
    self, bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = 0
) -> float:
    """


    C++ signature :double RusselSimilarity(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

    RusselSimilarity( (ExplicitBitVect)bv1, (str)pkl [, (bool)returnDistance=0]) -> float :B(bv1&bv2) / B(bv1)

    C++ signature :double RusselSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])
    """
    ...

@overload
def RusselSimilarity(
    self, bv1: SparseBitVect, pkl: str, returnDistance: bool = 0
) -> float:
    """
    B(bv1&bv2) / B(bv1)

    C++ signature :double RusselSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])
    """
    ...

@overload
def RusselSimilarity(
    self, bv1: ExplicitBitVect, pkl: str, returnDistance: bool = 0
) -> float: ...
def RusselSimilarityNeighbors(
    self, bvqueries: AtomPairsParameters, bvList: AtomPairsParameters
) -> list:
    """
    B(bv1&bv2) / B(bv1)

    C++ signature :boost::python::list RusselSimilarityNeighbors(boost::python::api::object,boost::python::api::object)
    """
    ...

def RusselSimilarityNeighbors_sparse(
    self, bvqueries: AtomPairsParameters, bvList: AtomPairsParameters
) -> list:
    """
    B(bv1&bv2) / B(bv1)

    C++ signature :boost::python::list RusselSimilarityNeighbors_sparse(boost::python::api::object,boost::python::api::object)
    """
    ...

@overload
def SokalSimilarity(
    self, bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0
) -> float:
    """
    B(bv1&bv2) / (2*B(bv1) + 2*B(bv2) - 3*B(bv1&bv2))

    C++ signature :double SokalSimilarity(ExplicitBitVect,ExplicitBitVect [,bool=0])

    SokalSimilarity( (SparseBitVect)bv1, (str)pkl [, (bool)returnDistance=0]) -> float :

    C++ signature :double SokalSimilarity(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

    SokalSimilarity( (ExplicitBitVect)bv1, (str)pkl [, (bool)returnDistance=0]) -> float :B(bv1&bv2) / (2*B(bv1) + 2*B(bv2) - 3*B(bv1&bv2))

    C++ signature :double SokalSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])
    """
    ...

@overload
def SokalSimilarity(
    self, bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = 0
) -> float:
    """


    C++ signature :double SokalSimilarity(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

    SokalSimilarity( (ExplicitBitVect)bv1, (str)pkl [, (bool)returnDistance=0]) -> float :B(bv1&bv2) / (2*B(bv1) + 2*B(bv2) - 3*B(bv1&bv2))

    C++ signature :double SokalSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])
    """
    ...

@overload
def SokalSimilarity(
    self, bv1: SparseBitVect, pkl: str, returnDistance: bool = 0
) -> float:
    """
    B(bv1&bv2) / (2*B(bv1) + 2*B(bv2) - 3*B(bv1&bv2))

    C++ signature :double SokalSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])
    """
    ...

@overload
def SokalSimilarity(
    self, bv1: ExplicitBitVect, pkl: str, returnDistance: bool = 0
) -> float: ...
def SokalSimilarityNeighbors(
    self, bvqueries: AtomPairsParameters, bvList: AtomPairsParameters
) -> list:
    """
    B(bv1&bv2) / (2*B(bv1) + 2*B(bv2) - 3*B(bv1&bv2))

    C++ signature :boost::python::list SokalSimilarityNeighbors(boost::python::api::object,boost::python::api::object)
    """
    ...

def SokalSimilarityNeighbors_sparse(
    self, bvqueries: AtomPairsParameters, bvList: AtomPairsParameters
) -> list:
    """
    B(bv1&bv2) / (2*B(bv1) + 2*B(bv2) - 3*B(bv1&bv2))

    C++ signature :boost::python::list SokalSimilarityNeighbors_sparse(boost::python::api::object,boost::python::api::object)
    """
    ...

@overload
def TanimotoSimilarity(
    self, bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0
) -> float:
    """
    B(bv1&bv2) / (B(bv1) + B(bv2) - B(bv1&bv2))

    C++ signature :double TanimotoSimilarity(ExplicitBitVect,ExplicitBitVect [,bool=0])

    TanimotoSimilarity( (SparseBitVect)bv1, (str)pkl [, (bool)returnDistance=0]) -> float :

    C++ signature :double TanimotoSimilarity(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

    TanimotoSimilarity( (ExplicitBitVect)bv1, (str)pkl [, (bool)returnDistance=0]) -> float :B(bv1&bv2) / (B(bv1) + B(bv2) - B(bv1&bv2))

    C++ signature :double TanimotoSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

    TanimotoSimilarity( (IntSparseIntVect)siv1, (IntSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tanimoto similarity between two vectors

    C++ signature :double TanimotoSimilarity(RDKit::SparseIntVect<int>,RDKit::SparseIntVect<int> [,bool=False [,double=0.0]])

    TanimotoSimilarity( (LongSparseIntVect)siv1, (LongSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tanimoto similarity between two vectors

    C++ signature :double TanimotoSimilarity(RDKit::SparseIntVect<long>,RDKit::SparseIntVect<long> [,bool=False [,double=0.0]])

    TanimotoSimilarity( (UIntSparseIntVect)siv1, (UIntSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tanimoto similarity between two vectors

    C++ signature :double TanimotoSimilarity(RDKit::SparseIntVect<unsigned int>,RDKit::SparseIntVect<unsigned int> [,bool=False [,double=0.0]])

    TanimotoSimilarity( (ULongSparseIntVect)siv1, (ULongSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tanimoto similarity between two vectors

    C++ signature :double TanimotoSimilarity(RDKit::SparseIntVect<unsigned long>,RDKit::SparseIntVect<unsigned long> [,bool=False [,double=0.0]])
    """
    ...

@overload
def TanimotoSimilarity(
    self, bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = 0
) -> float:
    """


    C++ signature :double TanimotoSimilarity(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

    TanimotoSimilarity( (ExplicitBitVect)bv1, (str)pkl [, (bool)returnDistance=0]) -> float :B(bv1&bv2) / (B(bv1) + B(bv2) - B(bv1&bv2))

    C++ signature :double TanimotoSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

    TanimotoSimilarity( (IntSparseIntVect)siv1, (IntSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tanimoto similarity between two vectors

    C++ signature :double TanimotoSimilarity(RDKit::SparseIntVect<int>,RDKit::SparseIntVect<int> [,bool=False [,double=0.0]])

    TanimotoSimilarity( (LongSparseIntVect)siv1, (LongSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tanimoto similarity between two vectors

    C++ signature :double TanimotoSimilarity(RDKit::SparseIntVect<long>,RDKit::SparseIntVect<long> [,bool=False [,double=0.0]])

    TanimotoSimilarity( (UIntSparseIntVect)siv1, (UIntSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tanimoto similarity between two vectors

    C++ signature :double TanimotoSimilarity(RDKit::SparseIntVect<unsigned int>,RDKit::SparseIntVect<unsigned int> [,bool=False [,double=0.0]])

    TanimotoSimilarity( (ULongSparseIntVect)siv1, (ULongSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tanimoto similarity between two vectors

    C++ signature :double TanimotoSimilarity(RDKit::SparseIntVect<unsigned long>,RDKit::SparseIntVect<unsigned long> [,bool=False [,double=0.0]])
    """
    ...

@overload
def TanimotoSimilarity(
    self, bv1: SparseBitVect, pkl: str, returnDistance: bool = 0
) -> float:
    """
    B(bv1&bv2) / (B(bv1) + B(bv2) - B(bv1&bv2))

    C++ signature :double TanimotoSimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=0])

    TanimotoSimilarity( (IntSparseIntVect)siv1, (IntSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tanimoto similarity between two vectors

    C++ signature :double TanimotoSimilarity(RDKit::SparseIntVect<int>,RDKit::SparseIntVect<int> [,bool=False [,double=0.0]])

    TanimotoSimilarity( (LongSparseIntVect)siv1, (LongSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tanimoto similarity between two vectors

    C++ signature :double TanimotoSimilarity(RDKit::SparseIntVect<long>,RDKit::SparseIntVect<long> [,bool=False [,double=0.0]])

    TanimotoSimilarity( (UIntSparseIntVect)siv1, (UIntSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tanimoto similarity between two vectors

    C++ signature :double TanimotoSimilarity(RDKit::SparseIntVect<unsigned int>,RDKit::SparseIntVect<unsigned int> [,bool=False [,double=0.0]])

    TanimotoSimilarity( (ULongSparseIntVect)siv1, (ULongSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tanimoto similarity between two vectors

    C++ signature :double TanimotoSimilarity(RDKit::SparseIntVect<unsigned long>,RDKit::SparseIntVect<unsigned long> [,bool=False [,double=0.0]])
    """
    ...

@overload
def TanimotoSimilarity(
    self, bv1: ExplicitBitVect, pkl: str, returnDistance: bool = 0
) -> float:
    """
    return the Tanimoto similarity between two vectors

    C++ signature :double TanimotoSimilarity(RDKit::SparseIntVect<int>,RDKit::SparseIntVect<int> [,bool=False [,double=0.0]])

    TanimotoSimilarity( (LongSparseIntVect)siv1, (LongSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tanimoto similarity between two vectors

    C++ signature :double TanimotoSimilarity(RDKit::SparseIntVect<long>,RDKit::SparseIntVect<long> [,bool=False [,double=0.0]])

    TanimotoSimilarity( (UIntSparseIntVect)siv1, (UIntSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tanimoto similarity between two vectors

    C++ signature :double TanimotoSimilarity(RDKit::SparseIntVect<unsigned int>,RDKit::SparseIntVect<unsigned int> [,bool=False [,double=0.0]])

    TanimotoSimilarity( (ULongSparseIntVect)siv1, (ULongSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tanimoto similarity between two vectors

    C++ signature :double TanimotoSimilarity(RDKit::SparseIntVect<unsigned long>,RDKit::SparseIntVect<unsigned long> [,bool=False [,double=0.0]])
    """
    ...

@overload
def TanimotoSimilarity(
    self,
    siv1: IntSparseIntVect,
    siv2: IntSparseIntVect,
    returnDistance: bool = False,
    bounds: float = 0.0,
) -> float:
    """
    return the Tanimoto similarity between two vectors

    C++ signature :double TanimotoSimilarity(RDKit::SparseIntVect<long>,RDKit::SparseIntVect<long> [,bool=False [,double=0.0]])

    TanimotoSimilarity( (UIntSparseIntVect)siv1, (UIntSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tanimoto similarity between two vectors

    C++ signature :double TanimotoSimilarity(RDKit::SparseIntVect<unsigned int>,RDKit::SparseIntVect<unsigned int> [,bool=False [,double=0.0]])

    TanimotoSimilarity( (ULongSparseIntVect)siv1, (ULongSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tanimoto similarity between two vectors

    C++ signature :double TanimotoSimilarity(RDKit::SparseIntVect<unsigned long>,RDKit::SparseIntVect<unsigned long> [,bool=False [,double=0.0]])
    """
    ...

@overload
def TanimotoSimilarity(
    self,
    siv1: LongSparseIntVect,
    siv2: LongSparseIntVect,
    returnDistance: bool = False,
    bounds: float = 0.0,
) -> float:
    """
    return the Tanimoto similarity between two vectors

    C++ signature :double TanimotoSimilarity(RDKit::SparseIntVect<unsigned int>,RDKit::SparseIntVect<unsigned int> [,bool=False [,double=0.0]])

    TanimotoSimilarity( (ULongSparseIntVect)siv1, (ULongSparseIntVect)siv2 [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tanimoto similarity between two vectors

    C++ signature :double TanimotoSimilarity(RDKit::SparseIntVect<unsigned long>,RDKit::SparseIntVect<unsigned long> [,bool=False [,double=0.0]])
    """
    ...

@overload
def TanimotoSimilarity(
    self,
    siv1: UIntSparseIntVect,
    siv2: UIntSparseIntVect,
    returnDistance: bool = False,
    bounds: float = 0.0,
) -> float:
    """
    return the Tanimoto similarity between two vectors

    C++ signature :double TanimotoSimilarity(RDKit::SparseIntVect<unsigned long>,RDKit::SparseIntVect<unsigned long> [,bool=False [,double=0.0]])
    """
    ...

@overload
def TanimotoSimilarity(
    self,
    siv1: ULongSparseIntVect,
    siv2: ULongSparseIntVect,
    returnDistance: bool = False,
    bounds: float = 0.0,
) -> float: ...
def TanimotoSimilarityNeighbors(
    self, bvqueries: AtomPairsParameters, bvList: AtomPairsParameters
) -> list:
    """
    B(bv1&bv2) / (B(bv1) + B(bv2) - B(bv1&bv2))

    C++ signature :boost::python::list TanimotoSimilarityNeighbors(boost::python::api::object,boost::python::api::object)
    """
    ...

def TanimotoSimilarityNeighbors_sparse(
    self, bvqueries: AtomPairsParameters, bvList: AtomPairsParameters
) -> list:
    """
    B(bv1&bv2) / (B(bv1) + B(bv2) - B(bv1&bv2))

    C++ signature :boost::python::list TanimotoSimilarityNeighbors_sparse(boost::python::api::object,boost::python::api::object)
    """
    ...

@overload
def TverskySimilarity(
    self,
    bv1: SparseBitVect,
    bv2: SparseBitVect,
    a: float,
    b: float,
    returnDistance: bool = 0,
) -> float:
    """
    B(bv1&bv2) / (a*B(bv1)+b*B(bv2)+(1-a-b)*B(bv1&bv2)

    C++ signature :double TverskySimilarity(ExplicitBitVect,ExplicitBitVect,double,double [,bool=0])

    TverskySimilarity( (SparseBitVect)bv1, (str)pkl, (float)a, (float)b [, (bool)returnDistance=0]) -> float :

    C++ signature :double TverskySimilarity(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,double,double [,bool=0])

    TverskySimilarity( (ExplicitBitVect)bv1, (str)pkl, (float)a, (float)b [, (bool)returnDistance=0]) -> float :B(bv1&bv2) / (a*B(bv1)+b*B(bv2)+(1-a-b)*B(bv1&bv2)

    C++ signature :double TverskySimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,double,double [,bool=0])

    TverskySimilarity( (IntSparseIntVect)siv1, (IntSparseIntVect)siv2, (float)a, (float)b [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tversky similarity between two vectors

    C++ signature :double TverskySimilarity(RDKit::SparseIntVect<int>,RDKit::SparseIntVect<int>,double,double [,bool=False [,double=0.0]])

    TverskySimilarity( (LongSparseIntVect)siv1, (LongSparseIntVect)siv2, (float)a, (float)b [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tversky similarity between two vectors

    C++ signature :double TverskySimilarity(RDKit::SparseIntVect<long>,RDKit::SparseIntVect<long>,double,double [,bool=False [,double=0.0]])

    TverskySimilarity( (UIntSparseIntVect)siv1, (UIntSparseIntVect)siv2, (float)a, (float)b [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tversky similarity between two vectors

    C++ signature :double TverskySimilarity(RDKit::SparseIntVect<unsigned int>,RDKit::SparseIntVect<unsigned int>,double,double [,bool=False [,double=0.0]])

    TverskySimilarity( (ULongSparseIntVect)siv1, (ULongSparseIntVect)siv2, (float)a, (float)b [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tversky similarity between two vectors

    C++ signature :double TverskySimilarity(RDKit::SparseIntVect<unsigned long>,RDKit::SparseIntVect<unsigned long>,double,double [,bool=False [,double=0.0]])
    """
    ...

@overload
def TverskySimilarity(
    self,
    bv1: ExplicitBitVect,
    bv2: ExplicitBitVect,
    a: float,
    b: float,
    returnDistance: bool = 0,
) -> float:
    """


    C++ signature :double TverskySimilarity(SparseBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,double,double [,bool=0])

    TverskySimilarity( (ExplicitBitVect)bv1, (str)pkl, (float)a, (float)b [, (bool)returnDistance=0]) -> float :B(bv1&bv2) / (a*B(bv1)+b*B(bv2)+(1-a-b)*B(bv1&bv2)

    C++ signature :double TverskySimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,double,double [,bool=0])

    TverskySimilarity( (IntSparseIntVect)siv1, (IntSparseIntVect)siv2, (float)a, (float)b [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tversky similarity between two vectors

    C++ signature :double TverskySimilarity(RDKit::SparseIntVect<int>,RDKit::SparseIntVect<int>,double,double [,bool=False [,double=0.0]])

    TverskySimilarity( (LongSparseIntVect)siv1, (LongSparseIntVect)siv2, (float)a, (float)b [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tversky similarity between two vectors

    C++ signature :double TverskySimilarity(RDKit::SparseIntVect<long>,RDKit::SparseIntVect<long>,double,double [,bool=False [,double=0.0]])

    TverskySimilarity( (UIntSparseIntVect)siv1, (UIntSparseIntVect)siv2, (float)a, (float)b [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tversky similarity between two vectors

    C++ signature :double TverskySimilarity(RDKit::SparseIntVect<unsigned int>,RDKit::SparseIntVect<unsigned int>,double,double [,bool=False [,double=0.0]])

    TverskySimilarity( (ULongSparseIntVect)siv1, (ULongSparseIntVect)siv2, (float)a, (float)b [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tversky similarity between two vectors

    C++ signature :double TverskySimilarity(RDKit::SparseIntVect<unsigned long>,RDKit::SparseIntVect<unsigned long>,double,double [,bool=False [,double=0.0]])
    """
    ...

@overload
def TverskySimilarity(
    self, bv1: SparseBitVect, pkl: str, a: float, b: float, returnDistance: bool = 0
) -> float:
    """
    B(bv1&bv2) / (a*B(bv1)+b*B(bv2)+(1-a-b)*B(bv1&bv2)

    C++ signature :double TverskySimilarity(ExplicitBitVect,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,double,double [,bool=0])

    TverskySimilarity( (IntSparseIntVect)siv1, (IntSparseIntVect)siv2, (float)a, (float)b [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tversky similarity between two vectors

    C++ signature :double TverskySimilarity(RDKit::SparseIntVect<int>,RDKit::SparseIntVect<int>,double,double [,bool=False [,double=0.0]])

    TverskySimilarity( (LongSparseIntVect)siv1, (LongSparseIntVect)siv2, (float)a, (float)b [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tversky similarity between two vectors

    C++ signature :double TverskySimilarity(RDKit::SparseIntVect<long>,RDKit::SparseIntVect<long>,double,double [,bool=False [,double=0.0]])

    TverskySimilarity( (UIntSparseIntVect)siv1, (UIntSparseIntVect)siv2, (float)a, (float)b [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tversky similarity between two vectors

    C++ signature :double TverskySimilarity(RDKit::SparseIntVect<unsigned int>,RDKit::SparseIntVect<unsigned int>,double,double [,bool=False [,double=0.0]])

    TverskySimilarity( (ULongSparseIntVect)siv1, (ULongSparseIntVect)siv2, (float)a, (float)b [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tversky similarity between two vectors

    C++ signature :double TverskySimilarity(RDKit::SparseIntVect<unsigned long>,RDKit::SparseIntVect<unsigned long>,double,double [,bool=False [,double=0.0]])
    """
    ...

@overload
def TverskySimilarity(
    self, bv1: ExplicitBitVect, pkl: str, a: float, b: float, returnDistance: bool = 0
) -> float:
    """
    return the Tversky similarity between two vectors

    C++ signature :double TverskySimilarity(RDKit::SparseIntVect<int>,RDKit::SparseIntVect<int>,double,double [,bool=False [,double=0.0]])

    TverskySimilarity( (LongSparseIntVect)siv1, (LongSparseIntVect)siv2, (float)a, (float)b [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tversky similarity between two vectors

    C++ signature :double TverskySimilarity(RDKit::SparseIntVect<long>,RDKit::SparseIntVect<long>,double,double [,bool=False [,double=0.0]])

    TverskySimilarity( (UIntSparseIntVect)siv1, (UIntSparseIntVect)siv2, (float)a, (float)b [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tversky similarity between two vectors

    C++ signature :double TverskySimilarity(RDKit::SparseIntVect<unsigned int>,RDKit::SparseIntVect<unsigned int>,double,double [,bool=False [,double=0.0]])

    TverskySimilarity( (ULongSparseIntVect)siv1, (ULongSparseIntVect)siv2, (float)a, (float)b [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tversky similarity between two vectors

    C++ signature :double TverskySimilarity(RDKit::SparseIntVect<unsigned long>,RDKit::SparseIntVect<unsigned long>,double,double [,bool=False [,double=0.0]])
    """
    ...

@overload
def TverskySimilarity(
    self,
    siv1: IntSparseIntVect,
    siv2: IntSparseIntVect,
    a: float,
    b: float,
    returnDistance: bool = False,
    bounds: float = 0.0,
) -> float:
    """
    return the Tversky similarity between two vectors

    C++ signature :double TverskySimilarity(RDKit::SparseIntVect<long>,RDKit::SparseIntVect<long>,double,double [,bool=False [,double=0.0]])

    TverskySimilarity( (UIntSparseIntVect)siv1, (UIntSparseIntVect)siv2, (float)a, (float)b [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tversky similarity between two vectors

    C++ signature :double TverskySimilarity(RDKit::SparseIntVect<unsigned int>,RDKit::SparseIntVect<unsigned int>,double,double [,bool=False [,double=0.0]])

    TverskySimilarity( (ULongSparseIntVect)siv1, (ULongSparseIntVect)siv2, (float)a, (float)b [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tversky similarity between two vectors

    C++ signature :double TverskySimilarity(RDKit::SparseIntVect<unsigned long>,RDKit::SparseIntVect<unsigned long>,double,double [,bool=False [,double=0.0]])
    """
    ...

@overload
def TverskySimilarity(
    self,
    siv1: LongSparseIntVect,
    siv2: LongSparseIntVect,
    a: float,
    b: float,
    returnDistance: bool = False,
    bounds: float = 0.0,
) -> float:
    """
    return the Tversky similarity between two vectors

    C++ signature :double TverskySimilarity(RDKit::SparseIntVect<unsigned int>,RDKit::SparseIntVect<unsigned int>,double,double [,bool=False [,double=0.0]])

    TverskySimilarity( (ULongSparseIntVect)siv1, (ULongSparseIntVect)siv2, (float)a, (float)b [, (bool)returnDistance=False [, (float)bounds=0.0]]) -> float :return the Tversky similarity between two vectors

    C++ signature :double TverskySimilarity(RDKit::SparseIntVect<unsigned long>,RDKit::SparseIntVect<unsigned long>,double,double [,bool=False [,double=0.0]])
    """
    ...

@overload
def TverskySimilarity(
    self,
    siv1: UIntSparseIntVect,
    siv2: UIntSparseIntVect,
    a: float,
    b: float,
    returnDistance: bool = False,
    bounds: float = 0.0,
) -> float:
    """
    return the Tversky similarity between two vectors

    C++ signature :double TverskySimilarity(RDKit::SparseIntVect<unsigned long>,RDKit::SparseIntVect<unsigned long>,double,double [,bool=False [,double=0.0]])
    """
    ...

@overload
def TverskySimilarity(
    self,
    siv1: ULongSparseIntVect,
    siv2: ULongSparseIntVect,
    a: float,
    b: float,
    returnDistance: bool = False,
    bounds: float = 0.0,
) -> float: ...
