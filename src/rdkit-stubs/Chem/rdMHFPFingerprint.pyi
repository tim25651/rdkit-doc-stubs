"""
rdkit.Chem.rdMHFPFingerprint module¶
"""
from typing import Any, ClassVar

import Boost.Python
from rdkit.Chem.rdchem import Mol
from rdkit.DataStructs.cDataStructs import ExplicitBitVect
from rdkit.rdBase import (
    _vectj,
    _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE,
    _vectSt6vectorIjSaIjEE,
)

class MHFPEncoder(Boost.Python.instance):
    """
    C++ signature :void __init__(_object* [,unsigned int [,unsigned int]])"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def CreateShinglingFromMol(
        self,
        arg1: MHFPEncoder,
        mol: Mol,
        radius: int = 3,
        rings: bool = True,
        isomeric: bool = False,
        kekulize: bool = False,
        min_radius: int = 1,
    ) -> _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE:
        """
        Creates a shingling (a list of circular n-grams / substructures) from a RDKit Mol instance.

        C++ signature :std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > CreateShinglingFromMol(RDKit::MHFPFingerprints::MHFPEncoder*,RDKit::ROMol [,unsigned char=3 [,bool=True [,bool=False [,bool=False [,unsigned char=1]]]]])
        """
        ...
    def CreateShinglingFromSmiles(
        self,
        arg1: MHFPEncoder,
        smiles: str,
        radius: int = 3,
        rings: bool = True,
        isomeric: bool = False,
        kekulize: bool = False,
        min_radius: int = 1,
    ) -> _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE:
        """
        Creates a shingling (a list of circular n-grams / substructures) from a SMILES string.

        C++ signature :std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > CreateShinglingFromSmiles(RDKit::MHFPFingerprints::MHFPEncoder*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,unsigned char=3 [,bool=True [,bool=False [,bool=False [,unsigned char=1]]]]])
        """
        ...
    def Distance(self, arg1: _vectj, arg2: _vectj) -> float:
        """
        C++ signature :double Distance(std::vector<unsigned int, std::allocator<unsigned int> >,std::vector<unsigned int, std::allocator<unsigned int> >)
        """
        ...
    def EncodeMol(
        self,
        arg1: MHFPEncoder,
        mol: Mol,
        radius: int = 3,
        rings: bool = True,
        isomeric: bool = False,
        kekulize: bool = False,
        min_radius: int = 1,
    ) -> _vectj:
        """
        Creates a MHFP vector from an RDKit Mol instance.

        C++ signature :std::vector<unsigned int, std::allocator<unsigned int> > EncodeMol(RDKit::MHFPFingerprints::MHFPEncoder*,RDKit::ROMol [,unsigned char=3 [,bool=True [,bool=False [,bool=False [,unsigned char=1]]]]])
        """
        ...
    def EncodeMolsBulk(
        self,
        arg1: MHFPEncoder,
        mols: list,
        radius: int = 3,
        rings: bool = True,
        isomeric: bool = False,
        kekulize: bool = False,
        min_radius: int = 1,
    ) -> _vectSt6vectorIjSaIjEE:
        """
        Creates a MHFP vector from a list of RDKit Mol instances.

        C++ signature :std::vector<std::vector<unsigned int, std::allocator<unsigned int> >, std::allocator<std::vector<unsigned int, std::allocator<unsigned int> > > > EncodeMolsBulk(RDKit::MHFPFingerprints::MHFPEncoder*,boost::python::list {lvalue} [,unsigned char=3 [,bool=True [,bool=False [,bool=False [,unsigned char=1]]]]])
        """
        ...
    def EncodeSECFPMol(
        self,
        arg1: MHFPEncoder,
        smiles: Mol,
        radius: int = 3,
        rings: bool = True,
        isomeric: bool = False,
        kekulize: bool = False,
        min_radius: int = 1,
        length: int = 2048,
    ) -> ExplicitBitVect:
        """
        Creates a SECFP binary vector from an RDKit Mol instance.

        C++ signature :ExplicitBitVect EncodeSECFPMol(RDKit::MHFPFingerprints::MHFPEncoder*,RDKit::ROMol [,unsigned char=3 [,bool=True [,bool=False [,bool=False [,unsigned char=1 [,unsigned long=2048]]]]]])
        """
        ...
    def EncodeSECFPMolsBulk(
        self,
        arg1: MHFPEncoder,
        smiles: list,
        radius: int = 3,
        rings: bool = True,
        isomeric: bool = False,
        kekulize: bool = False,
        min_radius: int = 1,
        length: int = 2048,
    ) -> object:
        """
        Creates a SECFP binary vector from a list of RDKit Mol instances.

        C++ signature :std::vector<ExplicitBitVect, std::allocator<ExplicitBitVect> > EncodeSECFPMolsBulk(RDKit::MHFPFingerprints::MHFPEncoder*,boost::python::list {lvalue} [,unsigned char=3 [,bool=True [,bool=False [,bool=False [,unsigned char=1 [,unsigned long=2048]]]]]])
        """
        ...
    def EncodeSECFPSmiles(
        self,
        arg1: MHFPEncoder,
        smiles: str,
        radius: int = 3,
        rings: bool = True,
        isomeric: bool = False,
        kekulize: bool = False,
        min_radius: int = 1,
        length: int = 2048,
    ) -> ExplicitBitVect:
        """
        Creates a SECFP binary vector from a SMILES string.

        C++ signature :ExplicitBitVect EncodeSECFPSmiles(RDKit::MHFPFingerprints::MHFPEncoder*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,unsigned char=3 [,bool=True [,bool=False [,bool=False [,unsigned char=1 [,unsigned long=2048]]]]]])
        """
        ...
    def EncodeSECFPSmilesBulk(
        self,
        arg1: MHFPEncoder,
        smiles: list,
        radius: int = 3,
        rings: bool = True,
        isomeric: bool = False,
        kekulize: bool = False,
        min_radius: int = 1,
        length: int = 2048,
    ) -> object:
        """
        Creates a SECFP binary vector from a list of SMILES strings.

        C++ signature :std::vector<ExplicitBitVect, std::allocator<ExplicitBitVect> > EncodeSECFPSmilesBulk(RDKit::MHFPFingerprints::MHFPEncoder*,boost::python::list {lvalue} [,unsigned char=3 [,bool=True [,bool=False [,bool=False [,unsigned char=1 [,unsigned long=2048]]]]]])
        """
        ...
    def EncodeSmiles(
        self,
        arg1: MHFPEncoder,
        smiles: str,
        radius: int = 3,
        rings: bool = True,
        isomeric: bool = False,
        kekulize: bool = False,
        min_radius: int = 1,
    ) -> _vectj:
        """
        Creates a MHFP vector from a SMILES string.

        C++ signature :std::vector<unsigned int, std::allocator<unsigned int> > EncodeSmiles(RDKit::MHFPFingerprints::MHFPEncoder*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,unsigned char=3 [,bool=True [,bool=False [,bool=False [,unsigned char=1]]]]])
        """
        ...
    def EncodeSmilesBulk(
        self,
        arg1: MHFPEncoder,
        smiles: list,
        radius: int = 3,
        rings: bool = True,
        isomeric: bool = False,
        kekulize: bool = False,
        min_radius: int = 1,
    ) -> _vectSt6vectorIjSaIjEE:
        """
        Creates a MHFP vector from a list of SMILES strings.

        C++ signature :std::vector<std::vector<unsigned int, std::allocator<unsigned int> >, std::allocator<std::vector<unsigned int, std::allocator<unsigned int> > > > EncodeSmilesBulk(RDKit::MHFPFingerprints::MHFPEncoder*,boost::python::list {lvalue} [,unsigned char=3 [,bool=True [,bool=False [,bool=False [,unsigned char=1]]]]])
        """
        ...
    def FromArray(self, vec: MHFPEncoder, mhfp: list[int]) -> _vectj:
        """
        C++ signature :std::vector<unsigned int, std::allocator<unsigned int> > FromArray(RDKit::MHFPFingerprints::MHFPEncoder*,boost::python::list {lvalue})
        """
        ...
    def FromStringArray(self, vec: MHFPEncoder, mhfp: list[str]) -> _vectj:
        """
        C++ signature :std::vector<unsigned int, std::allocator<unsigned int> > FromStringArray(RDKit::MHFPFingerprints::MHFPEncoder*,boost::python::list {lvalue})
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...
