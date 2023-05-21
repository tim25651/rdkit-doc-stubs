"""
rdkit.Chem.rdmolfiles module¶
Module containing RDKit functionality for working with molecular file formats.
"""
from typing import Any, ClassVar, overload

import Boost.Python
from rdkit import rdBase
from rdkit.Chem.rdchem import Atom, Bond, Mol
from rdkit.Chem.rdMolDescriptors import AtomPairsParameters
from rdkit.rdBase import (
    _vecti,
    _vectj,
    _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE,
)

class CXSmilesFields(Boost.Python.enum):
    CX_ALL: CXSmilesFields = ...
    CX_ATOM_LABELS: CXSmilesFields = ...
    CX_ATOM_PROPS: CXSmilesFields = ...
    CX_COORDS: CXSmilesFields = ...
    CX_ENHANCEDSTEREO: CXSmilesFields = ...
    CX_LINKNODES: CXSmilesFields = ...
    CX_MOLFILE_VALUES: CXSmilesFields = ...
    CX_NONE: CXSmilesFields = ...
    CX_POLYMER: CXSmilesFields = ...
    CX_RADICALS: CXSmilesFields = ...
    CX_SGROUPS: CXSmilesFields = ...
    names: dict[str, CXSmilesFields] = ...
    values: dict[int, CXSmilesFields] = ...
    __slots__: ClassVar[tuple] = ...

class ForwardSDMolSupplier(Boost.Python.instance):
    """
    A class which supplies molecules from file-like object containing SD data.
    Usage examples:

    Lazy evaluation: the molecules are not constructed until we ask for them:
    >>> suppl = ForwardSDMolSupplier(file('in.sdf'))
    >>> for mol in suppl:
    ...    if mol is not None: mol.GetNumAtoms()

    we can also read from compressed files:
    >>> import gzip
    >>> suppl = ForwardSDMolSupplier(gzip.open('in.sdf.gz'))
    >>> for mol in suppl:
    ...   if mol is not None: print mol.GetNumAtoms()

    Properties in the SD file are used to set properties on each molecule.
    The properties are accessible using the mol.GetProp(propName) method.

    C++ signature :void __init__(_object*,boost::python::api::object {lvalue} [,bool=True [,bool=True [,bool=True]]])

    __init__( (object)arg1, (streambuf)streambuf [, (bool)sanitize=True [, (bool)removeHs=True [, (bool)strictParsing=True]]]) -> None :

    C++ signature :void __init__(_object*,boost_adaptbx::python::streambuf {lvalue} [,bool=True [,bool=True [,bool=True]]])

    __init__( (object)arg1, (str)filename [, (bool)sanitize=True [, (bool)removeHs=True [, (bool)strictParsing=True]]]) -> None :

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=True [,bool=True [,bool=True]]])
    """

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def GetEOFHitOnRead(self, arg1: ForwardSDMolSupplier) -> bool:
        """
        Returns whether or EOF was hit while parsing the previous entry.

        C++ signature :bool GetEOFHitOnRead((anonymous namespace)::LocalForwardSDMolSupplier {lvalue})
        """
        ...
    def GetProcessPropertyLists(self, arg1: ForwardSDMolSupplier) -> bool:
        """
        returns whether or not any property lists that are present will be processed when reading molecules

        C++ signature :bool GetProcessPropertyLists((anonymous namespace)::LocalForwardSDMolSupplier {lvalue})
        """
        ...
    def SetProcessPropertyLists(self, arg1: ForwardSDMolSupplier, arg2: bool) -> None:
        """
        sets whether or not any property lists that are present will be processed when reading molecules

        C++ signature :void SetProcessPropertyLists((anonymous namespace)::LocalForwardSDMolSupplier {lvalue},bool)
        """
        ...
    def atEnd(self, arg1: ForwardSDMolSupplier) -> bool:
        """
        Returns whether or not we have hit EOF.

        C++ signature :bool atEnd((anonymous namespace)::LocalForwardSDMolSupplier {lvalue})
        """
        ...
    @classmethod
    def __enter__(cls, anonymousnamespace) -> Any: ...
    @classmethod
    def __exit__(cls, type, value, traceback) -> Any: ...
    @classmethod
    def __iter__(cls, anonymousnamespace) -> Any: ...
    @classmethod
    def __next__(cls, anonymousnamespace) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class MaeMolSupplier(Boost.Python.instance):
    """
    A class which supplies molecules from file-like object containing Maestro data.
    Usage examples:

    Lazy evaluation: the molecules are not constructed until we ask for them:
    >>> suppl = MaeMolSupplier(file('in.mae'))
    >>> for mol in suppl:
    ...    if mol is not None: mol.GetNumAtoms()

    we can also read from compressed files:
    >>> import gzip
    >>> suppl = MaeMolSupplier(gzip.open('in.maegz'))
    >>> for mol in suppl:
    ...   if mol is not None: print mol.GetNumAtoms()

    Properties in the Maestro file are used to set properties on each molecule.
    The properties are accessible using the mol.GetProp(propName) method.

    C++ signature :void __init__(_object*)

    __init__( (object)arg1, (AtomPairsParameters)fileobj [, (bool)sanitize=True [, (bool)removeHs=True]]) -> None :

    C++ signature :void __init__(_object*,boost::python::api::object {lvalue} [,bool=True [,bool=True]])

    __init__( (object)arg1, (streambuf)streambuf [, (bool)sanitize=True [, (bool)removeHs=True]]) -> None :

    C++ signature :void __init__(_object*,boost_adaptbx::python::streambuf {lvalue} [,bool=True [,bool=True]])

    __init__( (object)arg1, (str)filename [, (bool)sanitize=True [, (bool)removeHs=True]]) -> None :

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=True [,bool=True]])
    """

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def SetData(
        self,
        arg1: MaeMolSupplier,
        data: str,
        sanitize: bool = True,
        removeHs: bool = True,
    ) -> None:
        """
        Sets the text to be parsed

        C++ signature :void SetData((anonymous namespace)::LocalMaeMolSupplier {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=True [,bool=True]])
        """
        ...
    def atEnd(self, arg1: MaeMolSupplier) -> bool:
        """
        Returns whether or not we have hit EOF.

        C++ signature :bool atEnd((anonymous namespace)::LocalMaeMolSupplier {lvalue})
        """
        ...
    def reset(self, arg1: MaeMolSupplier) -> None:
        """
        Resets our position in the file to the beginning.

        C++ signature :void reset((anonymous namespace)::LocalMaeMolSupplier {lvalue})
        """
        ...
    @classmethod
    def __enter__(cls, anonymousnamespace) -> Any: ...
    @classmethod
    def __exit__(cls, type, value, traceback) -> Any: ...
    @classmethod
    def __getitem__(cls, anonymousnamespace, int) -> Any: ...
    @classmethod
    def __iter__(cls, anonymousnamespace) -> Any: ...
    @classmethod
    def __len__(cls, anonymousnamespace) -> Any: ...
    @classmethod
    def __next__(cls, anonymousnamespace) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class MaeWriter(Boost.Python.instance):
    """
    An experimental class for writing molecules to Maestro files.
    Usage examples:

    writing to a named file:
    >>> writer = MaeWriter('out.mae')
    >>> for mol in list_of_mols:
    ...    writer.write(mol)

    writing to a file-like object:
    >>> import gzip
    >>> outf=gzip.open('out.mae.gz','wt+')
    >>> writer = MaeWriter(outf)
    >>> for mol in list_of_mols:
    ...   writer.write(mol)
    >>> writer.close()
    >>> outf.close()

    By default all non-private molecule, atom and bond properties are written
    to the Maestro file. This can be changed using the SetProps method:
    >>> writer = MaeWriter('out.mae')
    >>> writer.SetProps(['prop1','prop2'])

    Properties that are specified, but are not present will be ignored.
    Kekulization is mandatory, as the Maestro format does not have
    the concept of an aromatic bond
    As this is an experimental writer, many features are not supported yet,
    e.g. chirality and bond stereo labels, stereo groups, substance groups,
    isotopes, or even dummy atoms. Note that these aren’t supported by
    MaeMolSupplier either.

    C++ signature :void __init__(_object*,boost::python::api::object {lvalue})

    __init__( (object)arg1, (streambuf)streambuf) -> None :

    C++ signature :void __init__(_object*,boost_adaptbx::python::streambuf {lvalue})

    __init__( (object)arg1, (str)filename) -> None :

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def GetText(
        self,
        mol: Mol,
        heavyAtomColor: str = "A0A0A0",
        confId: int = -1,
        props_list: rdBase._vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE = ...,
    ) -> str:
        """
        returns the Maestro ct block text for a molecule

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetText(RDKit::ROMol [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’A0A0A0’ [,int=-1 [,std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >=<rdkit.rdBase._vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE object at 0x7fd481d023c0>]]])
        """
        ...
    def NumMols(self, arg1: MaeWriter) -> int:
        """
        Returns the number of molecules written so far.

        C++ signature :unsigned int NumMols(RDKit::LocalMaeWriter {lvalue})"""
        ...
    def SetProps(
        self: MaeWriter,
        props_list: _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE,
    ) -> None:
        """
        Sets the atom and mol properties to be written to the output file

        ARGUMENTS:

        props: a list or tuple of atom and mol property names

        C++ signature :void SetProps(RDKit::LocalMaeWriter {lvalue},std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >)
        """
        ...
    def close(self, arg1: MaeWriter) -> None:
        """
        Flushes the output file and closes it. The Writer cannot be used after this.

        C++ signature :void close(RDKit::LocalMaeWriter {lvalue})"""
        ...
    def flush(self, arg1: MaeWriter) -> None:
        """
        Flushes the output file (forces the disk file to be updated).

        C++ signature :void flush(RDKit::LocalMaeWriter {lvalue})"""
        ...
    def write(
        self: MaeWriter, mol: Mol, heavyAtomColor: str = "A0A0A0", confId: int = -1
    ) -> None:
        """
        Writes a molecule to the output file.

        ARGUMENTS:

        mol: the Mol to be written
        heavyAtomColor: (optional) color which heavy atoms will have in Maestro
        confId: (optional) ID of the conformation to write

        C++ signature :void write(RDKit::LocalMaeWriter {lvalue},RDKit::ROMol [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’A0A0A0’ [,int=-1]])
        """
        ...
    @classmethod
    def __enter__(cls, RDKit) -> Any: ...
    @classmethod
    def __exit__(cls, type, value, traceback) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class MultithreadedSDMolSupplier(Boost.Python.instance):
    """
    A class which concurrently supplies molecules from a text file.
    Please note that this class is still a bit experimental and the API may
    change in future releases.
    Usage examples:

    Lazy evaluation: the molecules might not be constructed until we ask for them:
    >>> suppl = MultithreadedSDMolSupplier('in.sdf')
    >>> for mol in suppl:
    ...    if(mol):
    ...      mol.GetNumAtoms()

    Lazy evaluation 2:
    >>> suppl = MultithreadedSDMolSupplier('in.sdf')
    >>> while (!suppl.atEnd()):
    ...    mol = next(mol)
    ...    if(mol):
    ...      mol.GetNumAtoms()

    C++ signature :void __init__(_object*)

    __init__( (object)arg1, (str)fileName [, (bool)sanitize=True [, (bool)removeHs=True [, (bool)strictParsing=True [, (int)numWriterThreads=1 [, (int)sizeInputQueue=5 [, (int)sizeOutputQueue=5]]]]]]) -> None :Constructor

    ARGUMENTS:

    fileName: name of the file to be read
    sanitize: (optional) toggles sanitization of molecules as they are read.
    Defaults to true.
    removeHs: (optional) removes Hs. Defaults to true.
    strictParsing: (optional) allows strict or lax parsing. Defaults to true.
    numWriterThreads: (optional) number of writer threads. Defaults to 1.
    sizeInputQueue: (optional) size of input/reader queue. Defaults to 5.
    sizeOutputQueue: (optional) size of output/writer queue. Defaults to 5.

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=True [,bool=True [,bool=True [,unsigned int=1 [,unsigned long=5 [,unsigned long=5]]]]]])
    """

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def GetLastItemText(self, arg1: MultithreadedSDMolSupplier) -> str:
        """
        Returns the text for the last extracted item.

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetLastItemText(RDKit::MultithreadedSDMolSupplier*)
        """
        ...
    def GetLastRecordId(self, arg1: MultithreadedSDMolSupplier) -> int:
        """
        Returns the record id for the last extracted item.

        C++ signature :unsigned int GetLastRecordId(RDKit::MultithreadedSDMolSupplier*)
        """
        ...
    def GetProcessPropertyLists(self, arg1: MultithreadedSDMolSupplier) -> bool:
        """
        returns whether or not any property lists that are present will be processed when reading molecules

        C++ signature :bool GetProcessPropertyLists(RDKit::MultithreadedSDMolSupplier {lvalue})
        """
        ...
    def SetProcessPropertyLists(
        self, arg1: MultithreadedSDMolSupplier, arg2: bool
    ) -> None:
        """
        sets whether or not any property lists that are present will be processed when reading molecules

        C++ signature :void SetProcessPropertyLists(RDKit::MultithreadedSDMolSupplier {lvalue},bool)
        """
        ...
    def atEnd(self, arg1: MultithreadedSDMolSupplier) -> bool:
        """
        Returns true if we have read all records else false.

        C++ signature :bool atEnd(RDKit::MultithreadedSDMolSupplier {lvalue})"""
        ...
    @classmethod
    def __enter__(cls, RDKit) -> Any: ...
    @classmethod
    def __exit__(cls, type, value, traceback) -> Any: ...
    @classmethod
    def __iter__(cls, RDKit) -> Any: ...
    @classmethod
    def __next__(cls, RDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class MultithreadedSmilesMolSupplier(Boost.Python.instance):
    """
    A class which concurrently supplies molecules from a text file.
    Please note that this class is still a bit experimental and the API may
    change in future releases.
    Usage examples:

    Lazy evaluation: the molecules might not be constructed until we ask for them:
    >>> suppl = MultithreadedSmilesMolSupplier('in.smi')
    >>> for mol in suppl:
    ...    if(mol):
    ...      mol.GetNumAtoms()

    Lazy evaluation 2:
    >>> suppl = MultithreadedSmilesMolSupplier('in.smi')
    >>> while (!suppl.atEnd()):
    ...    mol = next(mol)
    ...    if(mol):
    ...      mol.GetNumAtoms()

    C++ signature :void __init__(_object*)

    __init__( (object)arg1, (str)fileName [, (str)delimiter=’ t’ [, (int)smilesColumn=0 [, (int)nameColumn=1 [, (bool)titleLine=True [, (bool)sanitize=True [, (int)numWriterThreads=1 [, (int)sizeInputQueue=5 [, (int)sizeOutputQueue=5]]]]]]]]) -> None :Constructor

    ARGUMENTS:

    fileName: name of the file to be read
    delimiter: (optional) text delimiter (a string).  Defauts to ‘        ‘.
    smilesColumn: (optional) index of the column containing the SMILES
    data.  Defaults to 0.
    nameColumn: (optional) index of the column containing molecule names.
    Defaults to 1.
    titleLine: (optional) set this toggle if the file contains a title line.
    Defaults to true.
    sanitize: (optional) toggles sanitization of molecules as they are read.
    Defaults to true.
    numWriterThreads: (optional) number of writer threads. Defaults to 1.
    sizeInputQueue: (optional) size of input/reader queue. Defaults to 5.
    sizeOutputQueue: (optional) size of output/writer queue. Defaults to 5.

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’ t’ [,int=0 [,int=1 [,bool=True [,bool=True [,unsigned int=1 [,unsigned long=5 [,unsigned long=5]]]]]]]])
    """

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def GetLastItemText(self, arg1: MultithreadedSmilesMolSupplier) -> str:
        """
        Returns the text for the last extracted item.

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetLastItemText(RDKit::MultithreadedSmilesMolSupplier*)
        """
        ...
    def GetLastRecordId(self, arg1: MultithreadedSmilesMolSupplier) -> int:
        """
        Returns the record id for the last extracted item.

        C++ signature :unsigned int GetLastRecordId(RDKit::MultithreadedSmilesMolSupplier*)
        """
        ...
    def atEnd(self, arg1: MultithreadedSmilesMolSupplier) -> bool:
        """
        Returns true if we have read all records else false.

        C++ signature :bool atEnd(RDKit::MultithreadedSmilesMolSupplier {lvalue})"""
        ...
    @classmethod
    def __enter__(cls, RDKit) -> Any: ...
    @classmethod
    def __exit__(cls, type, value, traceback) -> Any: ...
    @classmethod
    def __iter__(cls, RDKit) -> Any: ...
    @classmethod
    def __next__(cls, RDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class PDBWriter(Boost.Python.instance):
    """
    A class for writing molecules to PDB files.

    C++ signature :void* __init__(boost::python::api::object,boost::python::api::object {lvalue} [,unsigned int=0])

    __init__( (object)arg1, (str)fileName [, (int)flavor=0]) -> None :Constructor.

    ARGUMENTS:

    fileName: name of the output file. (‘-’ to write to stdout)
    flavor: (optional)

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,unsigned int=0])
    """

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def NumMols(self, arg1: PDBWriter) -> int:
        """
        Returns the number of molecules written so far.

        C++ signature :unsigned int NumMols(RDKit::PDBWriter {lvalue})"""
        ...
    def close(self, arg1: PDBWriter) -> None:
        """
        Flushes the output file and closes it. The Writer cannot be used after this.

        C++ signature :void close(RDKit::PDBWriter {lvalue})"""
        ...
    def flush(self, arg1: PDBWriter) -> None:
        """
        Flushes the output file (forces the disk file to be updated).

        C++ signature :void flush(RDKit::PDBWriter {lvalue})"""
        ...
    def write(self: PDBWriter, mol: Mol, confId: int = -1) -> None:
        """
        Writes a molecule to the output file.

        ARGUMENTS:

        mol: the Mol to be written
        confId: (optional) ignored

        C++ signature :void write(RDKit::PDBWriter {lvalue},RDKit::ROMol [,int=-1])"""
        ...
    @classmethod
    def __enter__(cls, RDKit) -> Any: ...
    @classmethod
    def __exit__(cls, type, value, traceback) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class SDMolSupplier(Boost.Python.instance):
    """
    A class which supplies molecules from an SD file.
    Usage examples:

    Lazy evaluation: the molecules are not constructed until we ask for them:
    >>> suppl = SDMolSupplier('in.sdf')
    >>> for mol in suppl:
    ...    mol.GetNumAtoms()

    Lazy evaluation 2:
    >>> suppl = SDMolSupplier('in.sdf')
    >>> mol1 = next(suppl)
    >>> mol2 = next(suppl)
    >>> suppl.reset()
    >>> mol3 = next(suppl)
    # mol3 and mol1 are the same:
    >>> MolToSmiles(mol3)==MolToSmiles(mol1)

    Random Access:
    >>> suppl = SDMolSupplier('in.sdf')
    >>> mol1 = suppl[0]
    >>> mol2 = suppl[1]
    # NOTE: this will generate an IndexError if the supplier doesn't have that many
    molecules.

    Random Access 2:  looping over all molecules
    >>> suppl = SDMolSupplier('in.sdf')
    >>> nMols = len(suppl)
    >>> for i in range(nMols):
    ...   suppl[i].GetNumAtoms()

    Properties in the SD file are used to set properties on each molecule.
    The properties are accessible using the mol.GetProp(propName) method.

    C++ signature :void __init__(_object*)

    __init__( (object)arg1, (str)fileName [, (bool)sanitize=True [, (bool)removeHs=True [, (bool)strictParsing=True]]]) -> None :

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=True [,bool=True [,bool=True]]])
    """

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def GetItemText(self: SDMolSupplier, index: int) -> str:
        """
        returns the text for an item

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetItemText(RDKit::SDMolSupplier {lvalue},unsigned int)
        """
        ...
    def GetProcessPropertyLists(self, arg1: SDMolSupplier) -> bool:
        """
        returns whether or not any property lists that are present will be processed when reading molecules

        C++ signature :bool GetProcessPropertyLists(RDKit::SDMolSupplier {lvalue})"""
        ...
    def SetData(
        self: SDMolSupplier,
        data: str,
        sanitize: bool = True,
        removeHs: bool = True,
        strictParsing: bool = True,
    ) -> None:
        """
        Sets the text to be parsed

        C++ signature :void SetData(RDKit::SDMolSupplier {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=True [,bool=True [,bool=True]]])
        """
        ...
    def SetProcessPropertyLists(self, arg1: SDMolSupplier, arg2: bool) -> None:
        """
        sets whether or not any property lists that are present will be processed when reading molecules

        C++ signature :void SetProcessPropertyLists(RDKit::SDMolSupplier {lvalue},bool)
        """
        ...
    @classmethod
    def _SetStreamIndices(cls, RDKit, boost) -> Any: ...
    def atEnd(self, arg1: SDMolSupplier) -> bool:
        """
        Returns whether or not we have hit EOF.

        C++ signature :bool atEnd(RDKit::SDMolSupplier {lvalue})"""
        ...
    def reset(self, arg1: SDMolSupplier) -> None:
        """
        Resets our position in the file to the beginning.

        C++ signature :void reset(RDKit::SDMolSupplier {lvalue})"""
        ...
    @classmethod
    def __enter__(cls, RDKit) -> Any: ...
    @classmethod
    def __exit__(cls, type, value, traceback) -> Any: ...
    @classmethod
    def __getitem__(cls, RDKit, int) -> Any: ...
    @classmethod
    def __iter__(cls, RDKit) -> Any: ...
    @classmethod
    def __len__(cls, RDKit) -> Any: ...
    @classmethod
    def __next__(cls, RDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class SDWriter(Boost.Python.instance):
    """
    A class for writing molecules to SD files.
    Usage examples:

    writing to a named file:
    >>> writer = SDWriter('out.sdf')
    >>> for mol in list_of_mols:
    ...    writer.write(mol)

    writing to a file-like object:
    >>> import gzip
    >>> outf=gzip.open('out.sdf.gz','wt+')
    >>> writer = SDWriter(outf)
    >>> for mol in list_of_mols:
    ...   writer.write(mol)
    >>> writer.close()
    >>> outf.close()

    By default all non-private molecular properties are written to the SD file.
    This can be changed using the SetProps method:

    >>> writer = SDWriter('out.sdf')
    >>> writer.SetProps(['prop1','prop2'])

    C++ signature :void* __init__(boost::python::api::object,boost::python::api::object {lvalue})

    __init__( (object)arg1, (str)fileName) -> None :Constructor.

    If a string argument is provided, it will be treated as the name of the output file.
    If a file-like object is provided, output will be sent there.

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def GetForceV3000(self, arg1: SDWriter) -> bool:
        """
        Returns whether or not V3000 mol file writing is being forced.

        C++ signature :bool GetForceV3000(RDKit::SDWriter {lvalue})"""
        ...
    def GetKekulize(self, arg1: SDWriter) -> bool:
        """
        Returns whether or not molecules are kekulized on writing.

        C++ signature :bool GetKekulize(RDKit::SDWriter {lvalue})"""
        ...
    def GetText(
        self,
        mol: Mol,
        confId: int = -1,
        kekulize: bool = True,
        force_v3000: bool = False,
        molid: int = -1,
    ) -> str:
        """
        returns the SD text for a molecule

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetText(RDKit::ROMol [,int=-1 [,bool=True [,bool=False [,int=-1]]]])
        """
        ...
    def NumMols(self, arg1: SDWriter) -> int:
        """
        Returns the number of molecules written so far.

        C++ signature :unsigned int NumMols(RDKit::SDWriter {lvalue})"""
        ...
    def SetForceV3000(self, arg1: SDWriter, arg2: bool) -> None:
        """
        Sets whether or not V3000 mol file writing is being forced.

        C++ signature :void SetForceV3000(RDKit::SDWriter {lvalue},bool)"""
        ...
    def SetKekulize(self, arg1: SDWriter, arg2: bool) -> None:
        """
        Sets whether or not molecules are kekulized on writing.

        C++ signature :void SetKekulize(RDKit::SDWriter {lvalue},bool)"""
        ...
    def SetProps(self, arg1: SDWriter, arg2: AtomPairsParameters) -> None:
        """
        Sets the properties to be written to the output file

        ARGUMENTS:

        props: a list or tuple of property names

        C++ signature :void SetProps(RDKit::SDWriter {lvalue},boost::python::api::object)
        """
        ...
    def close(self, arg1: SDWriter) -> None:
        """
        Flushes the output file and closes it. The Writer cannot be used after this.

        C++ signature :void close(RDKit::SDWriter {lvalue})"""
        ...
    def flush(self, arg1: SDWriter) -> None:
        """
        Flushes the output file (forces the disk file to be updated).

        C++ signature :void flush(RDKit::SDWriter {lvalue})"""
        ...
    def write(self: SDWriter, mol: Mol, confId: int = -1) -> None:
        """
        Writes a molecule to the output file.

        ARGUMENTS:

        mol: the Mol to be written
        confId: (optional) ID of the conformation to write

        C++ signature :void write(RDKit::SDWriter {lvalue},RDKit::ROMol {lvalue} [,int=-1])
        """
        ...
    @classmethod
    def __enter__(cls, RDKit) -> Any: ...
    @classmethod
    def __exit__(cls, type, value, traceback) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class SmartsParserParams(Boost.Python.instance):
    """
    Parameters controlling SMARTS Parsing

    C++ signature :void __init__(_object*)

    property allowCXSMILES¶
    controls whether or not the CXSMILES extensions are parsed

    property debugParse¶
    controls the amount of debugging information produced

    property mergeHs¶
    toggles merging H atoms in the SMARTS into neighboring atoms

    property parseName¶
    controls whether or not the molecule name is also parsed

    property strictCXSMILES¶
    controls whether or not problems in CXSMILES parsing causes molecule parsing to fail
    """

    ...
    __instance_size__: ClassVar[int] = ...
    allowCXSMILES: Any
    debugParse: Any
    mergeHs: Any
    parseName: Any
    strictCXSMILES: Any

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class SmilesMolSupplier(Boost.Python.instance):
    """
    A class which supplies molecules from a text file.
    Usage examples:

    Lazy evaluation: the molecules are not constructed until we ask for them:
    >>> suppl = SmilesMolSupplier('in.smi')
    >>> for mol in suppl:
    ...    mol.GetNumAtoms()

    Lazy evaluation 2:
    >>> suppl = SmilesMolSupplier('in.smi')
    >>> mol1 = next(suppl)
    >>> mol2 = next(suppl)
    >>> suppl.reset()
    >>> mol3 = next(suppl)
    # mol3 and mol1 are the same:
    >>> MolToSmiles(mol3)==MolToSmiles(mol1)

    Random Access:  all molecules are constructed as soon as we ask for the
    length:
    >>> suppl = SmilesMolSupplier('in.smi')
    >>> nMols = len(suppl)
    >>> for i in range(nMols):
    ...   suppl[i].GetNumAtoms()

    If the input file has a title line and more than two columns (smiles and id), the
    additional columns will be used to set properties on each molecule.  The properties
    are accessible using the mol.GetProp(propName) method.
    Constructor

    ARGUMENTS:

    fileName: name of the file to be read
    delimiter: (optional) text delimiter (a string).  Defauts to ‘ ‘.
    smilesColumn: (optional) index of the column containing the SMILES
    data.  Defaults to 0.
    nameColumn: (optional) index of the column containing molecule names.
    Defaults to 1.
    titleLine: (optional) set this toggle if the file contains a title line.
    Defaults to 1.
    sanitize: (optional) toggles sanitization of molecules as they are read.
    Defaults to 1.

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’ ‘ [,int=0 [,int=1 [,bool=True [,bool=True]]]]])

    __init__( (object)arg1) -> None :

    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def GetItemText(self: SmilesMolSupplier, index: int) -> str:
        """
        returns the text for an item

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetItemText(RDKit::SmilesMolSupplier {lvalue},unsigned int)
        """
        ...
    def SetData(
        self: SmilesMolSupplier,
        data: str,
        delimiter: str = " ",
        smilesColumn: int = 0,
        nameColumn: int = 1,
        titleLine: bool = True,
        sanitize: bool = True,
    ) -> None:
        """
        Sets the text to be parsed

        C++ signature :void SetData(RDKit::SmilesMolSupplier {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’ ‘ [,int=0 [,int=1 [,bool=True [,bool=True]]]]])
        """
        ...
    def reset(self, arg1: SmilesMolSupplier) -> None:
        """
        Resets our position in the file to the beginning.

        C++ signature :void reset(RDKit::SmilesMolSupplier {lvalue})"""
        ...
    @classmethod
    def __enter__(cls, RDKit) -> Any: ...
    @classmethod
    def __exit__(cls, type, value, traceback) -> Any: ...
    @classmethod
    def __getitem__(cls, RDKit, int) -> Any: ...
    @classmethod
    def __iter__(cls, RDKit) -> Any: ...
    @classmethod
    def __len__(cls, RDKit) -> Any: ...
    @classmethod
    def __next__(cls, RDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class SmilesParserParams(Boost.Python.instance):
    """
    Parameters controlling SMILES Parsing

    C++ signature :void __init__(_object*)

    property allowCXSMILES¶
    controls whether or not the CXSMILES extensions are parsed

    property debugParse¶
    controls the amount of debugging information produced

    property parseName¶
    controls whether or not the molecule name is also parsed

    property removeHs¶
    controls whether or not Hs are removed before the molecule is returned

    property sanitize¶
    controls whether or not the molecule is sanitized before being returned

    property strictCXSMILES¶
    controls whether or not problems in CXSMILES parsing causes molecule parsing to fail
    """

    ...
    __instance_size__: ClassVar[int] = ...
    allowCXSMILES: Any
    debugParse: Any
    parseName: Any
    removeHs: Any
    sanitize: Any
    strictCXSMILES: Any

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class SmilesWriteParams(Boost.Python.instance):
    """
    Parameters controlling SMILES writing

    C++ signature :void __init__(_object*)

    property allBondsExplicit¶
    include symbols for all bonds

    property allHsExplicit¶
    provide hydrogen counts for every atom

    property canonical¶
    generate canonical SMILES

    property doIsomericSmiles¶
    include stereochemistry and isotope information

    property doKekule¶
    kekulize the molecule before generating the SMILES and output single/double bonds. NOTE that the output is not canonical and that this will thrown an exception if the molecule cannot be kekulized

    property doRandom¶
    randomize the output order. The resulting SMILES is not canonical

    property rootedAtAtom¶
    make sure the SMILES starts at the specified atom. The resulting SMILES is not canonical
    """

    ...
    __instance_size__: ClassVar[int] = ...
    allBondsExplicit: Any
    allHsExplicit: Any
    canonical: Any
    doIsomericSmiles: Any
    doKekule: Any
    doRandom: Any
    rootedAtAtom: Any

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class SmilesWriter(Boost.Python.instance):
    """
    A class for writing molecules to text files.

    C++ signature :void* __init__(boost::python::api::object,boost::python::api::object {lvalue} [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’ ‘ [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’Name’ [,bool=True [,bool=True [,bool=False]]]]])

    __init__( (object)arg1, (str)fileName [, (str)delimiter=’ ‘ [, (str)nameHeader=’Name’ [, (bool)includeHeader=True [, (bool)isomericSmiles=True [, (bool)kekuleSmiles=False]]]]]) -> None :Constructor.

    ARGUMENTS:

    fileName: name of the output file. (‘-’ to write to stdout)
    delimiter: (optional) delimiter to be used to separate entries on each line.

    nameHeader: (optional) text to use for the name column in the header line.If this is blank, names will not be included in the output.

    includeHeader: (optional) toggles inclusion of a header line in the output file.
    isomericSmiles: (optional) toggles output of isomeric smiles (includes stereochem information).
    kekuleSmiles: (optional) toggles output of kekule smiles (no aromatic bonds for molecules that have been kekulized).

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’ ‘ [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’Name’ [,bool=True [,bool=True [,bool=False]]]]])
    """

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def NumMols(self, arg1: SmilesWriter) -> int:
        """
        Returns the number of molecules written so far.

        C++ signature :unsigned int NumMols(RDKit::SmilesWriter {lvalue})"""
        ...
    def SetProps(self, arg1: SmilesWriter, arg2: AtomPairsParameters) -> None:
        """
        Sets the properties to be written to the output file

        ARGUMENTS:

        props: a list or tuple of property names

        C++ signature :void SetProps(RDKit::SmilesWriter {lvalue},boost::python::api::object)
        """
        ...
    def close(self, arg1: SmilesWriter) -> None:
        """
        Flushes the output file and closes it. The Writer cannot be used after this.

        C++ signature :void close(RDKit::SmilesWriter {lvalue})"""
        ...
    def flush(self, arg1: SmilesWriter) -> None:
        """
        Flushes the output file (forces the disk file to be updated).

        C++ signature :void flush(RDKit::SmilesWriter {lvalue})"""
        ...
    def write(self: SmilesWriter, mol: Mol, confId: int = -1) -> None:
        """
        Writes a molecule to the output file.

        ARGUMENTS:

        mol: the Mol to be written
        confId: (optional) ignored

        C++ signature :void write(RDKit::SmilesWriter {lvalue},RDKit::ROMol [,int=-1])
        """
        ...
    @classmethod
    def __enter__(cls, RDKit) -> Any: ...
    @classmethod
    def __exit__(cls, type, value, traceback) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class TDTMolSupplier(Boost.Python.instance):
    """
    A class which supplies molecules from a TDT file.
    Usage examples:

    Lazy evaluation: the molecules are not constructed until we ask for them:
    >>> suppl = TDTMolSupplier('in.smi')
    >>> for mol in suppl:
    ...    mol.GetNumAtoms()

    Lazy evaluation 2:
    >>> suppl = TDTMolSupplier('in.smi')
    >>> mol1 = next(suppl)
    >>> mol2 = next(suppl)
    >>> suppl.reset()
    >>> mol3 = next(suppl)

    # mol3 and mol1 are the same:       >>> MolToSmiles(mol3)==MolToSmiles(mol1)

    Random Access:  all molecules are constructed as soon as we ask for the
    length:
    >>> suppl = TDTMolSupplier('in.smi')
    >>> nMols = len(suppl)
    >>> for i in range(nMols):
    ...   suppl[i].GetNumAtoms()

    Properties in the file are used to set properties on each molecule.
    The properties are accessible using the mol.GetProp(propName) method.

    C++ signature :void __init__(_object*)

    __init__( (object)arg1, (str)fileName [, (str)nameRecord=’’ [, (int)confId2D=-1 [, (int)confId3D=-1 [, (bool)sanitize=True]]]]) -> None :

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’’ [,int=-1 [,int=-1 [,bool=True]]]])
    """

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def GetItemText(self: TDTMolSupplier, index: int) -> str:
        """
        returns the text for an item

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetItemText(RDKit::TDTMolSupplier {lvalue},unsigned int)
        """
        ...
    def SetData(
        self: TDTMolSupplier,
        data: str,
        nameRecord: str = "",
        confId2D: int = -1,
        confId3D: int = -1,
        sanitize: bool = True,
    ) -> None:
        """
        Sets the text to be parsed

        C++ signature :void SetData(RDKit::TDTMolSupplier {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’’ [,int=-1 [,int=-1 [,bool=True]]]])
        """
        ...
    def reset(self, arg1: TDTMolSupplier) -> None:
        """
        Resets our position in the file to the beginning.

        C++ signature :void reset(RDKit::TDTMolSupplier {lvalue})"""
        ...
    @classmethod
    def __enter__(cls, RDKit) -> Any: ...
    @classmethod
    def __exit__(cls, type, value, traceback) -> Any: ...
    @classmethod
    def __getitem__(cls, RDKit, int) -> Any: ...
    @classmethod
    def __iter__(cls, RDKit) -> Any: ...
    @classmethod
    def __len__(cls, RDKit) -> Any: ...
    @classmethod
    def __next__(cls, RDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class TDTWriter(Boost.Python.instance):
    """
    A class for writing molecules to TDT files.

    C++ signature :void* __init__(boost::python::api::object,boost::python::api::object {lvalue})

    __init__( (object)arg1, (str)fileName) -> None :Constructor.

    If a string argument is provided, it will be treated as the name of the output file.
    If a file-like object is provided, output will be sent there.

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def GetNumDigits(self, arg1: TDTWriter) -> int:
        """
        C++ signature :unsigned int GetNumDigits(RDKit::TDTWriter {lvalue})"""
        ...
    def GetWrite2D(self, arg1: TDTWriter) -> bool:
        """
        C++ signature :bool GetWrite2D(RDKit::TDTWriter {lvalue})"""
        ...
    def GetWriteNames(self, arg1: TDTWriter) -> bool:
        """
        C++ signature :bool GetWriteNames(RDKit::TDTWriter {lvalue})"""
        ...
    def NumMols(self, arg1: TDTWriter) -> int:
        """
        Returns the number of molecules written so far.

        C++ signature :unsigned int NumMols(RDKit::TDTWriter {lvalue})"""
        ...
    def SetNumDigits(self, arg1: TDTWriter, arg2: int) -> None:
        """
        sets the number of digits to be written for coordinates

        C++ signature :void SetNumDigits(RDKit::TDTWriter {lvalue},unsigned int)"""
        ...
    def SetProps(self, arg1: TDTWriter, arg2: AtomPairsParameters) -> None:
        """
        Sets the properties to be written to the output file

        ARGUMENTS:

        props: a list or tuple of property names

        C++ signature :void SetProps(RDKit::TDTWriter {lvalue},boost::python::api::object)
        """
        ...
    def SetWrite2D(self: TDTWriter, state: bool = True) -> None:
        """
        causes 2D conformations to be written (default is 3D conformations)

        C++ signature :void SetWrite2D(RDKit::TDTWriter {lvalue} [,bool=True])"""
        ...
    def SetWriteNames(self: TDTWriter, state: bool = True) -> None:
        """
        causes names to be written to the output file as NAME records

        C++ signature :void SetWriteNames(RDKit::TDTWriter {lvalue} [,bool=True])"""
        ...
    def close(self, arg1: TDTWriter) -> None:
        """
        Flushes the output file and closes it. The Writer cannot be used after this.

        C++ signature :void close(RDKit::TDTWriter {lvalue})"""
        ...
    def flush(self, arg1: TDTWriter) -> None:
        """
        Flushes the output file (forces the disk file to be updated).

        C++ signature :void flush(RDKit::TDTWriter {lvalue})"""
        ...
    def write(self: TDTWriter, mol: Mol, confId: int = -1) -> None:
        """
        Writes a molecule to the output file.

        ARGUMENTS:

        mol: the Mol to be written
        confId: (optional) ID of the conformation to write

        C++ signature :void write(RDKit::TDTWriter {lvalue},RDKit::ROMol [,int=-1])"""
        ...
    @classmethod
    def __enter__(cls, RDKit) -> Any: ...
    @classmethod
    def __exit__(cls, type, value, traceback) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

def AddMetadataToPNGFile(self, metadata: dict, filename: AtomPairsParameters) -> object:
    """
    Adds metadata to PNG data read from a file.

    ARGUMENTS:

    metadata: dict with the metadata to be written(keys and values should be strings)

    filename: the PNG filename

    RETURNS:the updated PNG data

    C++ signature :boost::python::api::object AddMetadataToPNGFile(boost::python::dict,boost::python::api::object)
    """
    ...

def AddMetadataToPNGString(self, metadata: dict, png: AtomPairsParameters) -> object:
    """
    Adds metadata to a PNG string.

    ARGUMENTS:

    metadata: dict with the metadata to be written(keys and values should be strings)

    png: the PNG string

    RETURNS:the updated PNG data

    C++ signature :boost::python::api::object AddMetadataToPNGString(boost::python::dict,boost::python::api::object)
    """
    ...

def AtomFromSmarts(self, SMARTS: str) -> Atom:
    """
    Construct an atom from a SMARTS string

    C++ signature :RDKit::Atom* AtomFromSmarts(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
    ...

def AtomFromSmiles(self, SMILES: str) -> Atom:
    """
    Construct an atom from a SMILES string

    C++ signature :RDKit::Atom* AtomFromSmiles(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
    ...

def BondFromSmarts(self, SMILES: str) -> Bond:
    """
    Construct a bond from a SMARTS string

    C++ signature :RDKit::Bond* BondFromSmarts(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
    ...

def BondFromSmiles(self, SMILES: str) -> Bond:
    """
    Construct a bond from a SMILES string

    C++ signature :RDKit::Bond* BondFromSmiles(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
    ...

def CanonicalRankAtoms(
    self,
    mol: Mol,
    breakTies: bool = True,
    includeChirality: bool = True,
    includeIsotopes: bool = True,
) -> _vectj:
    """
    Returns the canonical atom ranking for each atom of a molecule fragment.If breakTies is False, this returns the symmetry class for each atom.  The symmetry
    class is used by the canonicalization routines to type each atom based on the whole
    chemistry of the molecular graph.  Any atom with the same rank (symmetry class) is
    indistinguishable.  For example:
    >>> mol = MolFromSmiles('C1NCN1')
    >>> list(CanonicalRankAtoms(mol, breakTies=False))
    [0,1,0,1]

    In this case the carbons have the same symmetry class and the nitrogens have the same
    symmetry class.  From the perspective of the Molecular Graph, they are identical.
    ARGUMENTS:

    mol: the molecule
    breakTies: (optional) force breaking of ranked ties [default=True]
    includeChirality: (optional) use chiral information when computing rank [default=True]
    includeIsotopes: (optional) use isotope information when computing rank [default=True]

    RETURNS:

    a string

    C++ signature :std::vector<unsigned int, std::allocator<unsigned int> > CanonicalRankAtoms(RDKit::ROMol [,bool=True [,bool=True [,bool=True]]])
    """
    ...

def CanonicalRankAtomsInFragment(
    self,
    mol: Mol,
    atomsToUse: AtomPairsParameters,
    bondsToUse: AtomPairsParameters = 0,
    atomSymbols: AtomPairsParameters = 0,
    breakTies: bool = True,
    includeChirality: bool = True,
    includeIsotopes: bool = True,
) -> _vecti:
    """
    Returns the canonical atom ranking for each atom of a molecule fragmentSee help(CanonicalRankAtoms) for more information.
    >>> mol = MolFromSmiles('C1NCN1.C1NCN1')
    >>> list(CanonicalRankAtomsInFragment(mol, atomsToUse=range(0,4), breakTies=False))
    [4,6,4,6,-1,-1,-1,-1]
    >>> list(CanonicalRankAtomsInFragment(mol, atomsToUse=range(4,8), breakTies=False))
    [-1,-1,-1,-1,4,6,4,6]

    ARGUMENTS:

    mol: the molecule
    atomsToUse : a list of atoms to include in the fragment
    bondsToUse : (optional) a list of bonds to include in the fragment
    if not provided, all bonds between the atoms provided
    will be included.
    atomSymbols : (optional) a list with the symbols to use for the atoms
    in the SMILES. This should have be mol.GetNumAtoms() long.
    breakTies: (optional) force breaking of ranked ties
    includeChirality: (optional) use chiral information when computing rank [default=True]
    includeIsotopes: (optional) use isotope information when computing rank [default=True]

    RETURNS:

    a string

    C++ signature :std::vector<int, std::allocator<int> > CanonicalRankAtomsInFragment(RDKit::ROMol,boost::python::api::object [,boost::python::api::object=0 [,boost::python::api::object=0 [,bool=True [,bool=True [,bool=True]]]]])
    """
    ...

def CanonicalizeEnhancedStereo(self, mol: Mol) -> None:
    """
    C++ signature :void CanonicalizeEnhancedStereo(RDKit::ROMol {lvalue})"""
    ...

def CreateAtomBoolPropertyList(
    self, mol: Mol, propName: str, missingValueMarker: str = "", lineSize: int = 190
) -> None:
    """
    creates a list property on the molecule from individual atom property values

    C++ signature :void CreateAtomBoolPropertyList(RDKit::ROMol {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’’ [,unsigned int=190]])
    """
    ...

def CreateAtomDoublePropertyList(
    self, mol: Mol, propName: str, missingValueMarker: str = "", lineSize: int = 190
) -> None:
    """
    creates a list property on the molecule from individual atom property values

    C++ signature :void CreateAtomDoublePropertyList(RDKit::ROMol {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’’ [,unsigned int=190]])
    """
    ...

def CreateAtomIntPropertyList(
    self, mol: Mol, propName: str, missingValueMarker: str = "", lineSize: int = 190
) -> None:
    """
    creates a list property on the molecule from individual atom property values

    C++ signature :void CreateAtomIntPropertyList(RDKit::ROMol {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’’ [,unsigned int=190]])
    """
    ...

def CreateAtomStringPropertyList(
    self, mol: Mol, propName: str, missingValueMarker: str = "", lineSize: int = 190
) -> None:
    """
    creates a list property on the molecule from individual atom property values

    C++ signature :void CreateAtomStringPropertyList(RDKit::ROMol {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’’ [,unsigned int=190]])
    """
    ...

def MetadataFromPNGFile(self, filename: AtomPairsParameters) -> dict:
    """
    Returns a dict with all metadata from the PNG file. Keys are strings, values are bytes.

    C++ signature :boost::python::dict MetadataFromPNGFile(boost::python::api::object)
    """
    ...

def MetadataFromPNGString(self, png: AtomPairsParameters) -> dict:
    """
    Returns a dict with all metadata from the PNG string. Keys are strings, values are bytes.

    C++ signature :boost::python::dict MetadataFromPNGString(boost::python::api::object)
    """
    ...

def MolFragmentToCXSmarts(
    self,
    mol: Mol,
    atomsToUse: AtomPairsParameters,
    bondsToUse: AtomPairsParameters = 0,
    isomericSmarts: bool = True,
) -> str:
    """
    Returns a SMARTS string for a fragment of a moleculeARGUMENTS:

    mol: the molecule
    atomsToUse: indices of atoms to include in the SMARTS string
    bondsToUse: indices of bonds to include in the SMARTS string (optional)
    isomericSmarts: (optional) include information about stereochemistry in
    the SMARTS.  Defaults to true.

    RETURNS:

    a string

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolFragmentToCXSmarts(RDKit::ROMol,boost::python::api::object [,boost::python::api::object=0 [,bool=True]])
    """
    ...

@overload
def MolFragmentToCXSmiles(
    self,
    mol: Mol,
    params: SmilesWriteParams,
    atomsToUse: AtomPairsParameters,
    bondsToUse: AtomPairsParameters = 0,
    atomSymbols: AtomPairsParameters = 0,
    bondSymbols: AtomPairsParameters = 0,
) -> str:
    """

    Returns the CXSMILES string for a fragment of a moleculeARGUMENTS:

    mol: the molecule
    atomsToUse : a list of atoms to include in the fragment
    bondsToUse : (optional) a list of bonds to include in the fragment
    if not provided, all bonds between the atoms provided
    will be included.
    atomSymbols : (optional) a list with the symbols to use for the atoms
    in the SMILES. This should have be mol.GetNumAtoms() long.
    bondSymbols : (optional) a list with the symbols to use for the bonds
    in the SMILES. This should have be mol.GetNumBonds() long.
    isomericSmiles: (optional) include information about stereochemistry in
    the SMILES.  Defaults to true.
    kekuleSmiles: (optional) use the Kekule form (no aromatic bonds) in
    the SMILES.  Defaults to false.
    rootedAtAtom: (optional) if non-negative, this forces the SMILES
    to start at a particular atom. Defaults to -1.
    canonical: (optional) if false no attempt will be made to canonicalize
    the molecule. Defaults to true.
    allBondsExplicit: (optional) if true, all bond orders will be explicitly indicated
    in the output SMILES. Defaults to false.
    allHsExplicit: (optional) if true, all H counts will be explicitly indicated
    in the output SMILES. Defaults to false.
    doRandom: (optional) if true, randomized the DFS transversal graph,
    so we can generate random smiles. Defaults to false.

    RETURNS:

    a string

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolFragmentToCXSmiles(RDKit::ROMol,boost::python::api::object [,boost::python::api::object=0 [,boost::python::api::object=0 [,boost::python::api::object=0 [,bool=True [,bool=False [,int=-1 [,bool=True [,bool=False [,bool=False]]]]]]]]])
    """
    ...

@overload
def MolFragmentToCXSmiles(
    self,
    mol: Mol,
    atomsToUse: AtomPairsParameters,
    bondsToUse: AtomPairsParameters = 0,
    atomSymbols: AtomPairsParameters = 0,
    bondSymbols: AtomPairsParameters = 0,
    isomericSmiles: bool = True,
    kekuleSmiles: bool = False,
    rootedAtAtom: int = -1,
    canonical: bool = True,
    allBondsExplicit: bool = False,
    allHsExplicit: bool = False,
) -> str: ...
def MolFragmentToSmarts(
    self,
    mol: Mol,
    atomsToUse: AtomPairsParameters,
    bondsToUse: AtomPairsParameters = 0,
    isomericSmarts: bool = True,
) -> str:
    """
    Returns a SMARTS string for a fragment of a moleculeARGUMENTS:

    mol: the molecule
    atomsToUse: indices of atoms to include in the SMARTS string
    bondsToUse: indices of bonds to include in the SMARTS string (optional)
    isomericSmarts: (optional) include information about stereochemistry in
    the SMARTS.  Defaults to true.

    RETURNS:

    a string

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolFragmentToSmarts(RDKit::ROMol,boost::python::api::object [,boost::python::api::object=0 [,bool=True]])
    """
    ...

@overload
def MolFragmentToSmiles(
    self,
    mol: Mol,
    params: SmilesWriteParams,
    atomsToUse: AtomPairsParameters,
    bondsToUse: AtomPairsParameters = 0,
    atomSymbols: AtomPairsParameters = 0,
    bondSymbols: AtomPairsParameters = 0,
) -> str:
    """

    Returns the canonical SMILES string for a fragment of a moleculeARGUMENTS:

    mol: the molecule
    atomsToUse : a list of atoms to include in the fragment
    bondsToUse : (optional) a list of bonds to include in the fragment
    if not provided, all bonds between the atoms provided
    will be included.
    atomSymbols : (optional) a list with the symbols to use for the atoms
    in the SMILES. This should have be mol.GetNumAtoms() long.
    bondSymbols : (optional) a list with the symbols to use for the bonds
    in the SMILES. This should have be mol.GetNumBonds() long.
    isomericSmiles: (optional) include information about stereochemistry in
    the SMILES.  Defaults to true.
    kekuleSmiles: (optional) use the Kekule form (no aromatic bonds) in
    the SMILES.  Defaults to false.
    rootedAtAtom: (optional) if non-negative, this forces the SMILES
    to start at a particular atom. Defaults to -1.
    canonical: (optional) if false no attempt will be made to canonicalize
    the molecule. Defaults to true.
    allBondsExplicit: (optional) if true, all bond orders will be explicitly indicated
    in the output SMILES. Defaults to false.
    allHsExplicit: (optional) if true, all H counts will be explicitly indicated
    in the output SMILES. Defaults to false.
    doRandom: (optional) if true, randomized the DFS transversal graph,
    so we can generate random smiles. Defaults to false.

    RETURNS:

    a string

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolFragmentToSmiles(RDKit::ROMol,boost::python::api::object [,boost::python::api::object=0 [,boost::python::api::object=0 [,boost::python::api::object=0 [,bool=True [,bool=False [,int=-1 [,bool=True [,bool=False [,bool=False]]]]]]]]])
    """
    ...

@overload
def MolFragmentToSmiles(
    self,
    mol: Mol,
    atomsToUse: AtomPairsParameters,
    bondsToUse: AtomPairsParameters = 0,
    atomSymbols: AtomPairsParameters = 0,
    bondSymbols: AtomPairsParameters = 0,
    isomericSmiles: bool = True,
    kekuleSmiles: bool = False,
    rootedAtAtom: int = -1,
    canonical: bool = True,
    allBondsExplicit: bool = False,
    allHsExplicit: bool = False,
) -> str: ...
def MolFromFASTA(
    self, text: AtomPairsParameters, sanitize: bool = True, flavor: int = 0
) -> Mol:
    """
    Construct a molecule from a FASTA string (currently only supports peptides).

    ARGUMENTS:

    text: string containing the FASTA
    sanitize: (optional) toggles sanitization of the molecule.
    Defaults to True.

    flavor: (optional)
    0 Protein, L amino acids (default)
    1 Protein, D amino acids
    2 RNA, no cap
    3 RNA, 5’ cap
    4 RNA, 3’ cap
    5 RNA, both caps
    6 DNA, no cap
    7 DNA, 5’ cap
    8 DNA, 3’ cap
    9 DNA, both caps

    RETURNS:

    a Mol object, None on failure.

    C++ signature :RDKit::ROMol* MolFromFASTA(boost::python::api::object [,bool=True [,int=0]])
    """
    ...

def MolFromHELM(self, text: AtomPairsParameters, sanitize: bool = True) -> Mol:
    """
    Construct a molecule from a HELM string (currently only supports peptides).

    ARGUMENTS:

    text: string containing the HELM
    sanitize: (optional) toggles sanitization of the molecule.
    Defaults to true.

    RETURNS:

    a Mol object, None on failure.

    C++ signature :RDKit::ROMol* MolFromHELM(boost::python::api::object [,bool=True])"""
    ...

def MolFromMol2Block(
    self,
    molBlock: str,
    sanitize: bool = True,
    removeHs: bool = True,
    cleanupSubstructures: bool = True,
) -> Mol:
    """
    Construct a molecule from a Tripos Mol2 block.

    NOTE:The parser expects the atom-typing scheme used by Corina.
    Atom types from Tripos’ dbtranslate are less supported.
    Other atom typing schemes are unlikely to work.

    ARGUMENTS:

    mol2Block: string containing the Mol2 block
    sanitize: (optional) toggles sanitization of the molecule.
    Defaults to True.
    removeHs: (optional) toggles removing hydrogens from the molecule.
    This only make sense when sanitization is done.
    Defaults to true.
    cleanupSubstructures: (optional) toggles standardizing some
    substructures found in mol2 files.
    Defaults to true.

    RETURNS:

    a Mol object, None on failure.

    C++ signature :RDKit::ROMol* MolFromMol2Block(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=True [,bool=True [,bool=True]]])
    """
    ...

def MolFromMol2File(
    self,
    molFileName: str,
    sanitize: bool = True,
    removeHs: bool = True,
    cleanupSubstructures: bool = True,
) -> Mol:
    """
    Construct a molecule from a Tripos Mol2 file.

    NOTE:The parser expects the atom-typing scheme used by Corina.
    Atom types from Tripos’ dbtranslate are less supported.
    Other atom typing schemes are unlikely to work.

    ARGUMENTS:

    fileName: name of the file to read
    sanitize: (optional) toggles sanitization of the molecule.
    Defaults to true.
    removeHs: (optional) toggles removing hydrogens from the molecule.
    This only make sense when sanitization is done.
    Defaults to true.
    cleanupSubstructures: (optional) toggles standardizing some
    substructures found in mol2 files.
    Defaults to true.

    RETURNS:

    a Mol object, None on failure.

    C++ signature :RDKit::ROMol* MolFromMol2File(char const* [,bool=True [,bool=True [,bool=True]]])
    """
    ...

@overload
def MolFromMolBlock(
    self,
    molBlock: AtomPairsParameters,
    sanitize: bool = True,
    removeHs: bool = True,
    strictParsing: bool = True,
) -> Mol:
    """
    Construct a molecule from a Mol block.

    ARGUMENTS:

    molBlock: string containing the Mol block
    sanitize: (optional) toggles sanitization of the molecule.
    Defaults to True.
    removeHs: (optional) toggles removing hydrogens from the molecule.
    This only make sense when sanitization is done.
    Defaults to true.
    strictParsing: (optional) if this is false, the parser is more lax about.
    correctness of the content.
    Defaults to true.

    RETURNS:

    a Mol object, None on failure.

    C++ signature :RDKit::ROMol* MolFromMolBlock(boost::python::api::object [,bool=True [,bool=True [,bool=True]]])
    """
    ...

@overload
def MolFromMolBlock(
    self,
    molBlock: AtomPairsParameters,
    sanitize: bool = True,
    removeHs: bool = True,
    strictParsing: bool = True,
) -> Mol: ...
@overload
def MolFromMolFile(
    self,
    molFileName: str,
    sanitize: bool = True,
    removeHs: bool = True,
    strictParsing: bool = True,
) -> Mol:
    """
    Construct a molecule from a Mol file.

    ARGUMENTS:

    fileName: name of the file to read
    sanitize: (optional) toggles sanitization of the molecule.
    Defaults to true.
    removeHs: (optional) toggles removing hydrogens from the molecule.
    This only make sense when sanitization is done.
    Defaults to true.
    strictParsing: (optional) if this is false, the parser is more lax about.
    correctness of the content.
    Defaults to true.

    RETURNS:

    a Mol object, None on failure.

    C++ signature :RDKit::ROMol* MolFromMolFile(char const* [,bool=True [,bool=True [,bool=True]]])
    """
    ...

@overload
def MolFromMolFile(
    self,
    molFileName: str,
    sanitize: bool = True,
    removeHs: bool = True,
    strictParsing: bool = True,
) -> Mol: ...
def MolFromPDBBlock(
    self,
    molBlock: AtomPairsParameters,
    sanitize: bool = True,
    removeHs: bool = True,
    flavor: int = 0,
    proximityBonding: bool = True,
) -> Mol:
    """
    Construct a molecule from a PDB block.

    ARGUMENTS:

    molBlock: string containing the PDB block
    sanitize: (optional) toggles sanitization of the molecule.
    Defaults to True.
    removeHs: (optional) toggles removing hydrogens from the molecule.
    This only make sense when sanitization is done.
    Defaults to true.
    flavor: (optional)
    proximityBonding: (optional) toggles automatic proximity bonding

    RETURNS:

    a Mol object, None on failure.

    C++ signature :RDKit::ROMol* MolFromPDBBlock(boost::python::api::object [,bool=True [,bool=True [,unsigned int=0 [,bool=True]]]])
    """
    ...

def MolFromPDBFile(
    self,
    molFileName: str,
    sanitize: bool = True,
    removeHs: bool = True,
    flavor: int = 0,
    proximityBonding: bool = True,
) -> Mol:
    """
    Construct a molecule from a PDB file.

    ARGUMENTS:

    fileName: name of the file to read
    sanitize: (optional) toggles sanitization of the molecule.
    Defaults to true.
    removeHs: (optional) toggles removing hydrogens from the molecule.
    This only make sense when sanitization is done.
    Defaults to true.
    flavor: (optional)
    proximityBonding: (optional) toggles automatic proximity bonding

    RETURNS:

    a Mol object, None on failure.

    C++ signature :RDKit::ROMol* MolFromPDBFile(char const* [,bool=True [,bool=True [,unsigned int=0 [,bool=True]]]])
    """
    ...

def MolFromPNGFile(self, filename: str, params: AtomPairsParameters = None) -> Mol:
    """
    Construct a molecule from metadata in a PNG file.

    ARGUMENTS:

    filename: the PNG filename
    params: used to provide optional parameters for the metadata parsing

    RETURNS:a Mol object, None on failure.

    C++ signature :RDKit::ROMol* MolFromPNGFile(char const* [,boost::python::api::object=None])
    """
    ...

def MolFromPNGString(
    self, png: AtomPairsParameters, params: AtomPairsParameters = None
) -> Mol:
    """
    Construct a molecule from metadata in a PNG string.

    ARGUMENTS:

    png: the PNG string
    params: used to provide optional parameters for the metadata parsing

    RETURNS:a Mol object, None on failure.

    C++ signature :RDKit::ROMol* MolFromPNGString(boost::python::api::object [,boost::python::api::object=None])
    """
    ...

def MolFromRDKitSVG(
    self, molBlock: AtomPairsParameters, sanitize: bool = True, removeHs: bool = True
) -> Mol:
    """
    Construct a molecule from an RDKit-generate SVG string.

    ARGUMENTS:

    svg: string containing the SVG data (must include molecule metadata)
    sanitize: (optional) toggles sanitization of the molecule.
    Defaults to True.
    removeHs: (optional) toggles removing hydrogens from the molecule.
    This only make sense when sanitization is done.
    Defaults to true.

    RETURNS:

    a Mol object, None on failure.

    NOTE: this functionality should be considered beta.

    C++ signature :RDKit::ROMol* MolFromRDKitSVG(boost::python::api::object [,bool=True [,bool=True]])
    """
    ...

def MolFromSequence(
    self, text: AtomPairsParameters, sanitize: bool = True, flavor: int = 0
) -> Mol:
    """
    Construct a molecule from a sequence string (currently only supports peptides).

    ARGUMENTS:

    text: string containing the sequence
    sanitize: (optional) toggles sanitization of the molecule.
    Defaults to True.

    flavor: (optional)
    0 Protein, L amino acids (default)
    1 Protein, D amino acids
    2 RNA, no cap
    3 RNA, 5’ cap
    4 RNA, 3’ cap
    5 RNA, both caps
    6 DNA, no cap
    7 DNA, 5’ cap
    8 DNA, 3’ cap
    9 DNA, both caps

    RETURNS:

    a Mol object, None on failure.

    C++ signature :RDKit::ROMol* MolFromSequence(boost::python::api::object [,bool=True [,int=0]])
    """
    ...

@overload
def MolFromSmarts(
    self, SMARTS: AtomPairsParameters, mergeHs: bool = False, replacements: dict = {}
) -> Mol:
    """
    Construct a molecule from a SMARTS string.

    ARGUMENTS:

    SMARTS: the smarts string
    params: used to provide optional parameters for the SMARTS parsing

    RETURNS:

    a Mol object, None on failure.

    C++ signature :RDKit::ROMol* MolFromSmarts(boost::python::api::object,RDKit::SmartsParserParams)
    """
    ...

@overload
def MolFromSmarts(
    self, SMARTS: AtomPairsParameters, params: SmartsParserParams
) -> Mol: ...
@overload
def MolFromSmiles(self, SMILES: AtomPairsParameters, params: SmilesParserParams) -> Mol:
    """
    Construct a molecule from a SMILES string.

    ARGUMENTS:

    SMILES: the smiles string
    sanitize: (optional) toggles sanitization of the molecule.
    Defaults to True.
    replacements: (optional) a dictionary of replacement strings (see below)
    Defaults to {}.

    RETURNS:

    a Mol object, None on failure.

    The optional replacements dict can be used to do string substitution of abbreviations
    in the input SMILES. The set of substitutions is repeatedly looped through until
    the string no longer changes. It is the responsibility of the caller to make sure
    that substitutions results in legal and sensible SMILES.
    Examples of replacements:

    CC{Q}C with {‘{Q}’:’OCCO’} -> CCOCCOC
    C{A}C{Q}C with {‘{Q}’:’OCCO’, ‘{A}’:’C1(CC1)’} -> CC1(CC1)COCCOC
    C{A}C{Q}C with {‘{Q}’:’{X}CC{X}’, ‘{A}’:’C1CC1’, ‘{X}’:’N’} -> CC1CC1CNCCNC

    C++ signature :RDKit::ROMol* MolFromSmiles(boost::python::api::object [,bool=True [,boost::python::dict={}]])
    """
    ...

@overload
def MolFromSmiles(
    self, SMILES: AtomPairsParameters, sanitize: bool = True, replacements: dict = {}
) -> Mol: ...
def MolFromTPLBlock(
    self,
    tplBlock: AtomPairsParameters,
    sanitize: bool = True,
    skipFirstConf: bool = False,
) -> Mol:
    """
    Construct a molecule from a TPL block.

    ARGUMENTS:

    fileName: name of the file to read
    sanitize: (optional) toggles sanitization of the molecule.
    Defaults to True.
    skipFirstConf: (optional) skips reading the first conformer.
    Defaults to False.
    This should be set to True when reading TPLs written by
    the CombiCode.

    RETURNS:

    a Mol object, None on failure.

    C++ signature :RDKit::ROMol* MolFromTPLBlock(boost::python::api::object [,bool=True [,bool=False]])
    """
    ...

def MolFromTPLFile(
    self, fileName: str, sanitize: bool = True, skipFirstConf: bool = False
) -> Mol:
    """
    Construct a molecule from a TPL file.

    ARGUMENTS:

    fileName: name of the file to read
    sanitize: (optional) toggles sanitization of the molecule.
    Defaults to True.
    skipFirstConf: (optional) skips reading the first conformer.
    Defaults to False.
    This should be set to True when reading TPLs written by
    the CombiCode.

    RETURNS:

    a Mol object, None on failure.

    C++ signature :RDKit::ROMol* MolFromTPLFile(char const* [,bool=True [,bool=False]])
    """
    ...

def MolFromXYZBlock(self, xyzFileName: AtomPairsParameters) -> Mol:
    """
    Construct a molecule from an XYZ string.

    ARGUMENTS:

    xyzBlock: the XYZ data to read

    RETURNS:

    a Mol object, None on failure.

    C++ signature :RDKit::ROMol* MolFromXYZBlock(boost::python::api::object)"""
    ...

def MolFromXYZFile(self, xyzFileName: str) -> Mol:
    """
    Construct a molecule from an XYZ file.

    ARGUMENTS:

    xyzname: name of the file to read

    RETURNS:

    a Mol object, None on failure.

    C++ signature :RDKit::ROMol* MolFromXYZFile(char const*)"""
    ...

def MolMetadataToPNGFile(
    self,
    mol: Mol,
    filename: AtomPairsParameters,
    includePkl: bool = True,
    includeSmiles: bool = True,
    includeMol: bool = False,
) -> object:
    """
    Adds molecular metadata to PNG data read from a file.

    ARGUMENTS:

    mol: the molecule
    filename: the PNG filename
    includePkl: include the RDKit’s internal binary format in the output
    includeSmiles: include CXSmiles in the output
    includeMol: include CTAB (Mol) in the output

    RETURNS:the updated PNG data

    C++ signature :boost::python::api::object MolMetadataToPNGFile(RDKit::ROMol,boost::python::api::object [,bool=True [,bool=True [,bool=False]]])
    """
    ...

def MolMetadataToPNGString(
    self,
    mol: Mol,
    png: AtomPairsParameters,
    includePkl: bool = True,
    includeSmiles: bool = True,
    includeMol: bool = False,
) -> object:
    """
    Adds molecular metadata to a PNG string.

    ARGUMENTS:

    mol: the molecule
    png: the PNG string
    includePkl: include the RDKit’s internal binary format in the output
    includeSmiles: include CXSmiles in the output
    includeMol: include CTAB (Mol) in the output

    RETURNS:the updated PNG data

    C++ signature :boost::python::api::object MolMetadataToPNGString(RDKit::ROMol,boost::python::api::object [,bool=True [,bool=True [,bool=False]]])
    """
    ...

def MolToCMLBlock(self, mol: Mol, confId: int = -1, kekulize: bool = True) -> str:
    """
    Writes a CML block for a moleculeARGUMENTS:

    mol: the molecule
    confId: (optional) selects which conformation to output
    kekulize: (optional) triggers kekulization of the molecule before it’s written

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToCMLBlock(RDKit::ROMol [,int=-1 [,bool=True]])
    """
    ...

def MolToCMLFile(
    self, mol: Mol, filename: str, confId: int = -1, kekulize: bool = True
) -> None:
    """
    Writes a CML file for a moleculeARGUMENTS:

    mol: the molecule
    filename: the file to write to
    confId: (optional) selects which conformation to output
    kekulize: (optional) triggers kekulization of the molecule before it’s written

    C++ signature :void MolToCMLFile(RDKit::ROMol,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,int=-1 [,bool=True]])
    """
    ...

def MolToCXSmarts(self, mol: Mol, isomericSmiles: bool = True) -> str:
    """
    Returns a SMARTS string for a moleculeARGUMENTS:

    mol: the molecule
    isomericSmiles: (optional) include information about stereochemistry in
    the SMARTS.  Defaults to true.

    RETURNS:

    a string

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToCXSmarts(RDKit::ROMol [,bool=True])
    """
    ...

@overload
def MolToCXSmiles(
    self, mol: Mol, params: SmilesWriteParams, flags: int = CXSmilesFields.CX_ALL
) -> str:
    """

    Returns the CXSMILES string for a moleculeARGUMENTS:

    mol: the molecule
    isomericSmiles: (optional) include information about stereochemistry in
    the SMILES.  Defaults to true.
    kekuleSmiles: (optional) use the Kekule form (no aromatic bonds) in
    the SMILES.  Defaults to false.
    rootedAtAtom: (optional) if non-negative, this forces the SMILES
    to start at a particular atom. Defaults to -1.
    canonical: (optional) if false no attempt will be made to canonicalize
    the molecule. Defaults to true.
    allBondsExplicit: (optional) if true, all bond orders will be explicitly indicated
    in the output SMILES. Defaults to false.
    allHsExplicit: (optional) if true, all H counts will be explicitly indicated
    in the output SMILES. Defaults to false.

    RETURNS:

    a string

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToCXSmiles(RDKit::ROMol [,bool=True [,bool=False [,int=-1 [,bool=True [,bool=False [,bool=False [,bool=False]]]]]]])
    """
    ...

@overload
def MolToCXSmiles(
    self,
    mol: Mol,
    isomericSmiles: bool = True,
    kekuleSmiles: bool = False,
    rootedAtAtom: int = -1,
    canonical: bool = True,
    allBondsExplicit: bool = False,
    allHsExplicit: bool = False,
    doRandom: bool = False,
) -> str: ...
def MolToFASTA(self, mol: Mol) -> str:
    """
    Returns the FASTA string for a moleculeARGUMENTS:

    mol: the molecule

    NOTE: the molecule should contain monomer information in AtomMonomerInfo structures
    RETURNS:

    a string

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToFASTA(RDKit::ROMol)
    """
    ...

def MolToHELM(self, mol: Mol) -> str:
    """
    Returns the HELM string for a moleculeARGUMENTS:

    mol: the molecule

    NOTE: the molecule should contain monomer information in AtomMonomerInfo structures
    RETURNS:

    a string

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToHELM(RDKit::ROMol)
    """
    ...

def MolToMolBlock(
    self,
    mol: Mol,
    includeStereo: bool = True,
    confId: int = -1,
    kekulize: bool = True,
    forceV3000: bool = False,
) -> str:
    """
    Returns a Mol block for a moleculeARGUMENTS:

    mol: the molecule
    includeStereo: (optional) toggles inclusion of stereochemical
    information in the output
    confId: (optional) selects which conformation to output (-1 = default)
    kekulize: (optional) triggers kekulization of the molecule before it’s written,
    as suggested by the MDL spec.
    forceV3000 (optional) force generation a V3000 mol block (happens automatically with
    more than 999 atoms or bonds)

    RETURNS:

    a string

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToMolBlock(RDKit::ROMol [,bool=True [,int=-1 [,bool=True [,bool=False]]]])
    """
    ...

def MolToMolFile(
    self,
    mol: Mol,
    filename: str,
    includeStereo: bool = True,
    confId: int = -1,
    kekulize: bool = True,
    forceV3000: bool = False,
) -> None:
    """
    Writes a Mol file for a moleculeARGUMENTS:

    mol: the molecule
    filename: the file to write to
    includeStereo: (optional) toggles inclusion of stereochemical
    information in the output
    confId: (optional) selects which conformation to output (-1 = default)
    kekulize: (optional) triggers kekulization of the molecule before it’s written,
    as suggested by the MDL spec.
    forceV3000 (optional) force generation a V3000 mol block (happens automatically with
    more than 999 atoms or bonds)

    RETURNS:

    a string

    C++ signature :void MolToMolFile(RDKit::ROMol,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=True [,int=-1 [,bool=True [,bool=False]]]])
    """
    ...

def MolToPDBBlock(self, mol: Mol, confId: int = -1, flavor: int = 0) -> str:
    """
    Returns a PDB block for a moleculeARGUMENTS:

    mol: the molecule
    confId: (optional) selects which conformation to output (-1 = default)

    flavor: (optional)
    flavor & 1 : Write MODEL/ENDMDL lines around each record
    flavor & 2 : Don’t write any CONECT records
    flavor & 4 : Write CONECT records in both directions
    flavor & 8 : Don’t use multiple CONECTs to encode bond order
    flavor & 16 : Write MASTER record
    flavor & 32 : Write TER record

    RETURNS:

    a string

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToPDBBlock(RDKit::ROMol [,int=-1 [,unsigned int=0]])
    """
    ...

def MolToPDBFile(
    self, mol: Mol, filename: str, confId: int = -1, flavor: int = 0
) -> None:
    """
    Writes a PDB file for a moleculeARGUMENTS:

    mol: the molecule
    filename: name of the file to write
    confId: (optional) selects which conformation to output (-1 = default)

    flavor: (optional)
    flavor & 1 : Write MODEL/ENDMDL lines around each record
    flavor & 2 : Don’t write any CONECT records
    flavor & 4 : Write CONECT records in both directions
    flavor & 8 : Don’t use multiple CONECTs to encode bond order
    flavor & 16 : Write MASTER record
    flavor & 32 : Write TER record

    RETURNS:

    a string

    C++ signature :void MolToPDBFile(RDKit::ROMol,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,int=-1 [,unsigned int=0]])
    """
    ...

def MolToRandomSmilesVect(
    self,
    mol: Mol,
    numSmiles: int,
    randomSeed: int = 0,
    isomericSmiles: bool = True,
    kekuleSmiles: bool = False,
    allBondsExplicit: bool = False,
    allHsExplicit: bool = False,
) -> list:
    """
    returns a list of SMILES generated using the randomSmiles algorithm

    C++ signature :boost::python::list MolToRandomSmilesVect(RDKit::ROMol,unsigned int [,unsigned int=0 [,bool=True [,bool=False [,bool=False [,bool=False]]]]])
    """
    ...

def MolToSequence(self, mol: Mol) -> str:
    """
    Returns the sequence string for a moleculeARGUMENTS:

    mol: the molecule

    NOTE: the molecule should contain monomer information in AtomMonomerInfo structures
    RETURNS:

    a string

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToSequence(RDKit::ROMol)
    """
    ...

def MolToSmarts(self, mol: Mol, isomericSmiles: bool = True) -> str:
    """
    Returns a SMARTS string for a moleculeARGUMENTS:

    mol: the molecule
    isomericSmiles: (optional) include information about stereochemistry in
    the SMARTS.  Defaults to true.

    RETURNS:

    a string

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToSmarts(RDKit::ROMol [,bool=True])
    """
    ...

@overload
def MolToSmiles(self, mol: Mol, params: SmilesWriteParams) -> str:
    """

    Returns the canonical SMILES string for a moleculeARGUMENTS:

    mol: the molecule
    isomericSmiles: (optional) include information about stereochemistry in
    the SMILES.  Defaults to true.
    kekuleSmiles: (optional) use the Kekule form (no aromatic bonds) in
    the SMILES.  Defaults to false.
    rootedAtAtom: (optional) if non-negative, this forces the SMILES
    to start at a particular atom. Defaults to -1.
    canonical: (optional) if false no attempt will be made to canonicalize
    the molecule. Defaults to true.
    allBondsExplicit: (optional) if true, all bond orders will be explicitly indicated
    in the output SMILES. Defaults to false.
    allHsExplicit: (optional) if true, all H counts will be explicitly indicated
    in the output SMILES. Defaults to false.

    RETURNS:

    a string

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToSmiles(RDKit::ROMol [,bool=True [,bool=False [,int=-1 [,bool=True [,bool=False [,bool=False [,bool=False]]]]]]])
    """
    ...

@overload
def MolToSmiles(
    self,
    mol: Mol,
    isomericSmiles: bool = True,
    kekuleSmiles: bool = False,
    rootedAtAtom: int = -1,
    canonical: bool = True,
    allBondsExplicit: bool = False,
    allHsExplicit: bool = False,
    doRandom: bool = False,
) -> str: ...
def MolToTPLBlock(
    self,
    mol: Mol,
    partialChargeProp: str = "_GasteigerCharge",
    writeFirstConfTwice: bool = False,
) -> str:
    """
    Returns the Tpl block for a molecule.

    ARGUMENTS:

    mol: the molecule
    partialChargeProp: name of the property to use for partial charges
    Defaults to ‘_GasteigerCharge’.
    writeFirstConfTwice: Defaults to False.
    This should be set to True when writing TPLs to be read by
    the CombiCode.

    RETURNS:

    a string

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToTPLBlock(RDKit::ROMol [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’_GasteigerCharge’ [,bool=False]])
    """
    ...

def MolToTPLFile(
    self,
    mol: Mol,
    fileName: str,
    partialChargeProp: str = "_GasteigerCharge",
    writeFirstConfTwice: bool = False,
) -> None:
    """
    Writes a molecule to a TPL file.

    ARGUMENTS:

    mol: the molecule
    fileName: name of the file to write
    partialChargeProp: name of the property to use for partial charges
    Defaults to ‘_GasteigerCharge’.
    writeFirstConfTwice: Defaults to False.
    This should be set to True when writing TPLs to be read by
    the CombiCode.

    C++ signature :void MolToTPLFile(RDKit::ROMol,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’_GasteigerCharge’ [,bool=False]])
    """
    ...

def MolToV3KMolBlock(
    self, mol: Mol, includeStereo: bool = True, confId: int = -1, kekulize: bool = True
) -> str:
    """
    Returns a V3000 Mol block for a moleculeARGUMENTS:

    mol: the molecule
    includeStereo: (optional) toggles inclusion of stereochemical
    information in the output
    confId: (optional) selects which conformation to output (-1 = default)
    kekulize: (optional) triggers kekulization of the molecule before it’s written,
    as suggested by the MDL spec.

    RETURNS:

    a string

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToV3KMolBlock(RDKit::ROMol [,bool=True [,int=-1 [,bool=True]]])
    """
    ...

def MolToV3KMolFile(
    self,
    mol: Mol,
    filename: str,
    includeStereo: bool = True,
    confId: int = -1,
    kekulize: bool = True,
) -> None:
    """
    Writes a V3000 Mol file for a moleculeARGUMENTS:

    mol: the molecule
    filename: the file to write to
    includeStereo: (optional) toggles inclusion of stereochemical
    information in the output
    confId: (optional) selects which conformation to output (-1 = default)
    kekulize: (optional) triggers kekulization of the molecule before it’s written,
    as suggested by the MDL spec.

    RETURNS:

    a string

    C++ signature :void MolToV3KMolFile(RDKit::ROMol,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=True [,int=-1 [,bool=True]]])
    """
    ...

def MolToXYZBlock(self, mol: Mol, confId: int = -1) -> str:
    """
    Returns a XYZ block for a moleculeARGUMENTS:

    mol: the molecule
    confId: (optional) selects which conformation to output (-1 = default)

    RETURNS:

    a string

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToXYZBlock(RDKit::ROMol [,int=-1])
    """
    ...

def MolToXYZFile(self, mol: Mol, filename: str, confId: int = -1) -> None:
    """
    Writes a XYZ file for a moleculeARGUMENTS:

    mol: the molecule
    filename: the file to write to
    confId: (optional) selects which conformation to output (-1 = default)

    C++ signature :void MolToXYZFile(RDKit::ROMol,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,int=-1])
    """
    ...

def MolsFromCDXML(
    self, cdxml: AtomPairsParameters, sanitize: bool = True, removeHs: bool = True
) -> tuple:
    """
    Construct a molecule from a cdxml string.

    Note that the CDXML format is large and complex, the RDKit doesn’t support
    full functionality, just the base ones required for molecule and
    reaction parsing.
    ARGUMENTS:

    filename: the cdxml string
    sanitize: if True, sanitize the molecules [default True]
    removeHs: if True, convert explicit Hs into implicit Hs. [default True]

    RETURNS:an iterator of parsed Mol objects.

    C++ signature :boost::python::tuple MolsFromCDXML(boost::python::api::object [,bool=True [,bool=True]])
    """
    ...

def MolsFromCDXMLFile(
    self, filename: str, sanitize: bool = True, removeHs: bool = True
) -> object:
    """
    Construct a molecule from a cdxml file.

    Note that the CDXML format is large and complex, the RDKit doesn’t support
    full functionality, just the base ones required for molecule and
    reaction parsing.
    ARGUMENTS:

    filename: the cdxml filename
    sanitize: if True, sanitize the molecules [default True]
    removeHs: if True, convert explicit Hs into implicit Hs. [default True]

    RETURNS:an iterator of parsed Mol objects.

    C++ signature :boost::python::api::object MolsFromCDXMLFile(char const* [,bool=True [,bool=True]])
    """
    ...

def MolsFromPNGFile(
    self, filename: str, tag: str = "rdkitPKL", params: AtomPairsParameters = None
) -> object:
    """
    returns a tuple of molecules constructed from the PNG file

    C++ signature :boost::python::api::object MolsFromPNGFile(char const* [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’rdkitPKL’ [,boost::python::api::object=None]])
    """
    ...

def MolsFromPNGString(
    self,
    png: AtomPairsParameters,
    tag: str = "rdkitPKL",
    params: AtomPairsParameters = None,
) -> tuple:
    """
    returns a tuple of molecules constructed from the PNG string

    C++ signature :boost::python::tuple MolsFromPNGString(boost::python::api::object [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’rdkitPKL’ [,boost::python::api::object=None]])
    """
    ...

def SmilesMolSupplierFromText(
    self,
    text: str,
    delimiter: str = " ",
    smilesColumn: int = 0,
    nameColumn: int = 1,
    titleLine: bool = True,
    sanitize: bool = True,
) -> SmilesMolSupplier:
    """
    C++ signature :RDKit::SmilesMolSupplier* SmilesMolSupplierFromText(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’ ‘ [,int=0 [,int=1 [,bool=True [,bool=True]]]]])
    """
    ...
