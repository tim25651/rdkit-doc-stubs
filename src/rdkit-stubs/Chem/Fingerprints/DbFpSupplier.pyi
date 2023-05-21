"""
rdkit.Chem.Fingerprints.DbFpSupplier module¶
Supplies a class for working with fingerprints from databases
#DOC
"""
from _typeshed import Incomplete
from rdkit import DataStructs as DataStructs
from rdkit.VLib import Node
from rdkit.VLib.Node import VLibNode as VLibNode

class DbFpSupplier(Node.VLibNode):
    """
    new fps come back with all additional fields from the
    database set in a “_fieldsFromDb” data member
    DbResults should be a subclass of Dbase.DbResultSet.DbResultBase"""

    fpCol: Incomplete

    def __init__(
        self, dbResults, fpColName: str = ..., usePickles: bool = ...
    ) -> None: ...
    def GetColumnNames(self): ...
    def next(self):
        """
        part of the iterator interface
        raises StopIteration on failure"""
        ...
    __next__ = next

class ForwardDbFpSupplier(DbFpSupplier):
    """
    DbFp supplier supporting only forward iteration
    >>> from rdkit import RDConfig
    >>> from rdkit.Dbase.DbConnection import DbConnect
    >>> fName = RDConfig.RDTestDatabase
    >>> conn = DbConnect(fName,'simple_combined')
    >>> suppl = ForwardDbFpSupplier(conn.GetData())

    we can loop over the supplied fingerprints:
    >>> fps = []
    >>> for fp in suppl:
    ...   fps.append(fp)
    >>> len(fps)
    12

    DbResults should be a subclass of Dbase.DbResultSet.DbResultBase"""

    def NextItem(self):
        """
        NOTE: this has side effects"""
        ...
    def __init__(self, *args, **kwargs) -> None: ...
    def reset(self):
        """
        resets our iteration state"""
        ...
    def NextItem(self):
        """
        NOTE: this has side effects"""
        ...

class RandomAccessDbFpSupplier(DbFpSupplier):
    """
    DbFp supplier supporting random access:
    >>> import os.path
    >>> from rdkit import RDConfig
    >>> from rdkit.Dbase.DbConnection import DbConnect
    >>> fName = RDConfig.RDTestDatabase
    >>> conn = DbConnect(fName,'simple_combined')
    >>> suppl = RandomAccessDbFpSupplier(conn.GetData())
    >>> len(suppl)
    12

    we can pull individual fingerprints:
    >>> fp = suppl[5]
    >>> fp.GetNumBits()
    128
    >>> fp.GetNumOnBits()
    54

    a standard loop over the fingerprints:
    >>> fps = []
    >>> for fp in suppl:
    ...   fps.append(fp)
    >>> len(fps)
    12

    or we can use an indexed loop:
    >>> fps = [None] * len(suppl)
    >>> for i in range(len(suppl)):
    ...   fps[i] = suppl[i]
    >>> len(fps)
    12

    DbResults should be a subclass of Dbase.DbResultSet.DbResultBase"""

    def NextItem(self): ...
    def __init__(self, *args, **kwargs) -> None: ...
    def __len__(self) -> int: ...
    def __getitem__(self, idx): ...
    def reset(self):
        """
        resets our iteration state"""
        ...
    def NextItem(self): ...
