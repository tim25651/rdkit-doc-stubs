"""
rdkit.VLib.NodeLib.DbPickleSupplier module¶
"""
from _typeshed import Incomplete
from rdkit import RDConfig as RDConfig
from rdkit.VLib import Supply
from rdkit.VLib.Supply import SupplyNode as SupplyNode

class _lazyDataSeq:
    cursor: Incomplete
    cmd: Incomplete

    def __init__(
        self,
        cursor,
        cmd,
        pickleCol: int = ...,
        depickle: int = ...,
        klass: Incomplete | None = ...,
    ) -> None: ...
    def __iter__(self): ...
    def next(self): ...

class _dataSeq(_lazyDataSeq):
    cursor: Incomplete
    cmd: Incomplete
    res: Incomplete
    rowCount: int
    idx: int

    def __init__(
        self, cursor, cmd, pickleCol: int = ..., depickle: int = ...
    ) -> None: ...
    def __iter__(self): ...
    def next(self): ...
    def __len__(self) -> int: ...
    def __getitem__(self, idx): ...

class DbPickleSupplyNode(Supply.SupplyNode):
    """
    Supplies pickled objects from a db result set:

    Sample Usage:>>> from rdkit.Dbase.DbConnection import DbConnect"""

    def next(self): ...
    def __init__(self, cursor, cmd, binaryCol, **kwargs) -> None: ...
    def reset(self):
        """
        resets our iteration state"""
        ...
    def next(self): ...

def GetNode(self, dbName, tableName): ...
