"""
rdkit.Chem.Suppliers.DbMolSupplier module¶
Supplies a class for working with molecules from databases
"""
import _io
from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit.Chem.Suppliers import MolSupplier
from rdkit.Chem.Suppliers.MolSupplier import MolSupplier as MolSupplier

def warning(self, msg, dest: _io.TextIOWrapper = ...): ...

class DbMolSupplier(MolSupplier.MolSupplier):
    """
    new molecules come back with all additional fields from the
    database set in a “_fieldsFromDb” data member
    DbResults should be a subclass of Dbase.DbResultSet.DbResultBase"""

    molCol: int
    transformFunc: Incomplete
    nameCol: Incomplete
    molFmt: Incomplete

    def __init__(
        self,
        dbResults,
        molColumnFormats=...,
        nameCol: str = ...,
        transformFunc: Incomplete | None = ...,
        **kwargs,
    ) -> None: ...
    def GetColumnNames(self): ...

class ForwardDbMolSupplier(DbMolSupplier):
    """
    DbMol supplier supporting only forward iteration
    new molecules come back with all additional fields from the
    database set in a “_fieldsFromDb” data member
    DbResults should be an iterator for Dbase.DbResultSet.DbResultBase"""

    def NextMol(self):
        """
        NOTE: this has side effects"""
        ...
    def __init__(self, dbResults, **kwargs) -> None: ...
    def Reset(self): ...
    def NextMol(self):
        """
        NOTE: this has side effects"""
        ...

class RandomAccessDbMolSupplier(DbMolSupplier):
    """
    DbResults should be a Dbase.DbResultSet.RandomAccessDbResultSet"""

    def NextMol(self):
        """
        Must be implemented in child class"""
        ...
    def __init__(self, dbResults, **kwargs) -> None: ...
    def __len__(self) -> int: ...
    def __getitem__(self, idx): ...
    def Reset(self): ...
    def NextMol(self):
        """
        Must be implemented in child class"""
        ...

def warning(self, msg, dest: _io.TextIOWrapper = ...): ...
