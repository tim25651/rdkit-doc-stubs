"""
rdkit.Dbase.DbResultSet module¶
defines class _DbResultSet_ for lazy interactions with Db query results
Note

this uses the Python iterator interface, so you’ll need python 2.2 or above.
"""
from _typeshed import Incomplete
from rdkit.Dbase import DbInfo as DbInfo

class DbResultBase(object):
    cursor: Incomplete
    removeDups: Incomplete
    transform: Incomplete
    cmd: Incomplete
    conn: Incomplete
    extras: Incomplete

    def __init__(
        self,
        cursor,
        conn,
        cmd,
        removeDups: int = ...,
        transform: Incomplete | None = ...,
        extras: Incomplete | None = ...,
    ) -> None: ...
    def Reset(self):
        """
        implement in subclasses"""
        ...
    def __iter__(self): ...
    def GetColumnNames(self): ...
    def GetColumnNamesAndTypes(self): ...
    def GetColumnTypes(self): ...
    def Reset(self):
        """
        implement in subclasses"""
        ...
    def GetColumnNamesAndTypes(self): ...

class DbResultSet(DbResultBase):
    """
    Only supports forward iteration"""

    seen: Incomplete

    def __init__(self, *args, **kwargs) -> None: ...
    def Reset(self):
        """
        implement in subclasses"""
        ...
    def next(self): ...
    __next__ = next

class RandomAccessDbResultSet(DbResultBase):
    """
    Supports random access"""

    results: Incomplete
    seen: Incomplete

    def __init__(self, *args, **kwargs) -> None: ...
    def Reset(self):
        """
        implement in subclasses"""
        ...
    cursor: Incomplete

    def __getitem__(self, idx): ...
    def __len__(self) -> int: ...
    def next(self): ...
    __next__ = next
