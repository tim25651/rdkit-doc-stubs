"""
rdkit.Dbase.StorageUtils module¶
Various storage (molecular and otherwise) functionality
"""
from _typeshed import Incomplete
from rdkit import RDConfig as RDConfig
from rdkit.Dbase import DbModule as DbModule

def ValidateRDId(self, ID):
    """
    returns whether or not an RDId is valid
    >>> ValidateRDId('RDCmpd-000-009-9')
    1
    >>> ValidateRDId('RDCmpd-009-000-009-8')
    1
    >>> ValidateRDId('RDCmpd-009-000-109-8')
    0
    >>> ValidateRDId('bogus')
    0"""
    ...

def RDIdToInt(self, ID, validate=1):
    """
    Returns the integer index for a given RDId
    Throws a ValueError on error
    >>> RDIdToInt('RDCmpd-000-009-9')
    9
    >>> RDIdToInt('RDCmpd-009-000-009-8')
    9000009
    >>> RDIdToInt('RDData_000_009_9')
    9
    >>> try:
    ...   RDIdToInt('RDCmpd-009-000-109-8')
    ... except ValueError:
    ...   print('ok')
    ... else:
    ...   print('failed')
    ok
    >>> try:
    ...   RDIdToInt('bogus')
    ... except ValueError:
    ...   print('ok')
    ... else:
    ...   print('failed')
    ok"""
    ...

def IndexToRDId(self, idx, leadText="RDCmpd"):
    """
    Converts an integer index into an RDId

    The format of the ID is:leadText-xxx-xxx-xxx-y

    The number blocks are zero padded and the final digit (y)
    is a checksum:
    >>> str(IndexToRDId(9))
    'RDCmpd-000-009-9'
    >>> str(IndexToRDId(9009))
    'RDCmpd-009-009-8'

    A millions block is included if it’s nonzero:
    >>> str(IndexToRDId(9000009))
    'RDCmpd-009-000-009-8'

    The text at the beginning can be altered:
    >>> str(IndexToRDId(9,leadText='RDAlt'))
    'RDAlt-000-009-9'

    Negative indices are errors:
    >>> try:
    ...   IndexToRDId(-1)
    ... except ValueError:
    ...   print('ok')
    ... else:
    ...   print('failed')
    ok"""
    ...

def GetNextId(self, conn, table, idColName="Id"):
    """
    returns the next available Id in the database
    see RegisterItem for testing/documentation"""
    ...

def GetNextRDId(self, conn, table, idColName="Id", leadText=""):
    """
    returns the next available RDId in the database
    see RegisterItem for testing/documentation"""
    ...

def IndexToRDId(self, idx, leadText="RDCmpd"):
    """
    Converts an integer index into an RDId

    The format of the ID is:leadText-xxx-xxx-xxx-y

    The number blocks are zero padded and the final digit (y)
    is a checksum:
    >>> str(IndexToRDId(9))
    'RDCmpd-000-009-9'
    >>> str(IndexToRDId(9009))
    'RDCmpd-009-009-8'

    A millions block is included if it’s nonzero:
    >>> str(IndexToRDId(9000009))
    'RDCmpd-009-000-009-8'

    The text at the beginning can be altered:
    >>> str(IndexToRDId(9,leadText='RDAlt'))
    'RDAlt-000-009-9'

    Negative indices are errors:
    >>> try:
    ...   IndexToRDId(-1)
    ... except ValueError:
    ...   print('ok')
    ... else:
    ...   print('failed')
    ok"""
    ...

def RDIdToInt(self, ID, validate=1):
    """
    Returns the integer index for a given RDId
    Throws a ValueError on error
    >>> RDIdToInt('RDCmpd-000-009-9')
    9
    >>> RDIdToInt('RDCmpd-009-000-009-8')
    9000009
    >>> RDIdToInt('RDData_000_009_9')
    9
    >>> try:
    ...   RDIdToInt('RDCmpd-009-000-109-8')
    ... except ValueError:
    ...   print('ok')
    ... else:
    ...   print('failed')
    ok
    >>> try:
    ...   RDIdToInt('bogus')
    ... except ValueError:
    ...   print('ok')
    ... else:
    ...   print('failed')
    ok"""
    ...

def RegisterItem(
    self,
    conn,
    table,
    value,
    columnName,
    data=None,
    id="",
    idColName="Id",
    leadText="RDCmpd",
):
    """
    >>> from rdkit.Dbase.DbConnection import DbConnect
    >>> conn = DbConnect(tempDbName)
    >>> tblName = 'StorageTest'
    >>> conn.AddTable(tblName,'id varchar(32) not null primary key,label varchar(40),val int')
    >>> RegisterItem(conn,tblName,'label1','label',['label1',1])==(1, 'RDCmpd-000-001-1')
    True
    >>> RegisterItem(conn,tblName,'label2','label',['label2',1])==(1, 'RDCmpd-000-002-2')
    True
    >>> RegisterItem(conn,tblName,'label1','label',['label1',1])==(0, 'RDCmpd-000-001-1')
    True
    >>> str(GetNextRDId(conn,tblName))
    'RDCmpd-000-003-3'
    >>> tuple(conn.GetData(table=tblName)[0])==('RDCmpd-000-001-1', 'label1', 1)
    True

    It’s also possible to provide ids by hand:
    >>> RegisterItem(conn,tblName,'label10','label',['label10',1],
    ...              id='RDCmpd-000-010-1')==(1, 'RDCmpd-000-010-1')
    True
    >>> str(GetNextRDId(conn,tblName))
    'RDCmpd-000-011-2'"""
    ...

def RegisterItems(
    self,
    conn,
    table,
    values,
    columnName,
    rows,
    startId="",
    idColName="Id",
    leadText="RDCmpd",
): ...
def ValidateRDId(self, ID):
    """
    returns whether or not an RDId is valid
    >>> ValidateRDId('RDCmpd-000-009-9')
    1
    >>> ValidateRDId('RDCmpd-009-000-009-8')
    1
    >>> ValidateRDId('RDCmpd-009-000-109-8')
    0
    >>> ValidateRDId('bogus')
    0"""
    ...

__test__: Incomplete