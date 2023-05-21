"""
rdkit.Dbase.DbModule module¶
"""
from _typeshed import Incomplete
from pyPgSQL.PgSQL import *
from rdkit import RDConfig as RDConfig

sqlTextTypes: Incomplete
sqlIntTypes: Incomplete
sqlFloatTypes: Incomplete
sqlBinTypes: Incomplete
getTablesSql: str
getTablesAndViewsSql: str
getDbSql: str
fileWildcard = dbFileWildcard
placeHolder: str
binaryTypeName: str
binaryHolder = memoryview
RDTestDatabase: str
dbFileWildcard: str
fileWildcard = dbFileWildcard

def connect(self, x, *args): ...
