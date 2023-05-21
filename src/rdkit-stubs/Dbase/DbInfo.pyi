"""
rdkit.Dbase.DbInfo module¶
"""
from _typeshed import Incomplete
from rdkit import RDConfig as RDConfig
from rdkit.Dbase import DbModule as DbModule

def GetColumnInfoFromCursor(self, cursor): ...
def GetColumnNames(
    self, dBase, table, user="sysdba", password="masterkey", join="", what="*", cn=None
):
    """
    gets a list of columns available in a DB table
    Arguments

    dBase: the name of the DB file to be used
    table: the name of the table to query
    user: the username for DB access
    password: the password to be used for DB access
    join: an optional join clause  (omit the verb ‘join’)
    what: an optional clause indicating what to select

    Returns

    a list of column names"""
    ...

def GetColumnNamesAndTypes(
    self, dBase, table, user="sysdba", password="masterkey", join="", what="*", cn=None
):
    """
    gets a list of columns available in a DB table along with their types
    Arguments

    dBase: the name of the DB file to be used
    table: the name of the table to query
    user: the username for DB access
    password: the password to be used for DB access
    join: an optional join clause (omit the verb ‘join’)
    what: an optional clause indicating what to select

    Returns

    a list of 2-tuples containing:

    column name
    column type"""
    ...

sqlTextTypes: Incomplete
sqlIntTypes: Incomplete
sqlFloatTypes: Incomplete
sqlBinTypes: Incomplete

def GetDbNames(
    self, user="sysdba", password="masterkey", dirName=".", dBase="::template1", cn=None
):
    """
    returns a list of databases that are available
    Arguments

    user: the username for DB access
    password: the password to be used for DB access

    Returns

    a list of db names (strings)"""
    ...

def GetTableNames(
    self, dBase, user="sysdba", password="masterkey", includeViews=0, cn=None
):
    """
    returns a list of tables available in a database
    Arguments

    dBase: the name of the DB file to be used
    user: the username for DB access
    password: the password to be used for DB access
    includeViews: if this is non-null, the views in the db will
    also be returned

    Returns

    a list of table names (strings)"""
    ...

def GetColumnInfoFromCursor(self, cursor): ...
def GetColumnNamesAndTypes(
    self, dBase, table, user="sysdba", password="masterkey", join="", what="*", cn=None
):
    """
    gets a list of columns available in a DB table along with their types
    Arguments

    dBase: the name of the DB file to be used
    table: the name of the table to query
    user: the username for DB access
    password: the password to be used for DB access
    join: an optional join clause (omit the verb ‘join’)
    what: an optional clause indicating what to select

    Returns

    a list of 2-tuples containing:

    column name
    column type"""
    ...

def GetColumnNames(
    self, dBase, table, user="sysdba", password="masterkey", join="", what="*", cn=None
):
    """
    gets a list of columns available in a DB table
    Arguments

    dBase: the name of the DB file to be used
    table: the name of the table to query
    user: the username for DB access
    password: the password to be used for DB access
    join: an optional join clause  (omit the verb ‘join’)
    what: an optional clause indicating what to select

    Returns

    a list of column names"""
    ...
