"""
rdkit.Dbase.DbUtils module¶
a set of functions for interacting with databases
When possible, it’s probably preferable to use a _DbConnection.DbConnect_ object
"""
from _typeshed import Incomplete
from rdkit.Dbase import DbInfo as DbInfo
from rdkit.Dbase import DbModule as DbModule
from rdkit.Dbase.DbResultSet import DbResultSet as DbResultSet
from rdkit.Dbase.DbResultSet import RandomAccessDbResultSet as RandomAccessDbResultSet

def DatabaseToDatabase(
    self,
    fromDb,
    fromTbl,
    toDb,
    toTbl,
    fields="*",
    join="",
    where="",
    user="sysdba",
    password="masterkey",
    keyCol=None,
    nullMarker="None",
):
    """
    FIX: at the moment this is a hack"""
    ...

def DatabaseToText(
    self,
    dBase,
    table,
    fields="*",
    join="",
    where="",
    user="sysdba",
    password="masterkey",
    delim=",",
    cn=None,
):
    """
    Pulls the contents of a database and makes a deliminted text file from them

    Arguments
    dBase: the name of the DB file to be used
    table: the name of the table to query
    fields: the fields to select with the SQL query
    join: the join clause of the SQL query
    (e.g. ‘join foo on foo.bar=base.bar’)
    where: the where clause of the SQL query
    (e.g. ‘where foo = 2’ or ‘where bar > 17.6’)
    user: the username for DB access
    password: the password to be used for DB access

    Returns

    the CSV data (as text)"""
    ...

def GetColumns(
    self,
    dBase,
    table,
    fieldString,
    user="sysdba",
    password="masterkey",
    join="",
    cn=None,
):
    """
    gets a set of data from a table
    Arguments

    dBase: database name
    table: table name

    fieldString: a string with the names of the fields to be extracted,this should be a comma delimited list

    user and  password:
    join: a join clause (omit the verb ‘join’)

    Returns

    a list of the data"""
    ...

def GetData(
    self,
    dBase,
    table,
    fieldString="*",
    whereString="",
    user="sysdba",
    password="masterkey",
    removeDups=-1,
    join="",
    forceList=0,
    transform=None,
    randomAccess=1,
    extras=None,
    cn=None,
):
    """
    a more flexible method to get a set of data from a table
    Arguments

    fields: a string with the names of the fields to be extracted,this should be a comma delimited list

    where: the SQL where clause to be used with the DB query

    removeDups indicates the column which should be used to screenout duplicates.  Only the first appearance of a duplicate will
    be left in the dataset.

    Returns

    a list of the data

    Notes

    EFF: this isn’t particularly efficient"""
    ...

def DatabaseToText(
    self,
    dBase,
    table,
    fields="*",
    join="",
    where="",
    user="sysdba",
    password="masterkey",
    delim=",",
    cn=None,
):
    """
    Pulls the contents of a database and makes a deliminted text file from them

    Arguments
    dBase: the name of the DB file to be used
    table: the name of the table to query
    fields: the fields to select with the SQL query
    join: the join clause of the SQL query
    (e.g. ‘join foo on foo.bar=base.bar’)
    where: the where clause of the SQL query
    (e.g. ‘where foo = 2’ or ‘where bar > 17.6’)
    user: the username for DB access
    password: the password to be used for DB access

    Returns

    the CSV data (as text)"""
    ...

def TypeFinder(self, data, nRows, nCols, nullMarker=None):
    """
    finds the types of the columns in _data_

    if nullMarker is not None, elements of the data table which areequal to nullMarker will not count towards setting the type of
    their columns."""
    ...

def GetTypeStrings(self, colHeadings, colTypes, keyCol=None):
    """
    returns a list of SQL type strings"""
    ...

def TextFileToDatabase(
    self,
    dBase,
    table,
    inF,
    delim=",",
    user="sysdba",
    password="masterkey",
    maxColLabelLen=31,
    keyCol=None,
    nullMarker=None,
):
    """
    loads the contents of the text file into a database.
    Arguments

    dBase: the name of the DB to use
    table: the name of the table to create/overwrite
    inF: the file like object from which the data should
    be pulled (must support readline())
    delim: the delimiter used to separate fields
    user: the user name to use in connecting to the DB
    password: the password to use in connecting to the DB
    maxColLabelLen: the maximum length a column label should be
    allowed to have (truncation otherwise)
    keyCol: the column to be used as an index for the db

    Notes

    if _table_ already exists, it is destroyed before we write
    the new data
    we assume that the first row of the file contains the column names"""
    ...

def TypeFinder(self, data, nRows, nCols, nullMarker=None):
    """
    finds the types of the columns in _data_

    if nullMarker is not None, elements of the data table which areequal to nullMarker will not count towards setting the type of
    their columns."""
    ...

def DatabaseToDatabase(
    self,
    fromDb,
    fromTbl,
    toDb,
    toTbl,
    fields="*",
    join="",
    where="",
    user="sysdba",
    password="masterkey",
    keyCol=None,
    nullMarker="None",
):
    """
    FIX: at the moment this is a hack"""
    ...
