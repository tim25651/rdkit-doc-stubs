"""
rdkit.ML.CompositeRun module¶
contains a class to store parameters for and results from
Composite building
"""
from _typeshed import Incomplete
from rdkit import RDConfig as RDConfig
from rdkit.Dbase import DbModule as DbModule
from rdkit.Dbase.DbConnection import DbConnect as DbConnect

class CompositeRun(object):
    """
    class to store parameters for and results from Composite building
    This class has a default set of fields which are added to the database.

    By default these fields are stored in a tuple, so they are immutable.  Thisis probably what you want.
    """

    fields: tuple[tuple[str]] = ...

    def Store(
        self, db="models.gdb", table="results", user="sysdba", password="masterkey"
    ):
        """
        adds the result to a database
        Arguments

        db: name of the database to use
        table: name of the table to use
        user&password: connection information"""
        ...
    def GetDataSet(self, **kwargs):
        """
        Returns a MLDataSet pulled from a database using our stored
        values."""
        ...
    def GetDataSetInfo(self, **kwargs):
        """
        Returns a MLDataSet pulled from a database using our stored
        values."""
        ...
    def Store(
        self, db="models.gdb", table="results", user="sysdba", password="masterkey"
    ):
        """
        adds the result to a database
        Arguments

        db: name of the database to use
        table: name of the table to use
        user&password: connection information"""
        ...

def SetDefaults(self, runDetails):
    """
    initializes a details object with default values
    Arguments

    details:  (optional) a _CompositeRun.CompositeRun_ object.
    If this is not provided, the global _runDetails will be used.

    Returns

    the initialized _CompositeRun_ object."""
    ...

class CompositeRun(object):
    """
    class to store parameters for and results from Composite building
    This class has a default set of fields which are added to the database.

    By default these fields are stored in a tuple, so they are immutable.  Thisis probably what you want.
    """

    fields: tuple[tuple[str]] = ...

    def Store(
        self, db="models.gdb", table="results", user="sysdba", password="masterkey"
    ):
        """
        adds the result to a database
        Arguments

        db: name of the database to use
        table: name of the table to use
        user&password: connection information"""
        ...
    def GetDataSet(self, **kwargs):
        """
        Returns a MLDataSet pulled from a database using our stored
        values."""
        ...
    def GetDataSetInfo(self, **kwargs):
        """
        Returns a MLDataSet pulled from a database using our stored
        values."""
        ...
    def Store(
        self, db="models.gdb", table="results", user="sysdba", password="masterkey"
    ):
        """
        adds the result to a database
        Arguments

        db: name of the database to use
        table: name of the table to use
        user&password: connection information"""
        ...
