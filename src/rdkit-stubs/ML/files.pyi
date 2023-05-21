"""
rdkit.ML.files module¶
Generic file manipulation stuff
"""
from _typeshed import Incomplete

class ReFile(object):
    """
    convenience class for dealing with files with comments

    blank (all whitespace) lines, and lines beginning with commentcharacters are skipped.

    anything following a comment character on a line is stripped off"""

    def readline(self):
        """
        read the next line and return it.
        return ‘’ on EOF"""
        ...
    def readlines(self):
        """
        return a list of all the lines left in the file
        return [] if there are none"""
        ...
    def rewind(self):
        """
        rewinds the file (seeks to the beginning)"""
        ...
    regExp: Incomplete
    inFile: Incomplete

    def __init__(
        self, fileName, mode: str = ..., comment: str = ..., trailer: str = ...
    ) -> None: ...

def ReadDataFile(self, fileName, comment="#", depVarCol=0, dataType=float):
    """
    read in the data file and return a tuple of two Numeric arrays:
    (independent variables, dependant variables).
    ARGUMENTS:

    fileName: the fileName
    comment: the comment character for the file
    depVarcol: the column number containing the dependant variable
    dataType: the Numeric short-hand for the data type

    RETURNS:

    a tuple of two Numeric arrays:

    (independent variables, dependant variables)."""
    ...
