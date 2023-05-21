"""
rdkit.ML.Data.DataUtils module¶
Utilities for data manipulation
FILE FORMATS:

.qdat files contain quantized data suitable for

feeding to learning algorithms.
The .qdat file, written by _DecTreeGui_, is structured as follows:

Any number of lines which are ignored.
A line containing the string ‘Variable Table’
any number of variable definitions in the format:
‘# Variable_name [quant_bounds]’

where ‘[quant_bounds]’ is a list of the boundaries used for quantizingthat variable.  If the variable is inherently integral (i.e. not
quantized), this can be an empty list.

A line beginning with ‘# —-’ which signals the end of the variable list
Any number of lines containing data points, in the format:
‘Name_of_point var1 var2 var3 …. varN’
all variable values should be integers

Throughout, it is assumed that varN is the result

.dat files contain the same information as .qdat files, but the variable
values can be anything (floats, ints, strings).  These files should
still contain quant_bounds!
.qdat.pkl file contain a pickled (binary) representation of
the data read in.  They stores, in order:

A python list of the variable names
A python list of lists with the quantization bounds
A python list of the point names
A python list of lists with the data points
"""
from _typeshed import Incomplete
from rdkit.DataStructs import BitUtils as BitUtils
from rdkit.ML.Data import MLData as MLData
from rdkit.utils import fileutils as fileutils

def BuildDataSet(self, fileName):
    """
    builds a data set from a .dat file
    Arguments

    fileName: the name of the .dat file

    Returns

    an _MLData.MLDataSet_"""
    ...

def permutation(self, nToDo): ...
def WriteData(self, outFile, varNames, qBounds, examples):
    """
    writes out a .qdat file
    Arguments

    outFile: a file object
    varNames: a list of variable names

    qBounds: the list of quantization bounds (should be the same lengthas _varNames_)

    examples: the data to be written"""
    ...

def ReadVars(self, inFile):
    """
    reads the variables and quantization bounds from a .qdat or .dat file
    Arguments

    inFile: a file object

    Returns

    a 2-tuple containing:

    varNames: a list of the variable names
    qbounds: the list of quantization bounds for each variable"""
    ...

def ReadQuantExamples(self, inFile):
    """
    reads the examples from a .qdat file
    Arguments

    inFile: a file object

    Returns

    a 2-tuple containing:

    the names of the examples
    a list of lists containing the examples themselves

    Note

    because this is reading a .qdat file, it assumed that all variable values
    are integers"""
    ...

def ReadGeneralExamples(self, inFile):
    """
    reads the examples from a .dat file
    Arguments

    inFile: a file object

    Returns

    a 2-tuple containing:

    the names of the examples
    a list of lists containing the examples themselves

    Note

    this attempts to convert variable values to ints, then floats.
    if those both fail, they are left as strings"""
    ...

def BuildQuantDataSet(self, fileName):
    """
    builds a data set from a .qdat file
    Arguments

    fileName: the name of the .qdat file

    Returns

    an _MLData.MLQuantDataSet_"""
    ...

def BuildDataSet(self, fileName):
    """
    builds a data set from a .dat file
    Arguments

    fileName: the name of the .dat file

    Returns

    an _MLData.MLDataSet_"""
    ...

def CalcNPossibleUsingMap(self, data, order, qBounds, nQBounds=None, silent=True):
    """
    calculates the number of possible values for each variable in a data set
    Arguments

    data: a list of examples
    order: the ordering map between the variables in _data_ and _qBounds_
    qBounds: the quantization bounds for the variables

    Returns

    a list with the number of possible values each variable takes on in the data set

    Notes

    variables present in _qBounds_ will have their _nPossible_ number read
    from _qbounds
    _nPossible_ for other numeric variables will be calculated"""
    ...

def CountResults(self, inData, col=-1, bounds=None):
    """
    #DOC"""
    ...

def WritePickledData(self, outName, data):
    """
    writes either a .qdat.pkl or a .dat.pkl file
    Arguments

    outName: the name of the file to be used
    data: either an _MLData.MLDataSet_ or an _MLData.MLQuantDataSet_"""
    ...

def TakeEnsemble(self, vect, ensembleIds, isDataVect=False):
    """
    >>> v = [10,20,30,40,50]
    >>> TakeEnsemble(v,(1,2,3))
    [20, 30, 40]
    >>> v = ['foo',10,20,30,40,50,1]
    >>> TakeEnsemble(v,(1,2,3),isDataVect=True)
    ['foo', 20, 30, 40, 1]"""
    ...

def DBToData(
    self,
    dbName,
    tableName,
    user="sysdba",
    password="masterkey",
    dupCol=-1,
    what="*",
    where="",
    join="",
    pickleCol=-1,
    pickleClass=None,
    ensembleIds=None,
):
    """
    constructs  an _MLData.MLDataSet_ from a database
    Arguments

    dbName: the name of the database to be opened
    tableName: the table name containing the data in the database
    user: the user name to be used to connect to the database
    password: the password to be used to connect to the database
    dupCol: if nonzero specifies which column should be used to recognize
    duplicates.

    Returns

    an _MLData.MLDataSet_

    Notes

    this uses Dbase.DataUtils functionality"""
    ...

def FilterData(self, inData, val, frac, col=-1, indicesToUse=None, indicesOnly=0):
    """
    #DOC"""
    ...

def TextToData(self, reader, ignoreCols=[], onlyCols=None):
    """
    constructs  an _MLData.MLDataSet_ from a bunch of text
    #DOC

    Arguments
    reader needs to be iterable and return lists of elements
    (like a csv.reader)

    Returns

    an _MLData.MLDataSet_"""
    ...

def TextFileToData(self, fName, onlyCols=None):
    """
    #DOC"""
    ...

def InitRandomNumbers(self, seed):
    """
    Seeds the random number generators
    Arguments

    seed: a 2-tuple containing integers to be used as the random number seeds

    Notes

    this seeds both the RDRandom generator and the one in the standard
    Python _random_ module"""
    ...

def FilterData(self, inData, val, frac, col=-1, indicesToUse=None, indicesOnly=0):
    """
    #DOC"""
    ...

def CountResults(self, inData, col=-1, bounds=None):
    """
    #DOC"""
    ...

def RandomizeActivities(self, dataSet, shuffle=0, runDetails=None):
    """
    randomizes the activity values of a dataset
    Arguments

    dataSet: a _ML.Data.MLQuantDataSet_, the activities here will be randomized
    shuffle: an optional toggle. If this is set, the activity values
    will be shuffled (so the number in each class remains constant)
    runDetails: an optional CompositeRun object

    Note

    _examples_ are randomized in place"""
    ...

def ReadGeneralExamples(self, inFile):
    """
    reads the examples from a .dat file
    Arguments

    inFile: a file object

    Returns

    a 2-tuple containing:

    the names of the examples
    a list of lists containing the examples themselves

    Note

    this attempts to convert variable values to ints, then floats.
    if those both fail, they are left as strings"""
    ...

def ReadQuantExamples(self, inFile):
    """
    reads the examples from a .qdat file
    Arguments

    inFile: a file object

    Returns

    a 2-tuple containing:

    the names of the examples
    a list of lists containing the examples themselves

    Note

    because this is reading a .qdat file, it assumed that all variable values
    are integers"""
    ...

def ReadVars(self, inFile):
    """
    reads the variables and quantization bounds from a .qdat or .dat file
    Arguments

    inFile: a file object

    Returns

    a 2-tuple containing:

    varNames: a list of the variable names
    qbounds: the list of quantization bounds for each variable"""
    ...

def TakeEnsemble(self, vect, ensembleIds, isDataVect=False):
    """
    >>> v = [10,20,30,40,50]
    >>> TakeEnsemble(v,(1,2,3))
    [20, 30, 40]
    >>> v = ['foo',10,20,30,40,50,1]
    >>> TakeEnsemble(v,(1,2,3),isDataVect=True)
    ['foo', 20, 30, 40, 1]"""
    ...

def TextFileToData(self, fName, onlyCols=None):
    """
    #DOC"""
    ...

def TextToData(self, reader, ignoreCols=[], onlyCols=None):
    """
    constructs  an _MLData.MLDataSet_ from a bunch of text
    #DOC

    Arguments
    reader needs to be iterable and return lists of elements
    (like a csv.reader)

    Returns

    an _MLData.MLDataSet_"""
    ...

def WriteData(self, outFile, varNames, qBounds, examples):
    """
    writes out a .qdat file
    Arguments

    outFile: a file object
    varNames: a list of variable names

    qBounds: the list of quantization bounds (should be the same lengthas _varNames_)

    examples: the data to be written"""
    ...

def WritePickledData(self, outName, data):
    """
    writes either a .qdat.pkl or a .dat.pkl file
    Arguments

    outName: the name of the file to be used
    data: either an _MLData.MLDataSet_ or an _MLData.MLQuantDataSet_"""
    ...

def permutation(self, nToDo): ...
