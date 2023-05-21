"""
rdkit.ML.Data.MLData module¶
classes to be used to help work with data sets
"""
from _typeshed import Incomplete

numericTypes: Incomplete

class MLDataSet(object):
    """
    A data set for holding general data (floats, ints, and strings)

    Notethis is intended to be a read-only data structure
    (i.e. after calling the constructor you cannot touch it)

    Constructor
    Arguments

    data: a list of lists containing the data. The data are copied, so don’t worryabout us overwriting them.

    nVars: the number of variables
    nPts: the number of points

    nPossibleVals: an list containing the number of possible valuesfor each variable (should contain 0 when not relevant)
    This is _nVars_ long

    qBounds: a list of lists containing quantization bounds for variableswhich are to be quantized (note, this class does not quantize
    the variables itself, it merely stores quantization bounds.
    an empty sublist indicates no quantization for a given variable
    This is _nVars_ long

    varNames: a list of the names of the variables.This is _nVars_ long

    ptNames: the names (labels) of the individual data pointsThis is _nPts_ long

    nResults: the number of results columns in the data lists.  This is usually1, but can be higher.
    """

    def AddPoint(self, pt): ...
    def AddPoints(self, pts, names): ...
    def GetAllData(self):
        """
        returns a copy of the data"""
        ...
    def GetInputData(self):
        """
        returns the input data
        Note

        _inputData_ means the examples without their result fields(the last _NResults_ entries)
        """
        ...
    def GetNPossibleVals(self): ...
    def GetNPts(self): ...
    data: Incomplete
    nResults: Incomplete
    nVars: Incomplete
    nPts: Incomplete
    qBounds: Incomplete
    nPossibleVals: Incomplete
    varNames: Incomplete
    ptNames: Incomplete

    def __init__(
        self,
        data,
        nVars: Incomplete | None = ...,
        nPts: Incomplete | None = ...,
        nPossibleVals: Incomplete | None = ...,
        qBounds: Incomplete | None = ...,
        varNames: Incomplete | None = ...,
        ptNames: Incomplete | None = ...,
        nResults: int = ...,
    ) -> None: ...
    def GetNResults(self): ...
    def GetNVars(self): ...
    def GetNamedData(self):
        """
        returns a list of named examples
        Note

        a named example is the result of prepending the examplename to the data list"""
        ...
    def GetPtNames(self): ...
    def GetNPts(self): ...
    def GetNPossibleVals(self): ...
    def GetQuantBounds(self): ...
    def __getitem__(self, idx): ...
    def __setitem__(self, idx, val): ...
    def GetNamedData(self):
        """
        returns a list of named examples
        Note

        a named example is the result of prepending the examplename to the data list"""
        ...
    def GetAllData(self):
        """
        returns a copy of the data"""
        ...
    def GetInputData(self):
        """
        returns the input data
        Note

        _inputData_ means the examples without their result fields(the last _NResults_ entries)
        """
        ...
    def GetResults(self):
        """
        Returns the result fields from each example"""
        ...
    def GetVarNames(self): ...
    def GetPtNames(self): ...
    def AddPoint(self, pt): ...
    def AddPoints(self, pts, names): ...

class MLQuantDataSet(MLDataSet):
    """
    a data set for holding quantized data
    Note

    this is intended to be a read-only data structure
    (i.e. after calling the constructor you cannot touch it)

    Big differences to MLDataSet

    data are stored in a numpy array since they are homogenous
    results are assumed to be quantized (i.e. no qBounds entry is required)

    Constructor
    Arguments

    data: a list of lists containing the data. The data are copied, so don’t worryabout us overwriting them.

    nVars: the number of variables
    nPts: the number of points

    nPossibleVals: an list containing the number of possible valuesfor each variable (should contain 0 when not relevant)
    This is _nVars_ long

    qBounds: a list of lists containing quantization bounds for variableswhich are to be quantized (note, this class does not quantize
    the variables itself, it merely stores quantization bounds.
    an empty sublist indicates no quantization for a given variable
    This is _nVars_ long

    varNames: a list of the names of the variables.This is _nVars_ long

    ptNames: the names (labels) of the individual data pointsThis is _nPts_ long

    nResults: the number of results columns in the data lists.  This is usually1, but can be higher.
    """

    def GetNamedData(self):
        """
        returns a list of named examples
        Note

        a named example is the result of prepending the examplename to the data list"""
        ...
    def GetAllData(self):
        """
        returns a copy of the data"""
        ...
    def GetInputData(self):
        """
        returns the input data
        Note

        _inputData_ means the examples without their result fields(the last _NResults_ entries)
        """
        ...
    def GetNamedData(self):
        """
        returns a list of named examples
        Note

        a named example is the result of prepending the examplename to the data list"""
        ...
    def GetResults(self):
        """
        Returns the result fields from each example"""
        ...
    data: Incomplete
    nResults: Incomplete
    nVars: Incomplete
    nPts: Incomplete
    qBounds: Incomplete
    nPossibleVals: Incomplete
    varNames: Incomplete
    ptNames: Incomplete

    def __init__(
        self,
        data,
        nVars: Incomplete | None = ...,
        nPts: Incomplete | None = ...,
        nPossibleVals: Incomplete | None = ...,
        qBounds: Incomplete | None = ...,
        varNames: Incomplete | None = ...,
        ptNames: Incomplete | None = ...,
        nResults: int = ...,
    ) -> None: ...
