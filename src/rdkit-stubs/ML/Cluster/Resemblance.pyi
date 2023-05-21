"""
rdkit.ML.Cluster.Resemblance module¶
code for dealing with resemblance (metric) matrices
Here’s how the matrices are stored:

‘[(0,1),(0,2),(1,2),(0,3),(1,3),(2,3)…]  (row,col), col>row’
or, alternatively the matrix can be drawn, with indices as:

|| - || 0 || 1 || 3
|| - || - || 2 || 4
|| - || - || - || 5
|| - || - || - || -

the index of a given (row,col) pair is:‘(col*(col-1))/2 + row’
"""
from _typeshed import Incomplete

def CalcMetricMatrix(self, inData, metricFunc):
    """
    generates a metric matrix

    Arguments
    inData is assumed to be a list of clusters (or anything with
    a GetPosition() method)
    metricFunc is the function to be used to generate the matrix

    Returns

    the metric matrix as a Numeric array"""
    ...

def EuclideanDistance(self, inData):
    """
    returns the euclidean metricMat between the points in _inData_
    Arguments

    inData: a Numeric array of data points

    Returns

    a Numeric array with the metric matrix.  See the module documentation
    for the format."""
    ...

def CalcMetricMatrix(self, inData, metricFunc):
    """
    generates a metric matrix

    Arguments
    inData is assumed to be a list of clusters (or anything with
    a GetPosition() method)
    metricFunc is the function to be used to generate the matrix

    Returns

    the metric matrix as a Numeric array"""
    ...

def FindMinValInList(self, mat, nObjs, minIdx=None):
    """
    finds the minimum value in a metricMatrix and returns it and its indices
    Arguments

    mat: the metric matrix
    nObjs: the number of objects to be considered
    minIdx: the index of the minimum value (value, row and column still need
    to be calculated

    Returns

    a 3-tuple containing:

    the row
    the column
    the minimum value itself

    Notes

    -this probably ain’t the speediest thing on earth"""
    ...

def ShowMetricMat(self, metricMat, nObjs):
    """
    displays a metric matrix
    Arguments

    metricMat: the matrix to be displayed
    nObjs: the number of objects to display"""
    ...

methods: Incomplete