"""
rdkit.ML.Data.Stats module¶
various statistical operations on data
"""
from _typeshed import Incomplete

def FormCorrelationMatrix(self, mat):
    """
    form and return the covariance matrix"""
    ...

def StandardizeMatrix(self, mat):
    """
    This is the standard subtract off the average and divide by the deviation
    standardization function.

    Arguments

    mat: a numpy array

    Notes

    in addition to being returned, _mat_ is modified in place, so beware"""
    ...

def FormCovarianceMatrix(self, mat):
    """
    form and return the covariance matrix"""
    ...

def GetConfidenceInterval(self, sd, n, level=95): ...
def MeanAndDev(self, vect, sampleSD=1):
    """
    returns the mean and standard deviation of a vector"""
    ...

def FormCorrelationMatrix(self, mat):
    """
    form and return the covariance matrix"""
    ...

def PrincipalComponents(self, mat, reverseOrder=1):
    """
    do a principal components analysis"""
    ...

def R2(self, orig, residSum):
    """
    returns the R2 value for a set of predictions"""
    ...

def StandardizeMatrix(self, mat):
    """
    This is the standard subtract off the average and divide by the deviation
    standardization function.

    Arguments

    mat: a numpy array

    Notes

    in addition to being returned, _mat_ is modified in place, so beware"""
    ...

def TransformPoints(self, tFormMat, pts):
    """
    transforms a set of points using tFormMat
    Arguments

    tFormMat: a numpy array
    pts: a list of numpy arrays (or a 2D array)

    Returns

    a list of numpy arrays"""
    ...

def MeanAndDev(self, vect, sampleSD=1):
    """
    returns the mean and standard deviation of a vector"""
    ...

def R2(self, orig, residSum):
    """
    returns the R2 value for a set of predictions"""
    ...

tConfs: Incomplete
tTable: Incomplete

def GetConfidenceInterval(self, sd, n, level=95): ...
