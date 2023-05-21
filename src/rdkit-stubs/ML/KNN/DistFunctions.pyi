"""
rdkit.ML.KNN.DistFunctions module¶
"""

def EuclideanDist(self, ex1, ex2, attrs):
    """
    >>> v1 = [0,1,0,1]
    >>> v2 = [1,0,1,0]
    >>> EuclideanDist(v1,v2,range(4))
    2.0
    >>> EuclideanDist(v1,v1,range(4))
    0.0
    >>> v2 = [0,0,0,1]
    >>> EuclideanDist(v1,v2,range(4))
    1.0
    >>> v2 = [0,.5,0,.5]
    >>> abs(EuclideanDist(v1,v2,range(4))-1./math.sqrt(2))<1e-4
    1"""
    ...

def TanimotoDist(self, ex1, ex2, attrs):
    """
    >>> v1 = [0,1,0,1]
    >>> v2 = [1,0,1,0]
    >>> TanimotoDist(v1,v2,range(4))
    1.0
    >>> v2 = [1,0,1,1]
    >>> TanimotoDist(v1,v2,range(4))
    0.75
    >>> TanimotoDist(v2,v2,range(4))
    0.0

    # this tests Issue 122
    >>> v3 = [0,0,0,0]
    >>> TanimotoDist(v3,v3,range(4))
    1.0"""
    ...
