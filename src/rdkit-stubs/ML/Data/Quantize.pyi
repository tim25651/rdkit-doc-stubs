"""
rdkit.ML.Data.Quantize module¶
Automatic search for quantization bounds
This uses the expected informational gain to determine where quantization bounds should
lie.
Notes:

bounds are less than, so if the bounds are [1.,2.],
[0.9,1.,1.1,2.,2.2] -> [0,1,1,2,2]
"""
from rdkit.ML.Data import cQuantize as cQuantize
from rdkit.ML.InfoTheory import entropy as entropy

def FindVarMultQuantBounds(self, vals, nBounds, results, nPossibleRes):
    """
    finds multiple quantization bounds for a single variable
    Arguments

    vals: sequence of variable values (assumed to be floats)
    nBounds: the number of quantization bounds to find
    results: a list of result codes (should be integers)
    nPossibleRes: an integer with the number of possible values of the
    result variable

    Returns

    a 2-tuple containing:

    a list of the quantization bounds (floats)
    the information gain associated with this quantization"""
    ...

def FindVarQuantBound(self, vals, results, nPossibleRes):
    """
    Uses FindVarMultQuantBounds, only here for historic reasons"""
    ...

hascQuantize: int

def feq(self, v1, v2, tol=1e-08):
    """
    floating point equality with a tolerance factor
    Arguments

    v1: a float
    v2: a float
    tol: the tolerance for comparison

    Returns

    0 or 1"""
    ...

def FindVarQuantBound(self, vals, results, nPossibleRes):
    """
    Uses FindVarMultQuantBounds, only here for historic reasons"""
    ...

def FindVarMultQuantBounds(self, vals, nBounds, results, nPossibleRes):
    """
    finds multiple quantization bounds for a single variable
    Arguments

    vals: sequence of variable values (assumed to be floats)
    nBounds: the number of quantization bounds to find
    results: a list of result codes (should be integers)
    nPossibleRes: an integer with the number of possible values of the
    result variable

    Returns

    a 2-tuple containing:

    a list of the quantization bounds (floats)
    the information gain associated with this quantization"""
    ...
