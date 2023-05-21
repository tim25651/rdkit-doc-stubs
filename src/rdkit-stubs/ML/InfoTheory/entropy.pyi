"""
rdkit.ML.InfoTheory.entropy module¶
Informational Entropy functions
The definitions used are the same as those in Tom Mitchell’s
book “Machine Learning”
"""
from _typeshed import Incomplete

hascEntropy: int

def PyInfoEntropy(self, results):
    """
    Calculates the informational entropy of a set of results.
    Arguments

    results is a 1D Numeric array containing the number of times a
    given set hits each possible result.
    For example, if a function has 3 possible results, and the

    variable in question hits them 5, 6 and 1 times each,
    results would be [5,6,1]

    Returns

    the informational entropy"""
    ...

def PyInfoGain(self, varMat):
    """
    calculates the information gain for a variable
    Arguments

    varMat is a Numeric array with the number of possible occurrencesof each result for reach possible value of the given variable.

    So, for a variable which adopts 4 possible values and a result whichhas 3 possible values, varMat would be 4x3

    Returns

    The expected information gain"""
    ...

InfoEntropy = PyInfoEntropy
InfoGain = PyInfoGain
InfoEntropy = PyInfoEntropy
InfoGain = PyInfoGain
