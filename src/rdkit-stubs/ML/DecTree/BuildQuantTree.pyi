"""
rdkit.ML.DecTree.BuildQuantTree module¶
"""
from _typeshed import Incomplete
from rdkit.ML.Data import Quantize as Quantize
from rdkit.ML.DecTree import ID3 as ID3
from rdkit.ML.DecTree import QuantTree as QuantTree
from rdkit.ML.InfoTheory import entropy as entropy

def BuildQuantTree(
    self,
    examples,
    target,
    attrs,
    nPossibleVals,
    nBoundsPerVar,
    depth=0,
    maxDepth=-1,
    exIndices=None,
    **kwargs,
):
    """
    Arguments

    examples: a list of lists (nInstances x nVariables+1) of variable
    values + instance values
    target: an int
    attrs: a list of ints indicating which variables can be used in the tree

    nPossibleVals: a list containing the number of possible values ofevery variable.

    nBoundsPerVar: the number of bounds to include for each variable
    depth: (optional) the current depth in the tree

    maxDepth: (optional) the maximum depth to which the treewill be grown

    Returns

    a QuantTree.QuantTreeNode with the decision tree

    NOTE: This code cannot bootstrap (start from nothing…)use _QuantTreeBoot_ (below) for that.
    """
    ...

def FindBest(
    self,
    resCodes,
    examples,
    nBoundsPerVar,
    nPossibleRes,
    nPossibleVals,
    attrs,
    exIndices=None,
    **kwargs,
): ...
def BuildQuantTree(
    self,
    examples,
    target,
    attrs,
    nPossibleVals,
    nBoundsPerVar,
    depth=0,
    maxDepth=-1,
    exIndices=None,
    **kwargs,
):
    """
    Arguments

    examples: a list of lists (nInstances x nVariables+1) of variable
    values + instance values
    target: an int
    attrs: a list of ints indicating which variables can be used in the tree

    nPossibleVals: a list containing the number of possible values ofevery variable.

    nBoundsPerVar: the number of bounds to include for each variable
    depth: (optional) the current depth in the tree

    maxDepth: (optional) the maximum depth to which the treewill be grown

    Returns

    a QuantTree.QuantTreeNode with the decision tree

    NOTE: This code cannot bootstrap (start from nothing…)use _QuantTreeBoot_ (below) for that.
    """
    ...

def QuantTreeBoot(
    self,
    examples,
    attrs,
    nPossibleVals,
    nBoundsPerVar,
    initialVar=None,
    maxDepth=-1,
    **kwargs,
):
    """
    Bootstrapping code for the QuantTree

    If _initialVar_ is not set, the algorithm will automaticallychoose the first variable in the tree (the standard greedy
    approach).  Otherwise, _initialVar_ will be used as the first
    split."""
    ...

def TestTree(self):
    """
    testing code for named trees"""
    ...

def TestQuantTree(self):
    """
    Testing code for named trees
    The created pkl file is required by the unit test code."""
    ...

def TestQuantTree2(self):
    """
    testing code for named trees
    The created pkl file is required by the unit test code."""
    ...

def TestTree(self):
    """
    testing code for named trees"""
    ...
