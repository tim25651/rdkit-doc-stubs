"""
rdkit.ML.DecTree.PruneTree module¶
Contains functionality for doing tree pruning
"""
from rdkit.ML.DecTree import CrossValidate as CrossValidate
from rdkit.ML.DecTree import DecTree as DecTree

def MaxCount(self, examples):
    """
    given a set of examples, returns the most common result code
    Arguments

    examples: a list of examples to be counted

    Returns

    the most common result code"""
    ...

def PruneTree(self, tree, trainExamples, testExamples, minimizeTestErrorOnly=1):
    """
    implements a reduced-error pruning of decision trees
    This algorithm is described on page 69 of Mitchell’s book.
    Pruning can be done using just the set of testExamples (the validation set)
    or both the testExamples and the trainExamples by setting minimizeTestErrorOnly
    to 0.
    Arguments

    tree: the initial tree to be pruned
    trainExamples: the examples used to train the tree
    testExamples: the examples held out for testing the tree
    minimizeTestErrorOnly: if this toggle is zero, all examples (i.e.
    _trainExamples_ + _testExamples_ will be used to evaluate the error.

    Returns

    a 2-tuple containing:

    the best tree
    the best error (the one which corresponds to that tree)"""
    ...
