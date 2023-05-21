"""
rdkit.ML.DecTree.CrossValidate module¶
handles doing cross validation with decision trees
This is, perhaps, a little misleading.  For the purposes of this module,
cross validation == evaluating the accuracy of a tree.
"""
from _typeshed import Incomplete
from rdkit.ML.Data import SplitData as SplitData
from rdkit.ML.DecTree import ID3
from rdkit.ML.DecTree import ID3 as ID3
from rdkit.ML.DecTree import randomtest as randomtest

def ChooseOptimalRoot(
    self,
    examples,
    trainExamples,
    testExamples,
    attrs,
    nPossibleVals,
    treeBuilder,
    nQuantBounds=[],
    **kwargs,
):
    """
    loops through all possible tree roots and chooses the one which produces the best tree
    Arguments

    examples: the full set of examples
    trainExamples: the training examples
    testExamples: the testing examples
    attrs: a list of attributes to consider in the tree building
    nPossibleVals: a list of the number of possible values each variable can adopt
    treeBuilder: the function to be used to actually build the tree
    nQuantBounds: an optional list.  If present, it’s assumed that the builder
    algorithm takes this argument as well (for building QuantTrees)

    Returns

    The best tree found

    Notes

    Trees are built using _trainExamples_
    Testing of each tree (to determine which is best) is done using _CrossValidate_ and
    the entire set of data (i.e. all of _examples_)
    _trainExamples_ is not used at all, which immediately raises the question of
    why it’s even being passed in"""
    ...

def CrossValidate(self, tree, testExamples, appendExamples=0):
    """
    Determines the classification error for the testExamples
    Arguments

    tree: a decision tree (or anything supporting a _ClassifyExample()_ method)
    testExamples: a list of examples to be used for testing
    appendExamples: a toggle which is passed along to the tree as it does
    the classification. The trees can use this to store the examples they
    classify locally.

    Returns

    a 2-tuple consisting of:

    the percent error of the tree
    a list of misclassified examples"""
    ...

def CrossValidationDriver(
    self,
    examples,
    attrs,
    nPossibleVals,
    holdOutFrac=0.3,
    silent=0,
    calcTotalError=0,
    treeBuilder=ID3.ID3Boot,
    lessGreedy=0,
    startAt=None,
    nQuantBounds=[],
    maxDepth=-1,
    **kwargs,
):
    """
    Driver function for building trees and doing cross validation
    Arguments

    examples: the full set of examples
    attrs: a list of attributes to consider in the tree building
    nPossibleVals: a list of the number of possible values each variable can adopt

    holdOutFrac: the fraction of the data which should be reserved for the hold-out set(used to calculate the error)

    silent: a toggle used to control how much visual noise this makes as it goes.
    calcTotalError: a toggle used to indicate whether the classification error
    of the tree should be calculated using the entire data set (when true) or just
    the training hold out set (when false)
    treeBuilder: the function to call to build the tree
    lessGreedy: toggles use of the less greedy tree growth algorithm (see
    _ChooseOptimalRoot_).
    startAt: forces the tree to be rooted at this descriptor
    nQuantBounds: an optional list.  If present, it’s assumed that the builder
    algorithm takes this argument as well (for building QuantTrees)
    maxDepth: an optional integer.  If present, it’s assumed that the builder
    algorithm takes this argument as well

    Returns

    a 2-tuple containing:

    the tree
    the cross-validation error of the tree"""
    ...

def TestRun(self):
    """
    testing code"""
    ...
