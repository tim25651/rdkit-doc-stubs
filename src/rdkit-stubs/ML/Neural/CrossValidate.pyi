"""
rdkit.ML.Neural.CrossValidate module¶
handles doing cross validation with neural nets
This is, perhaps, a little misleading.  For the purposes of this module,
cross validation == evaluating the accuracy of a net.
"""
from _typeshed import Incomplete
from rdkit.ML.Data import SplitData as SplitData
from rdkit.ML.Neural import Network as Network
from rdkit.ML.Neural import Trainers as Trainers

def CrossValidate(self, net, testExamples, tolerance, appendExamples=0):
    """
    Determines the classification error for the testExamplesArguments

    tree: a decision tree (or anything supporting a _ClassifyExample()_ method)
    testExamples: a list of examples to be used for testing

    appendExamples: a toggle which is ignored, it’s just here to maintainthe same API as the decision tree code.

    Returns

    a 2-tuple consisting of:

    the percent error of the net
    a list of misclassified examples

    NoteAt the moment, this is specific to nets with only one output"""
    ...

def CrossValidationDriver(
    self,
    examples,
    attrs=[],
    nPossibleVals=[],
    holdOutFrac=0.3,
    silent=0,
    tolerance=0.3,
    calcTotalError=0,
    hiddenSizes=None,
    **kwargs,
):
    """
    Arguments

    examples: the full set of examples

    attrs: a list of attributes to consider in the tree buildingThis argument is ignored

    nPossibleVals: a list of the number of possible values each variable can adoptThis argument is ignored

    holdOutFrac: the fraction of the data which should be reserved for the hold-out set(used to calculate the error)

    silent: a toggle used to control how much visual noise this makes as it goes.
    tolerance: the tolerance for convergence of the net

    calcTotalError: if this is true the entire data set is used to calculateaccuracy of the net

    hiddenSizes: a list containing the size(s) of the hidden layers in the network.if _hiddenSizes_ is None, one hidden layer containing the same number of nodes
    as the input layer will be used

    Returns

    a 2-tuple containing:

    the net
    the cross-validation error of the net

    NoteAt the moment, this is specific to nets with only one output"""
    ...
