"""
rdkit.ML.KNN.KNNClassificationModel module¶
Define the class _KNNClassificationModel_, used to represent a k-nearest neighbhors
classification model
Inherits from _KNNModel_
"""
from _typeshed import Incomplete
from rdkit.ML.KNN import KNNModel
from rdkit.ML.KNN import KNNModel as KNNModel

class KNNClassificationModel(KNNModel.KNNModel):
    """
    This is used to represent a k-nearest neighbor classifier"""

    def ClassifyExample(self, example, appendExamples=0, neighborList=None):
        """
        Classify a an example by looking at its closest neighbors
        The class assigned to this example is same as the class for the mojority of its
        _k neighbors
        Arguments

        examples: the example to be classified
        appendExamples: if this is nonzero then the example will be stored on this model
        neighborList: if provided, will be used to return the list of neighbors

        Returns

        the classification of _example_"""
        ...
    def __init__(self, k, attrs, dfunc, radius: Incomplete | None = ...) -> None: ...
    def type(self): ...
    def SetBadExamples(self, examples): ...
    def GetBadExamples(self): ...
    def NameModel(self, varNames): ...
    def SetBadExamples(self, examples): ...
    def type(self): ...
    def ClassifyExample(self, example, appendExamples=0, neighborList=None):
        """
        Classify a an example by looking at its closest neighbors
        The class assigned to this example is same as the class for the mojority of its
        _k neighbors
        Arguments

        examples: the example to be classified
        appendExamples: if this is nonzero then the example will be stored on this model
        neighborList: if provided, will be used to return the list of neighbors

        Returns

        the classification of _example_"""
        ...
