"""
rdkit.ML.KNN.KNNRegressionModel module¶
Define the class _KNNRegressionModel_, used to represent a k-nearest neighbhors
regression model

Inherits from _KNNModel_
"""
from _typeshed import Incomplete
from rdkit.ML.KNN import KNNModel
from rdkit.ML.KNN import KNNModel as KNNModel

class KNNRegressionModel(KNNModel.KNNModel):
    """
    This is used to represent a k-nearest neighbor classifier"""

    def __init__(self, k, attrs, dfunc, radius: Incomplete | None = ...) -> None: ...
    def type(self): ...
    def SetBadExamples(self, examples): ...
    def GetBadExamples(self): ...
    def NameModel(self, varNames): ...
    def PredictExample(
        self, example, appendExamples=0, weightedAverage=0, neighborList=None
    ):
        """
        Generates a prediction for an example by looking at its closest neighbors
        Arguments

        examples: the example to be classified
        appendExamples: if this is nonzero then the example will be stored on this model

        weightedAverage: if provided, the neighbors’ contributions to the value will beweighed by their reciprocal square distance

        neighborList: if provided, will be used to return the list of neighbors

        Returns

        the classification of _example_"""
        ...
    def SetBadExamples(self, examples): ...
    def type(self): ...
