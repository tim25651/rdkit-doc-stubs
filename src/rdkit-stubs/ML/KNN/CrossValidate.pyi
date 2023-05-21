"""
rdkit.ML.KNN.CrossValidate module¶
handles doing cross validation with k-nearest neighbors model
and evaluation of individual models
"""
from rdkit.ML.Data import SplitData as SplitData
from rdkit.ML.KNN import DistFunctions
from rdkit.ML.KNN import DistFunctions as DistFunctions
from rdkit.ML.KNN.KNNClassificationModel import (
    KNNClassificationModel as KNNClassificationModel,
)
from rdkit.ML.KNN.KNNRegressionModel import KNNRegressionModel as KNNRegressionModel

def CrossValidate(self, knnMod, testExamples, appendExamples=0):
    """
    Determines the classification error for the testExamples
    Arguments

    tree: a decision tree (or anything supporting a _ClassifyExample()_ method)
    testExamples: a list of examples to be used for testing
    appendExamples: a toggle which is passed along to the tree as it does
    the classification. The trees can use this to store the examples they
    classify locally.

    Returns

    a 2-tuple consisting of:"""
    ...

def CrossValidationDriver(
    self,
    examples,
    attrs,
    nPossibleValues,
    numNeigh,
    modelBuilder=makeClassificationModel,
    distFunc=DistFunctions.EuclideanDist,
    holdOutFrac=0.3,
    silent=0,
    calcTotalError=0,
    **kwargs,
):
    """
    Driver function for building a KNN model of a specified type
    Arguments

    examples: the full set of examples
    numNeigh: number of neighbors for the KNN model (basically k in k-NN)
    knnModel: the type of KNN model (a classification vs regression model)
    holdOutFrac: the fraction of the data which should be reserved for the hold-out set
    (used to calculate error)
    silent: a toggle used to control how much visual noise this makes as it goes
    calcTotalError: a toggle used to indicate whether the classification error
    of the tree should be calculated using the entire data set (when true) or just
    the training hold out set (when false)"""
    ...

def makeClassificationModel(self, numNeigh, attrs, distFunc): ...
def makeRegressionModel(self, numNeigh, attrs, distFunc): ...
def CrossValidate(self, knnMod, testExamples, appendExamples=0):
    """
    Determines the classification error for the testExamples
    Arguments

    tree: a decision tree (or anything supporting a _ClassifyExample()_ method)
    testExamples: a list of examples to be used for testing
    appendExamples: a toggle which is passed along to the tree as it does
    the classification. The trees can use this to store the examples they
    classify locally.

    Returns

    a 2-tuple consisting of:"""
    ...

def CrossValidationDriver(
    self,
    examples,
    attrs,
    nPossibleValues,
    numNeigh,
    modelBuilder=makeClassificationModel,
    distFunc=DistFunctions.EuclideanDist,
    holdOutFrac=0.3,
    silent=0,
    calcTotalError=0,
    **kwargs,
):
    """
    Driver function for building a KNN model of a specified type
    Arguments

    examples: the full set of examples
    numNeigh: number of neighbors for the KNN model (basically k in k-NN)
    knnModel: the type of KNN model (a classification vs regression model)
    holdOutFrac: the fraction of the data which should be reserved for the hold-out set
    (used to calculate error)
    silent: a toggle used to control how much visual noise this makes as it goes
    calcTotalError: a toggle used to indicate whether the classification error
    of the tree should be calculated using the entire data set (when true) or just
    the training hold out set (when false)"""
    ...
