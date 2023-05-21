"""
rdkit.ML.NaiveBayes.ClassificationModel module¶
Defines Naive Baysean classification model
Based on development in: Chapter 6 of “Machine Learning” by Tom Mitchell
"""
from _typeshed import Incomplete
from rdkit.ML.Data import Quantize as Quantize

class NaiveBayesClassifier(object):
    """
    _NaiveBayesClassifier_s can save the following pieces of internal state, accessible via
    standard setter/getter functions:

    _Examples_: a list of examples which have been predicted
    _TrainingExamples_: List of training examples - the descriptor value of these examples

    are quantized based on info gain using ML/Data/Quantize.py if necessary

    _TestExamples_: the list of examples used to test the model
    _BadExamples_ : list of examples that were incorrectly classified

    _QBoundVals_: Quant bound values for each varaible - a list of lists
    _QBounds_ : Number of bounds for each variable

    Constructor"""

    def ClassifyExample(self, example, appendExamples=0):
        """
        Classify an example by summing over the conditional probabilities
        The most likely class is the one with the largest probability"""
        ...
    def ClassifyExamples(self, examples, appendExamples=0): ...
    def GetBadExamples(self): ...
    def GetClassificationDetails(self):
        """
        returns the probability of the last prediction"""
        ...
    def GetExamples(self): ...
    mprob: Incomplete

    def __init__(
        self, attrs, nPossibleVals, nQuantBounds, mEstimateVal=..., useSigs: bool = ...
    ) -> None: ...
    def GetName(self): ...
    def GetTestExamples(self): ...
    def GetTrainingExamples(self): ...
    def SetName(self, name): ...
    def NameModel(self, varNames): ...
    def SetBadExamples(self, examples): ...
    def GetExamples(self): ...
    def SetExamples(self, examples): ...
    def SetName(self, name): ...
    def SetTestExamples(self, examples): ...
    def GetTrainingExamples(self): ...
    def SetTrainingExamples(self, examples): ...
    def GetTestExamples(self): ...
    def SetTestExamples(self, examples): ...
    def SetBadExamples(self, examples): ...
    def GetBadExamples(self): ...
    def trainModel(self):
        """
        We will assume at this point that the training examples have been set
        We have to estmate the conditional probabilities for each of the (binned) descriptor
        component give a outcome (or class). Also the probabilities for each class is estimated
        """
        ...
    def ClassifyExamples(self, examples, appendExamples=0): ...
    def GetClassificationDetails(self):
        """
        returns the probability of the last prediction"""
        ...
    def ClassifyExample(self, example, appendExamples=0):
        """
        Classify an example by summing over the conditional probabilities
        The most likely class is the one with the largest probability"""
        ...
