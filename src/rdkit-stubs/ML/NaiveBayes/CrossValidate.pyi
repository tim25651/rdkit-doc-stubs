"""
rdkit.ML.NaiveBayes.CrossValidate module¶
handles doing cross validation with naive bayes models
and evaluation of individual models
"""
from _typeshed import Incomplete
from rdkit.ML.Data import SplitData as SplitData
from rdkit.ML.FeatureSelect import CMIM as CMIM
from rdkit.ML.NaiveBayes.ClassificationModel import (
    NaiveBayesClassifier as NaiveBayesClassifier,
)

def makeNBClassificationModel(
    self,
    trainExamples,
    attrs,
    nPossibleValues,
    nQuantBounds,
    mEstimateVal=-1.0,
    useSigs=False,
    ensemble=None,
    useCMIM=0,
    **kwargs,
): ...
def CrossValidate(self, NBmodel, testExamples, appendExamples=0): ...
def CrossValidationDriver(
    self,
    examples,
    attrs,
    nPossibleValues,
    nQuantBounds,
    mEstimateVal=0.0,
    holdOutFrac=0.3,
    modelBuilder=makeNBClassificationModel,
    silent=0,
    calcTotalError=0,
    **kwargs,
): ...
def makeNBClassificationModel(
    self,
    trainExamples,
    attrs,
    nPossibleValues,
    nQuantBounds,
    mEstimateVal=-1.0,
    useSigs=False,
    ensemble=None,
    useCMIM=0,
    **kwargs,
): ...
