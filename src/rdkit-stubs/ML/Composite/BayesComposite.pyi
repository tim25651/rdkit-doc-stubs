"""
rdkit.ML.Composite.BayesComposite module¶
code for dealing with Bayesian composite models
For a model to be useable here, it should support the following API:

_ClassifyExample(example)_, returns a classification

Other compatibility notes:

To use _Composite.Grow_ there must be some kind of builder
functionality which returns a 2-tuple containing (model,percent accuracy).
The models should be pickleable
It would be very happy if the models support the __cmp__ method so that
membership tests used to make sure models are unique work.
"""
from _typeshed import Incomplete
from rdkit.ML.Composite import Composite
from rdkit.ML.Composite import Composite as Composite

class BayesComposite(Composite.Composite):
    """
    a composite model using Bayesian statistics in the Decision Proxy
    Notes

    typical usage:

    grow the composite with AddModel until happy with it
    call AverageErrors to calculate the average error values
    call SortModels to put things in order by either error or count
    call Train to update the Bayesian stats."""

    def ClassifyExample(self, example, threshold=0, verbose=0, appendExample=0):
        """
        classifies the given example using the entire composite
        Arguments

        example: the data to be classified

        threshold:  if this is a number greater than zero, then aclassification will only be returned if the confidence is
        above _threshold_.  Anything lower is returned as -1.

        Returns

        a (result,confidence) tuple"""
        ...
    resultProbs: Incomplete
    condProbs: Incomplete

    def Train(self, data, verbose=0): ...
    modelVotes: Incomplete

    def ClassifyExample(self, example, threshold=0, verbose=0, appendExample=0):
        """
        classifies the given example using the entire composite
        Arguments

        example: the data to be classified

        threshold:  if this is a number greater than zero, then aclassification will only be returned if the confidence is
        above _threshold_.  Anything lower is returned as -1.

        Returns

        a (result,confidence) tuple"""
        ...
    def __init__(self) -> None: ...

def BayesCompositeToComposite(self, obj):
    """
    converts a BayesComposite to a Composite.Composite"""
    ...

def CompositeToBayesComposite(self, obj):
    """
    converts a Composite to a BayesComposite

    if _obj_ is already a BayesComposite or if it is not a _Composite.Composite_ ,nothing will be done.
    """
    ...

def BayesCompositeToComposite(self, obj):
    """
    converts a BayesComposite to a Composite.Composite"""
    ...
