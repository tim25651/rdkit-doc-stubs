"""
rdkit.ML.MLUtils.VoteImg module¶
functionality for generating an image showing the results of a composite model
voting on a data set

Uses Numeric and PIL
"""

def BuildVoteImage(
    self,
    nModels,
    data,
    values,
    trueValues=[],
    sortTrueVals=0,
    xScale=10,
    yScale=2,
    addLine=1,
):
    """
    constructs the actual image
    Arguments

    nModels: the number of models in the composite
    data: the results of voting
    values: predicted values for each example
    trueValues: true values for each example
    sortTrueVals: if nonzero the votes will be sorted so
    that the _trueValues_ are in order, otherwise the sort
    is by _values_
    xScale: number of pixels per vote in the x direction
    yScale: number of pixels per example in the y direction

    addLine: if nonzero, a purple line is drawn separatingthe votes from the examples

    Returns

    a PIL image"""
    ...

def CollectVotes(self, composite, data, badOnly):
    """
    collects the votes from _composite_ for the examples in _data_

    Arguments

    composite: a composite model
    data: a list of examples to run through _composite_
    badOnly: if set only bad (misclassified) examples will be kept

    Returns

    a 4-tuple containing:

    the expanded list of vote details (see below)
    the list of predicted results
    the list of true results
    the number of miscounted examples

    Notes

    pp      - the expanded list of vote details consists of:

    ‘[ vote1, vote2, … voteN, 0, res, trueRes]’
    where _res_ is the predicted results and _trueRes_ is the actual result.
    The extra zero is included to allow a line to be drawn between the votes
    and the results."""
    ...

def Usage(self):
    """
    provides a list of arguments for when this is used from the command line"""
    ...

def BuildVoteImage(
    self,
    nModels,
    data,
    values,
    trueValues=[],
    sortTrueVals=0,
    xScale=10,
    yScale=2,
    addLine=1,
):
    """
    constructs the actual image
    Arguments

    nModels: the number of models in the composite
    data: the results of voting
    values: predicted values for each example
    trueValues: true values for each example
    sortTrueVals: if nonzero the votes will be sorted so
    that the _trueValues_ are in order, otherwise the sort
    is by _values_
    xScale: number of pixels per vote in the x direction
    yScale: number of pixels per example in the y direction

    addLine: if nonzero, a purple line is drawn separatingthe votes from the examples

    Returns

    a PIL image"""
    ...

def VoteAndBuildImage(
    self, composite, data, badOnly=0, sortTrueVals=0, xScale=10, yScale=2, addLine=1
):
    """
    collects votes on the examples and constructs an image
    Arguments

    composte: a composite model
    data: the examples to be voted upon
    badOnly: if nonzero only the incorrect votes will be shown
    sortTrueVals: if nonzero the votes will be sorted so
    that the _trueValues_ are in order, otherwise the sort
    is by _values_
    xScale: number of pixels per vote in the x direction
    yScale: number of pixels per example in the y direction

    addLine: if nonzero, a purple line is drawn separatingthe votes from the examples

    Returns

    a PIL image"""
    ...

def Usage(self):
    """
    provides a list of arguments for when this is used from the command line"""
    ...
