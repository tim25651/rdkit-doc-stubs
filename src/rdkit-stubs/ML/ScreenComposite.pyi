"""
rdkit.ML.ScreenComposite module¶
command line utility for screening composite models
Usage

_ScreenComposite [optional args] modelfile(s) datafile_

Unless indicated otherwise (via command line arguments), _modelfile_ is
a file containing a pickled composite model and _filename_ is a QDAT file.
Command Line Arguments

-t threshold value(s): use high-confidence predictions for the finalanalysis of the hold-out data.  The threshold value can be either a single
float or a list/tuple of floats.  All thresholds should be between
0.0 and 1.0

-D: do a detailed screen.

-d database name: instead of reading the data from a QDAT file,pull it from a database.  In this case, the _datafile_ argument
provides the name of the database table containing the data set.

-N note: use all models from the database which have this note.The modelfile argument should contain the name of the table
with the models.

-H: screen only the hold out set (works only if a version ofBuildComposite more recent than 1.2.2 was used).

-T: screen only the training set (works only if a version ofBuildComposite more recent than 1.2.2 was used).

-E: do a detailed Error analysis.  This shows each misclassifiedpoint and the number of times it was missed across all screened
composites.  If the –enrich argument is also provided, only compounds
that have true activity value equal to the enrichment value will be
used.

–enrich enrichVal: target “active” value to be used in calculatingenrichments.

-A: show All predictions.
-S: shuffle activity values before screening
-R: randomize activity values before screening

-F filter frac: filters the data before training to change thedistribution of activity values in the training set.  filter frac
is the fraction of the training set that should have the target value.
See note in BuildComposite help about data filtering

-v filter value: filters the data before training to change thedistribution of activity values in the training set. filter value
is the target value to use in filtering.
See note in BuildComposite help about data filtering

-V: be verbose when screening multiple models
-h: show this message and exit

–OOB: Do out an “out-of-bag” generalization error estimate.  This onlymakes sense when applied to the original data set.

–pickleCol colId: index of the column containing a pickled value(used primarily for cases where fingerprints are used as descriptors)

* Options for making Prediction (Hanneke) Plots *

–predPlot=<fileName>: triggers the generation of a Hanneke plot andsets the name of the .txt file which will hold the output data.
A Gnuplot control file, <fileName>.gnu, will also be generated.

–predActTable=<name> (optional):  name of the database tablecontaining activity values.  If this is not provided, activities
will be read from the same table containing the screening data

–predActCol=<name> (optional):  name of the activity column. If notprovided, the name of the last column in the activity table will
be used.

–predLogScale (optional):  If provided, the x axis of theprediction plot (the activity axis) will be plotted using a log
scale

–predShow: launch a gnuplot instance and display the predictionplot (the plot will still be written to disk).

* The following options are likely obsolete *

-P: read pickled data.  The datafile argument should containa pickled data set. relevant only to qdat files

-q: data are not quantized (the composite should take care ofquantization itself if it requires quantized data). relevant only to
qdat files
"""
from _typeshed import Incomplete
from rdkit import DataStructs as DataStructs
from rdkit.Dbase import DbModule as DbModule
from rdkit.Dbase.DbConnection import DbConnect as DbConnect
from rdkit.ML import CompositeRun as CompositeRun
from rdkit.ML.Data import DataUtils as DataUtils
from rdkit.ML.Data import SplitData as SplitData

hasPil: int

def message(self, msg, noRet=0):
    """
    emits messages to _sys.stdout_
    override this in modules which import this one to redirect output
    Arguments

    msg: the string to be displayed"""
    ...

def error(self, msg):
    """
    emits messages to _sys.stderr_
    override this in modules which import this one to redirect output
    Arguments

    msg: the string to be displayed"""
    ...

def CalcEnrichment(self, mat, tgt=1): ...
def CollectResults(
    self, indices, dataSet, composite, callback=None, appendExamples=0, errorEstimate=0
):
    """
    screens a set of examples through a composite and returns theresults

    #DOC

    Arguments

    examples: the examples to be screened (a sequence of sequences)it’s assumed that the last element in each example is it’s “value”

    composite:  the composite model to be used
    callback: (optional)  if provided, this should be a function
    taking a single argument that is called after each example is
    screened with the number of examples screened so far as the
    argument.
    appendExamples: (optional)  this value is passed on to the
    composite’s _ClassifyExample()_ method.
    errorEstimate: (optional) calculate the “out of bag” error
    estimate for the composite using Breiman’s definition.  This
    only makes sense when screening the original data set!
    [L. Breiman “Out-of-bag Estimation”, UC Berkeley Dept of
    Statistics Technical Report (1996)]

    Returns

    a list of 3-tuples _nExamples_ long:

    answer: the value from the example
    pred: the composite model’s prediction
    conf: the confidence of the composite"""
    ...

def DetailedScreen(
    self,
    indices,
    data,
    composite,
    threshold=0,
    screenResults=None,
    goodVotes=None,
    badVotes=None,
    noVotes=None,
    callback=None,
    appendExamples=0,
    errorEstimate=0,
):
    """
    screens a set of examples cross a composite and breaks thepredictions into correct,*incorrect* and unclassified sets.

    #DOCArguments

    examples: the examples to be screened (a sequence of sequences)it’s assumed that the last element in each example is its “value”

    composite:  the composite model to be used
    threshold: (optional) the threshold to be used to decide whether
    or not a given prediction should be kept
    screenResults: (optional) the results of screening the results
    (a sequence of 3-tuples in the format returned by
    _CollectResults()_).  If this is provided, the examples will not
    be screened again.
    goodVotes,badVotes,noVotes: (optional)  if provided these should
    be lists (or anything supporting an _append()_ method) which
    will be used to pass the screening results back.
    callback: (optional)  if provided, this should be a function
    taking a single argument that is called after each example is
    screened with the number of examples screened so far as the
    argument.
    appendExamples: (optional)  this value is passed on to the
    composite’s _ClassifyExample()_ method.
    errorEstimate: (optional) calculate the “out of bag” error
    estimate for the composite using Breiman’s definition.  This
    only makes sense when screening the original data set!
    [L. Breiman “Out-of-bag Estimation”, UC Berkeley Dept of
    Statistics Technical Report (1996)]

    Notes

    since this function doesn’t return anything, if one or more of
    the arguments _goodVotes_, _badVotes_, and _noVotes_ is not
    provided, there’s not much reason to call it"""
    ...

def GetScreenImage(self, nGood, nBad, nRej, size=None): ...
def Go(self, details): ...
def MakePredPlot(
    self, details, indices, data, goodVotes, badVotes, nRes, idCol=0, verbose=0
):
    """
    Arguments

    details:  a CompositeRun.RunDetails object
    indices: a sequence of integer indices into _data_
    data: the data set in question.  We assume that the ids for
    the data points are in the _idCol_ column
    goodVotes/badVotes: predictions where the model was correct/incorrect.
    These are sequences of 4-tuples:

    (answer,prediction,confidence,index into _indices_)"""
    ...

def ParseArgs(self, details): ...
def ShowVoteResults(
    self,
    indices,
    data,
    composite,
    nResultCodes,
    threshold,
    verbose=1,
    screenResults=None,
    callback=None,
    appendExamples=0,
    goodVotes=None,
    badVotes=None,
    noVotes=None,
    errorEstimate=0,
):
    """
    screens the results and shows a detailed workup

    The work of doing the screening and processing the results is
    handled by _DetailedScreen()_

    #DOC

    Arguments

    examples: the examples to be screened (a sequence of sequences)it’s assumed that the last element in each example is its “value”

    composite:  the composite model to be used
    nResultCodes: the number of possible results the composite can
    return
    threshold: the threshold to be used to decide whether or not a
    given prediction should be kept
    screenResults: (optional) the results of screening the results
    (a sequence of 3-tuples in the format returned by
    _CollectResults()_).  If this is provided, the examples will not
    be screened again.
    callback: (optional)  if provided, this should be a function
    taking a single argument that is called after each example is
    screened with the number of examples screened so far as the
    argument.
    appendExamples: (optional)  this value is passed on to the
    composite’s _ClassifyExample()_ method.
    goodVotes,badVotes,noVotes: (optional)  if provided these should
    be lists (or anything supporting an _append()_ method) which
    will be used to pass the screening results back.
    errorEstimate: (optional) calculate the “out of bag” error
    estimate for the composite using Breiman’s definition.  This
    only makes sense when screening the original data set!
    [L. Breiman “Out-of-bag Estimation”, UC Berkeley Dept of
    Statistics Technical Report (1996)]

    Returns

    a 7-tuple:

    the number of good (correct) predictions
    the number of bad (incorrect) predictions
    the number of predictions skipped due to the _threshold_
    the average confidence in the good predictions
    the average confidence in the bad predictions
    the average confidence in the skipped predictions
    the results table"""
    ...

def ScreenIt(
    self,
    composite,
    indices,
    data,
    partialVote=0,
    voteTol=0.0,
    verbose=1,
    screenResults=None,
    goodVotes=None,
    badVotes=None,
    noVotes=None,
):
    """
    screens a set of data using a composite model and prints outstatistics about the screen.

    #DOC
    The work of doing the screening and processing the results is
    handled by _DetailedScreen()_

    Arguments

    composite:  the composite model to be used

    data: the examples to be screened (a sequence of sequences)it’s assumed that the last element in each example is its “value”

    partialVote: (optional) toggles use of the threshold value in
    the screnning.
    voteTol: (optional) the threshold to be used to decide whether or not a
    given prediction should be kept
    verbose: (optional) sets degree of verbosity of the screening
    screenResults: (optional) the results of screening the results
    (a sequence of 3-tuples in the format returned by
    _CollectResults()_).  If this is provided, the examples will not
    be screened again.
    goodVotes,badVotes,noVotes: (optional)  if provided these should
    be lists (or anything supporting an _append()_ method) which
    will be used to pass the screening results back.

    Returns

    a 7-tuple:

    the number of good (correct) predictions
    the number of bad (incorrect) predictions
    the number of predictions skipped due to the _threshold_
    the average confidence in the good predictions
    the average confidence in the bad predictions
    the average confidence in the skipped predictions
    None"""
    ...

def PrepareDataFromDetails(self, model, details, data, verbose=0): ...
def ScreenFromDetails(
    self,
    models,
    details,
    callback=None,
    setup=None,
    appendExamples=0,
    goodVotes=None,
    badVotes=None,
    noVotes=None,
    data=None,
    enrichments=None,
):
    """
    Screens a set of data using a a _CompositeRun.CompositeRun_instance to provide parameters

    # DOC

    The actual data to be used are extracted from the database and
    table specified in _details_
    Aside from dataset construction,  _ShowVoteResults()_ does most of
    the heavy lifting here.

    Arguments

    model: a composite model
    details:  a _CompositeRun.CompositeRun_ object containing details
    (options, parameters, etc.) about the run
    callback: (optional)  if provided, this should be a function
    taking a single argument that is called after each example is
    screened with the number of examples screened so far as the
    argument.
    setup: (optional) a function taking a single argument which is
    called at the start of screening with the number of points to
    be screened as the argument.
    appendExamples: (optional)  this value is passed on to the
    composite’s _ClassifyExample()_ method.
    goodVotes,badVotes,noVotes: (optional)  if provided these should
    be lists (or anything supporting an _append()_ method) which
    will be used to pass the screening results back.

    Returns

    a 7-tuple:

    the number of good (correct) predictions
    the number of bad (incorrect) predictions
    the number of predictions skipped due to the _threshold_
    the average confidence in the good predictions
    the average confidence in the bad predictions
    the average confidence in the skipped predictions
    the results table"""
    ...

def ScreenIt(
    self,
    composite,
    indices,
    data,
    partialVote=0,
    voteTol=0.0,
    verbose=1,
    screenResults=None,
    goodVotes=None,
    badVotes=None,
    noVotes=None,
):
    """
    screens a set of data using a composite model and prints outstatistics about the screen.

    #DOC
    The work of doing the screening and processing the results is
    handled by _DetailedScreen()_

    Arguments

    composite:  the composite model to be used

    data: the examples to be screened (a sequence of sequences)it’s assumed that the last element in each example is its “value”

    partialVote: (optional) toggles use of the threshold value in
    the screnning.
    voteTol: (optional) the threshold to be used to decide whether or not a
    given prediction should be kept
    verbose: (optional) sets degree of verbosity of the screening
    screenResults: (optional) the results of screening the results
    (a sequence of 3-tuples in the format returned by
    _CollectResults()_).  If this is provided, the examples will not
    be screened again.
    goodVotes,badVotes,noVotes: (optional)  if provided these should
    be lists (or anything supporting an _append()_ method) which
    will be used to pass the screening results back.

    Returns

    a 7-tuple:

    the number of good (correct) predictions
    the number of bad (incorrect) predictions
    the number of predictions skipped due to the _threshold_
    the average confidence in the good predictions
    the average confidence in the bad predictions
    the average confidence in the skipped predictions
    None"""
    ...

def GetScreenImage(self, nGood, nBad, nRej, size=None): ...
def ScreenToHtml(
    self,
    nGood,
    nBad,
    nRej,
    avgGood,
    avgBad,
    avgSkip,
    voteTable,
    imgDir=".",
    fullPage=1,
    skipImg=0,
    includeDefs=1,
):
    """
    returns the text of a web page showing the screening details
    #DOC

    Arguments

    nGood: number of correct predictions
    nBad:  number of incorrect predictions
    nRej:  number of rejected predictions
    avgGood: average correct confidence
    avgBad:  average incorrect confidence
    avgSkip: average rejected confidence
    voteTable: vote table
    imgDir: (optional) the directory to be used to hold the vote
    image (if constructed)

    Returns

    a string containing HTML"""
    ...

def MakePredPlot(
    self, details, indices, data, goodVotes, badVotes, nRes, idCol=0, verbose=0
):
    """
    Arguments

    details:  a CompositeRun.RunDetails object
    indices: a sequence of integer indices into _data_
    data: the data set in question.  We assume that the ids for
    the data points are in the _idCol_ column
    goodVotes/badVotes: predictions where the model was correct/incorrect.
    These are sequences of 4-tuples:

    (answer,prediction,confidence,index into _indices_)"""
    ...

def Go(self, details): ...
def SetDefaults(self, details=None): ...
def ShowVersion(self, includeArgs=0):
    """
    prints the version number of the program"""
    ...

def ShowVoteResults(
    self,
    indices,
    data,
    composite,
    nResultCodes,
    threshold,
    verbose=1,
    screenResults=None,
    callback=None,
    appendExamples=0,
    goodVotes=None,
    badVotes=None,
    noVotes=None,
    errorEstimate=0,
):
    """
    screens the results and shows a detailed workup

    The work of doing the screening and processing the results is
    handled by _DetailedScreen()_

    #DOC

    Arguments

    examples: the examples to be screened (a sequence of sequences)it’s assumed that the last element in each example is its “value”

    composite:  the composite model to be used
    nResultCodes: the number of possible results the composite can
    return
    threshold: the threshold to be used to decide whether or not a
    given prediction should be kept
    screenResults: (optional) the results of screening the results
    (a sequence of 3-tuples in the format returned by
    _CollectResults()_).  If this is provided, the examples will not
    be screened again.
    callback: (optional)  if provided, this should be a function
    taking a single argument that is called after each example is
    screened with the number of examples screened so far as the
    argument.
    appendExamples: (optional)  this value is passed on to the
    composite’s _ClassifyExample()_ method.
    goodVotes,badVotes,noVotes: (optional)  if provided these should
    be lists (or anything supporting an _append()_ method) which
    will be used to pass the screening results back.
    errorEstimate: (optional) calculate the “out of bag” error
    estimate for the composite using Breiman’s definition.  This
    only makes sense when screening the original data set!
    [L. Breiman “Out-of-bag Estimation”, UC Berkeley Dept of
    Statistics Technical Report (1996)]

    Returns

    a 7-tuple:

    the number of good (correct) predictions
    the number of bad (incorrect) predictions
    the number of predictions skipped due to the _threshold_
    the average confidence in the good predictions
    the average confidence in the bad predictions
    the average confidence in the skipped predictions
    the results table"""
    ...

def Usage(self):
    """
    prints a list of arguments for when this is used from the
    command line and then exits"""
    ...

def error(self, msg):
    """
    emits messages to _sys.stderr_
    override this in modules which import this one to redirect output
    Arguments

    msg: the string to be displayed"""
    ...

def message(self, msg, noRet=0):
    """
    emits messages to _sys.stdout_
    override this in modules which import this one to redirect output
    Arguments

    msg: the string to be displayed"""
    ...

def ShowVersion(self, includeArgs=0):
    """
    prints the version number of the program"""
    ...

def ParseArgs(self, details): ...
