"""
rdkit.ML.EnrichPlot module¶
Command line tool to construct an enrichment plot from saved composite models
Usage:  EnrichPlot [optional args] -d dbname -t tablename <models>

Required Arguments:-d “dbName”: the name of the database for screening
-t “tablename”: provide the name of the table with the data to be screened

<models>: file name(s) of pickled composite model(s).If the -p argument is also provided (see below), this argument is ignored.

Optional Arguments:

-a “list”: the list of result codes to be considered active.  This will beeval’ed, so be sure that it evaluates as a list or sequence of
integers. For example, -a “[1,2]” will consider activity values 1 and 2
to be active

–enrich “list”: identical to the -a argument above.

–thresh: sets a threshold for the plot.  If the confidence falls belowthis value, picking will be terminated

-H: screen only the hold out set (works only if a version ofBuildComposite more recent than 1.2.2 was used).

-T: screen only the training set (works only if a version ofBuildComposite more recent than 1.2.2 was used).

-S: shuffle activity values before screening
-R: randomize activity values before screening

-F filter frac: filters the data before training to change thedistribution of activity values in the training set.  filter frac
is the fraction of the training set that should have the target value.
See note in BuildComposite help about data filtering

-v filter value: filters the data before training to change thedistribution of activity values in the training set. filter value
is the target value to use in filtering.
See note in BuildComposite help about data filtering

-p “tableName”: provides the name of a db table containing themodels to be screened.  If you use this argument, you should also
use the -N argument (below) to specify a note value.

-N “note”: provides a note to be used to pull models from a db table.
–plotFile “filename”: writes the data to an output text file (filename.dat)
and creates a gnuplot input file (filename.gnu) to plot it
–showPlot: causes the gnuplot plot constructed using –plotFile to be
displayed in gnuplot.
"""
import _io
from rdkit import DataStructs as DataStructs
from rdkit import RDConfig as RDConfig
from rdkit.Dbase.DbConnection import DbConnect as DbConnect
from rdkit.ML import CompositeRun as CompositeRun
from rdkit.ML.Data import DataUtils as DataUtils
from rdkit.ML.Data import SplitData as SplitData
from rdkit.ML.Data import Stats as Stats

def cmp(self, t1, t2): ...
def message(self, msg, noRet=0, dest: _io.TextIOWrapper = ...):
    """
    emits messages to _sys.stderr_
    override this in modules which import this one to redirect output
    Arguments

    msg: the string to be displayed"""
    ...

def error(self, msg, dest: _io.TextIOWrapper = ...):
    """
    emits messages to _sys.stderr_
    override this in modules which import this one to redirect output
    Arguments

    msg: the string to be displayed"""
    ...

def ScreenModel(self, mdl, descs, data, picking=[1], indices=[], errorEstimate=0):
    """
    collects the results of screening an individual composite model that matcha particular value

    Arguments

    mdl: the composite model
    descs: a list of descriptor names corresponding to the data set
    data: the data set, a list of points to be screened.
    picking: (Optional) a list of values that are to be collected.
    For examples, if you want an enrichment plot for picking the values
    1 and 2, you’d having picking=[1,2].

    Returns

    a list of 4-tuples containing:

    the id of the point
    the true result (from the data set)
    the predicted result
    the confidence value for the prediction"""
    ...

def AccumulateCounts(self, predictions, thresh=0, sortIt=1):
    """
    Accumulates the data for the enrichment plot for a single model
    Arguments

    predictions: a list of 3-tuples (as returned by _ScreenModels_)
    thresh: a threshold for the confidence level.  Anything below
    this threshold will not be considered
    sortIt: toggles sorting on confidence levels

    Returns

    a list of 3-tuples:

    the id of the active picked here
    num actives found so far
    number of picks made so far"""
    ...

def MakePlot(self, details, final, counts, pickVects, nModels, nTrueActs=-1): ...
def ScreenModel(self, mdl, descs, data, picking=[1], indices=[], errorEstimate=0):
    """
    collects the results of screening an individual composite model that matcha particular value

    Arguments

    mdl: the composite model
    descs: a list of descriptor names corresponding to the data set
    data: the data set, a list of points to be screened.
    picking: (Optional) a list of values that are to be collected.
    For examples, if you want an enrichment plot for picking the values
    1 and 2, you’d having picking=[1,2].

    Returns

    a list of 4-tuples containing:

    the id of the point
    the true result (from the data set)
    the predicted result
    the confidence value for the prediction"""
    ...

def Usage(self):
    """
    displays a usage message and exits"""
    ...

def cmp(self, t1, t2): ...
def error(self, msg, dest: _io.TextIOWrapper = ...):
    """
    emits messages to _sys.stderr_
    override this in modules which import this one to redirect output
    Arguments

    msg: the string to be displayed"""
    ...

def message(self, msg, noRet=0, dest: _io.TextIOWrapper = ...):
    """
    emits messages to _sys.stderr_
    override this in modules which import this one to redirect output
    Arguments

    msg: the string to be displayed"""
    ...
