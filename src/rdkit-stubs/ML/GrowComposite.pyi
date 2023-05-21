"""
rdkit.ML.GrowComposite module¶
command line utility for growing composite models
Usage

_GrowComposite [optional args] filename_

Command Line Arguments

-n count: number of new models to build
-C pickle file name:  name of file containing composite upon which to build.

–inNote note: note to be used in loading composite models from the databasefor growing

–balTable table name:  table from which to take the original data set(for balancing)

–balWeight weight: (between 0 and 1) weighting factor for the new data(for balancing). OR, weight can be a list of weights

–balCnt count: number of individual models in the balanced composite(for balancing)

–balH: use only the holdout set from the original data set in the balancing(for balancing)

–balT: use only the training set from the original data set in the balancing(for balancing)

-S: shuffle the original data set(for balancing)

-r: randomize the activities of the original data set(for balancing)

-N note: note to be attached to the grown composite when it’s saved in thedatabase

–outNote note: equivalent to -N

-o filename: name of an output file to hold the pickled composite afterit has been grown.
If multiple balance weights are used, the weights will be added to
the filenames.

-L limit: provide an (integer) limit on individual model complexity

-d database name: instead of reading the data from a QDAT file,pull it from a database.  In this case, the _filename_ argument
provides the name of the database table containing the data set.

-p tablename: store persistence data in the databasein table tablename

-l: locks the random number generator to give consistent setsof training and hold-out data.  This is primarily intended
for testing purposes.

-g: be less greedy when training the models.
-G number: force trees to be rooted at descriptor number.

-D: show a detailed breakdown of the composite model performanceacross the training and, when appropriate, hold-out sets.

-t threshold value: use high-confidence predictions for the finalanalysis of the hold-out data.

-q list string:  Add QuantTrees to the composite and use the listspecified in list string as the number of target quantization
bounds for each descriptor.  Don’t forget to include 0’s at the
beginning and end of list string for the name and value fields.
For example, if there are 4 descriptors and you want 2 quant bounds
apiece, you would use _-q “[0,2,2,2,2,0]”_.
Two special cases:

If you would like to ignore a descriptor in the model building,
use ‘-1’ for its number of quant bounds.
If you have integer valued data that should not be quantized
further, enter 0 for that descriptor.

-V: print the version number and exit
"""
from _typeshed import Incomplete
from rdkit.Dbase.DbConnection import DbConnect as DbConnect
from rdkit.ML import BuildComposite as BuildComposite
from rdkit.ML import CompositeRun as CompositeRun
from rdkit.ML import ScreenComposite as ScreenComposite
from rdkit.ML.Composite import AdjustComposite as AdjustComposite
from rdkit.ML.Data import DataUtils as DataUtils
from rdkit.ML.Data import SplitData as SplitData

def BalanceComposite(self, details, composite, data1=None, data2=None):
    """
    balances the composite using the parameters provided in details
    Arguments

    details a _CompositeRun.RunDetails_ object
    composite: the composite model to be balanced
    data1: (optional) if provided, this should be the
    data set used to construct the original models
    data2: (optional) if provided, this should be the
    data set used to construct the new individual models"""
    ...

def GetComposites(self, details): ...
def message(self, msg):
    """
    emits messages to _sys.stdout_
    override this in modules which import this one to redirect output
    Arguments

    msg: the string to be displayed"""
    ...

def GrowIt(
    self, details, composite, progressCallback=None, saveIt=1, setDescNames=0, data=None
):
    """
    does the actual work of building a composite model
    Arguments

    details:  a _CompositeRun.CompositeRun_ object containing details
    (options, parameters, etc.) about the run
    composite: the composite model to grow
    progressCallback: (optional) a function which is called with a single
    argument (the number of models built so far) after each model is built.
    saveIt: (optional) if this is nonzero, the resulting model will be pickled
    and dumped to the filename specified in _details.outName_
    setDescNames: (optional) if nonzero, the composite’s _SetInputOrder()_ method
    will be called using the results of the data set’s _GetVarNames()_ method;
    it is assumed that the details object has a _descNames attribute which
    is passed to the composites _SetDescriptorNames()_ method.  Otherwise
    (the default), _SetDescriptorNames()_ gets the results of _GetVarNames()_.
    data: (optional) the data set to be used.  If this is not provided, the
    data set described in details will be used.

    Returns

    the enlarged composite model"""
    ...

def ParseArgs(self, runDetails):
    """
    parses command line arguments and updates _runDetails_
    Arguments

    runDetails:  a _CompositeRun.CompositeRun_ object."""
    ...

def SetDefaults(self, runDetails=None):
    """
    initializes a details object with default values
    Arguments

    details:  (optional) a _CompositeRun.CompositeRun_ object.
    If this is not provided, the global _runDetails will be used.

    Returns

    the initialized _CompositeRun_ object."""
    ...

def GetComposites(self, details): ...
def BalanceComposite(self, details, composite, data1=None, data2=None):
    """
    balances the composite using the parameters provided in details
    Arguments

    details a _CompositeRun.RunDetails_ object
    composite: the composite model to be balanced
    data1: (optional) if provided, this should be the
    data set used to construct the original models
    data2: (optional) if provided, this should be the
    data set used to construct the new individual models"""
    ...

def ShowVersion(self, includeArgs=0):
    """
    prints the version number"""
    ...

def Usage(self):
    """
    provides a list of arguments for when this is used from the command line"""
    ...

def message(self, msg):
    """
    emits messages to _sys.stdout_
    override this in modules which import this one to redirect output
    Arguments

    msg: the string to be displayed"""
    ...

def SetDefaults(self, runDetails=None):
    """
    initializes a details object with default values
    Arguments

    details:  (optional) a _CompositeRun.CompositeRun_ object.
    If this is not provided, the global _runDetails will be used.

    Returns

    the initialized _CompositeRun_ object."""
    ...

def ParseArgs(self, runDetails):
    """
    parses command line arguments and updates _runDetails_
    Arguments

    runDetails:  a _CompositeRun.CompositeRun_ object."""
    ...
