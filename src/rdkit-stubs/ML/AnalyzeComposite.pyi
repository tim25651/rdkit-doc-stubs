"""
rdkit.ML.AnalyzeComposite module¶
command line utility to report on the contributions of descriptors to
tree-based composite models
Usage:  AnalyzeComposite [optional args] <models>

<models>: file name(s) of pickled composite model(s)(this is the name of the db table if using a database)

Optional Arguments:

-n number: the number of levels of each model to consider
-d dbname: the database from which to read the models
-N Note: the note string to search for to pull models from the database
-v: be verbose whilst screening
"""
from rdkit.Dbase.DbConnection import DbConnect as DbConnect
from rdkit.ML import ScreenComposite as ScreenComposite
from rdkit.ML.Data import Stats as Stats
from rdkit.ML.DecTree import Tree as Tree
from rdkit.ML.DecTree import TreeUtils as TreeUtils

def ErrorStats(self, conn, where, enrich=1): ...
def ProcessIt(self, composites, nToConsider=3, verbose=0): ...
def ErrorStats(self, conn, where, enrich=1): ...
def ShowStats(self, statD, enrich=1): ...
def Usage(self): ...
