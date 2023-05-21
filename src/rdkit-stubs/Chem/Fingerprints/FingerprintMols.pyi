"""
rdkit.Chem.Fingerprints.FingerprintMols module¶

utility functionality for fingerprinting sets of moleculesincludes a command line app for working with fingerprints
and databases

Sample Usage:

python FingerprintMols.py  -d data.gdb         -t ‘raw_dop_data’ –smilesName=”Structure” –idName=”Mol_ID”          –outTable=”daylight_sig”
"""
from typing import Callable

from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit import DataStructs as DataStructs
from rdkit.Chem import MACCSkeys as MACCSkeys
from rdkit.ML.Cluster import Murtagh as Murtagh

class FingerprinterDetails(object):
    """
    class for storing the details of a fingerprinting run,
    generates sensible defaults on construction"""

    def __init__(self) -> None: ...
    def GetMetricName(self): ...
    metric: Incomplete

    def SetMetricFromName(self, name): ...

def error(self, msg): ...
def message(self, msg): ...
def GetRDKFingerprint(self, mol):
    """
    uses default parameters"""
    ...

def FoldFingerprintToTargetDensity(self, fp, **fpArgs): ...
def FingerprintMol(self, mol, fingerprinter: Callable = ..., **fpArgs): ...
def FingerprintsFromDetails(self, details, reportFreq=10): ...
def FingerprintsFromSmiles(
    self,
    dataSource,
    idCol,
    smiCol,
    fingerprinter: Callable = ...,
    reportFreq=10,
    maxMols=-1,
    **fpArgs,
):
    """
    fpArgs are passed as keyword arguments to the fingerprinter
    Returns a list of 2-tuples: (ID,fp)"""
    ...

def FingerprintsFromMols(
    self, mols, fingerprinter: Callable = ..., reportFreq=10, maxMols=-1, **fpArgs
):
    """
    fpArgs are passed as keyword arguments to the fingerprinter
    Returns a list of 2-tuples: (ID,fp)"""
    ...

def FingerprintsFromPickles(
    self,
    dataSource,
    idCol,
    pklCol,
    fingerprinter: Callable = ...,
    reportFreq=10,
    maxMols=-1,
    **fpArgs,
):
    """
    fpArgs are passed as keyword arguments to the fingerprinter
    Returns a list of 2-tuples: (ID,fp)"""
    ...

def FingerprintsFromSmiles(
    self,
    dataSource,
    idCol,
    smiCol,
    fingerprinter: Callable = ...,
    reportFreq=10,
    maxMols=-1,
    **fpArgs,
):
    """
    fpArgs are passed as keyword arguments to the fingerprinter
    Returns a list of 2-tuples: (ID,fp)"""
    ...

def FoldFingerprintToTargetDensity(self, fp, **fpArgs): ...
def GetRDKFingerprint(self, mol):
    """
    uses default parameters"""
    ...

def ParseArgs(self, details=None):
    """
    parses the command line arguments and returns a
    _FingerprinterDetails_ instance with the results.
    Note:

    If you make modifications here, please update the global
    _usageDoc string so the Usage message is up to date.
    This routine is used by both the fingerprinter, the clusterer and the
    screener; not all arguments make sense for all applications."""
    ...

def FingerprintsFromDetails(self, details, reportFreq=10): ...

class FingerprinterDetails(object):
    """
    class for storing the details of a fingerprinting run,
    generates sensible defaults on construction"""

    def __init__(self) -> None: ...
    def GetMetricName(self): ...
    metric: Incomplete

    def SetMetricFromName(self, name): ...

def Usage(self):
    """
    prints a usage string and exits"""
    ...

def error(self, msg): ...
def message(self, msg): ...
def ParseArgs(self, details=None):
    """
    parses the command line arguments and returns a
    _FingerprinterDetails_ instance with the results.
    Note:

    If you make modifications here, please update the global
    _usageDoc string so the Usage message is up to date.
    This routine is used by both the fingerprinter, the clusterer and the
    screener; not all arguments make sense for all applications."""
    ...
