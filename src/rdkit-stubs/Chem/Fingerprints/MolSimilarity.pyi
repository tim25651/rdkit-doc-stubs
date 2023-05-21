"""
rdkit.Chem.Fingerprints.MolSimilarity module¶

utility functionality for molecular similarityincludes a command line app for screening databases

Sample Usage:

python MolSimilarity.py  -d data.gdb -t daylight_sig –idName=”Mol_ID”       –topN=100 –smiles=’c1(C=O)ccc(Oc2ccccc2)cc1’ –smilesTable=raw_dop_data       –smilesName=”structure” -o results.csv
"""
from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit import DataStructs as DataStructs
from rdkit.Chem.Fingerprints import DbFpSupplier as DbFpSupplier
from rdkit.Chem.Fingerprints import FingerprintMols as FingerprintMols
from rdkit.DataStructs.TopNContainer import TopNContainer as TopNContainer
from rdkit.Dbase import DbModule as DbModule
from rdkit.Dbase.DbConnection import DbConnect as DbConnect

def ScreenInDb(self, details, mol): ...
def GetFingerprints(self, details):
    """
    returns an iterable sequence of fingerprints
    each fingerprint will have a _fieldsFromDb member whose first entry is
    the id."""
    ...

def ScreenFingerprints(self, details, data, mol=None, probeFp=None):
    """
    Returns a list of results"""
    ...

def ScreenFromDetails(self, details, mol=None):
    """
    Returns a list of results"""
    ...

def ScreenInDb(self, details, mol): ...
