"""
identify acyclic bonds
enumerate all single cuts
make sure you chop off more that 1 atom
keeps bits which are >60% query mol
enumerate all double cuts
keeps bits with 1 attachment point (i.e throw middle bit away)
need to be >60% query mol
identify exocyclic bonds
enumerate all single “ring” cuts
Check if it results in more that one component
keep correct bit if >40% query mol
enumerate successful “rings” cuts with an acyclic cut
Check if it results in more that one component
keep correct if >60% query mol
"""
from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit import DataStructs as DataStructs
from rdkit.Chem import rdqueries as rdqueries

rdkitFpParams: Incomplete
FTYPE_ACYCLIC: str
FTYPE_CYCLIC: str
FTYPE_CYCLIC_ACYCLIC: str
ACYC_SMARTS: Incomplete
CYC_SMARTS: Incomplete
cSma1: Incomplete
cSma2: Incomplete
dummyAtomQuery: Incomplete

def delete_bonds(mol, bonds, ftype, hac): ...
def select_fragments(fragments, ftype, hac): ...
def isValidRingCut(mol): ...
def generate_fraggle_fragmentation(mol, verbose: bool = ...): ...
def atomContrib(subs, mol, tverskyThresh: float = ...): ...

modified_query_fps: Incomplete

def compute_fraggle_similarity_for_subs(
    inMol, qMol, qSmi, qSubs, tverskyThresh: float = ...
): ...
def GetFraggleSimilarity(queryMol, refMol, tverskyThresh: float = ...): ...
