"""
rdkit.Chem.MolKey.MolKey module¶
"""
from typing import NamedTuple, _Alias

from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit import RDConfig as RDConfig
from rdkit.Avalon import pyAvalonTools as pyAvalonTools
from rdkit.Chem.MolKey import InchiInfo as InchiInfo

class InchiResult(tuple):
    """
    Create new instance of InchiResult(error, inchi, fixed_ctab)"""

    error: _Alias = ...
    fixed_ctab: _Alias = ...
    inchi: _Alias = ...
    fixed_ctab: _Alias = ...

class MolKeyResult(tuple):
    """
    Create new instance of MolKeyResult(mol_key, error, inchi, fixed_ctab, stereo_code, stereo_comment)
    """

    mol_key: _Alias = ...
    error: _Alias = ...
    fixed_ctab: _Alias = ...
    inchi: _Alias = ...
    mol_key: _Alias = ...
    fixed_ctab: _Alias = ...
    stereo_code: _Alias = ...
    stereo_comment: _Alias = ...

class BadMoleculeException(Exception):
    ...
    ...

class MolIdentifierException(Exception):
    ...
    ...

class BadMoleculeException(Exception):
    ...
    ...

MOL_KEY_VERSION: str
ERROR_DICT: Incomplete
INCHI_COMPUTATION_ERROR: Incomplete
RDKIT_CONVERSION_ERROR: Incomplete
INCHI_READWRITE_ERROR: Incomplete
NULL_MOL: Incomplete
BAD_SET: Incomplete
GET_STEREO_RE: Incomplete
NULL_SMILES_RE: Incomplete
PATTERN_NULL_MOL: str
CHIRAL_POS: int
T_NULL_MOL: Incomplete
stereo_code_dict: Incomplete

def initStruchk(self, configDir=None, logFile=None): ...
def CheckCTAB(self, ctab, isSmiles=True): ...
def ErrorBitsToText(self, err):
    """
    returns a list of error bit descriptions for the error code provided"""
    ...

class InchiResult(tuple):
    """
    Create new instance of InchiResult(error, inchi, fixed_ctab)"""

    error: _Alias = ...
    fixed_ctab: _Alias = ...
    inchi: _Alias = ...
    fixed_ctab: _Alias = ...

def GetInchiForCTAB(self, ctab):
    """
    >>> from rdkit.Chem.MolKey import MolKey
    >>> from rdkit.Avalon import pyAvalonTools
    >>> res = MolKey.GetInchiForCTAB(pyAvalonTools.Generate2DCoords('c1cn[nH]c1C(Cl)Br',True))
    >>> res.inchi
    'InChI=1/C4H4BrClN2/c5-4(6)3-1-2-7-8-3/h1-2,4H,(H,7,8)/t4?/f/h8H'
    >>> res = MolKey.GetInchiForCTAB(pyAvalonTools.Generate2DCoords('c1c[nH]nc1C(Cl)Br',True))
    >>> res.inchi
    'InChI=1/C4H4BrClN2/c5-4(6)3-1-2-7-8-3/h1-2,4H,(H,7,8)/t4?/f/h7H'
    >>>"""
    ...

def ErrorBitsToText(self, err):
    """
    returns a list of error bit descriptions for the error code provided"""
    ...

class MolKeyResult(tuple):
    """
    Create new instance of MolKeyResult(mol_key, error, inchi, fixed_ctab, stereo_code, stereo_comment)
    """

    mol_key: _Alias = ...
    error: _Alias = ...
    fixed_ctab: _Alias = ...
    inchi: _Alias = ...
    mol_key: _Alias = ...
    fixed_ctab: _Alias = ...
    stereo_code: _Alias = ...
    stereo_comment: _Alias = ...

def GetKeyForCTAB(self, ctab, stereo_info=None, stereo_comment=None, logger=None):
    """
    >>> from rdkit.Chem.MolKey import MolKey
    >>> from rdkit.Avalon import pyAvalonTools
    >>> res = MolKey.GetKeyForCTAB(pyAvalonTools.Generate2DCoords('c1ccccc1C(F)Cl', True))
    >>> res.mol_key
    '1|L7676nfGsSIU33wkx//NCg=='
    >>> res.stereo_code
    'R_ONE'
    >>> res = MolKey.GetKeyForCTAB(pyAvalonTools.Generate2DCoords('c1ccccc1[C@H](F)Cl', True))
    >>> res.mol_key
    '1|Aj38EIxf13RuPDQG2A0UMw=='
    >>> res.stereo_code
    'S_ABS'
    >>> res = MolKey.GetKeyForCTAB(pyAvalonTools.Generate2DCoords('c1ccccc1[C@@H](F)Cl', True))
    >>> res.mol_key
    '1|9ypfMrhxn1w0ncRooN5HXw=='
    >>> res.stereo_code
    'S_ABS'
    >>> res = MolKey.GetKeyForCTAB(pyAvalonTools.Generate2DCoords('c1cccc(C(Br)Cl)c1[C@@H](F)Cl', True))
    >>> res.mol_key
    '1|c96jMSlbn7O9GW5d5uB9Mw=='
    >>> res.stereo_code
    'S_PART'
    >>> res = MolKey.GetKeyForCTAB(pyAvalonTools.Generate2DCoords('c1cccc([C@H](Br)Cl)c1[C@@H](F)Cl', True))
    >>> res.mol_key
    '1|+B+GCEardrJteE8xzYdGLA=='
    >>> res.stereo_code
    'S_ABS'
    >>> res = MolKey.GetKeyForCTAB(pyAvalonTools.Generate2DCoords('c1cccc(C(Br)Cl)c1C(F)Cl', True))
    >>> res.mol_key
    '1|5H9R3LvclagMXHp3Clrc/g=='
    >>> res.stereo_code
    'S_UNKN'
    >>> res = MolKey.GetKeyForCTAB(pyAvalonTools.Generate2DCoords('c1cccc(C(Br)Cl)c1C(F)Cl',True), stereo_info='S_REL')
    >>> res.mol_key
    '1|cqKWVsUEY6QNpGCbDaDTYA=='
    >>> res.stereo_code
    'S_REL'
    >>> res.inchi
    'InChI=1/C8H6BrCl2F/c9-7(10)5-3-1-2-4-6(5)8(11)12/h1-4,7-8H/t7?,8?'"""
    ...

def initStruchk(self, configDir=None, logFile=None): ...
