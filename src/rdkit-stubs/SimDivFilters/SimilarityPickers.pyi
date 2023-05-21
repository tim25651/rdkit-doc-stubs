"""
rdkit.SimDivFilters.SimilarityPickers module¶
"""
from _typeshed import Incomplete
from rdkit import DataStructs as DataStructs
from rdkit.DataStructs.TopNContainer import TopNContainer as TopNContainer

class GenericPicker(object):
    def MakePicks(self, force=False): ...
    def __len__(self) -> int: ...
    def __getitem__(self, which): ...

class SpreadPicker(GenericPicker):
    """
    A class for picking the best matches across a library
    Connect to a database:
    >>> from rdkit import Chem
    >>> from rdkit import RDConfig
    >>> import os.path
    >>> from rdkit.Dbase.DbConnection import DbConnect
    >>> dbName = RDConfig.RDTestDatabase
    >>> conn = DbConnect(dbName,'simple_mols1')
    >>> [x.upper() for x in conn.GetColumnNames()]
    ['SMILES', 'ID']
    >>> mols = []
    >>> for smi,id in conn.GetData():
    ...   mol = Chem.MolFromSmiles(str(smi))
    ...   mol.SetProp('_Name',str(id))
    ...   mols.append(mol)
    >>> len(mols)
    12

    Calculate fingerprints:
    >>> probefps = []
    >>> for mol in mols:
    ...   fp = Chem.RDKFingerprint(mol)
    ...   fp._id = mol.GetProp('_Name')
    ...   probefps.append(fp)

    Start by finding the top matches for a single probe.  This ether should pull
    other ethers from the db:
    >>> mol = Chem.MolFromSmiles('COC')
    >>> probeFp = Chem.RDKFingerprint(mol)
    >>> picker = SpreadPicker(numToPick=2,probeFps=[probeFp],dataSet=probefps)
    >>> len(picker)
    2
    >>> fp,score = picker[0]
    >>> id = fp._id
    >>> str(id)
    'ether-1'
    >>> score
    1.0

    The results come back in order:
    >>> fp,score = picker[1]
    >>> id = fp._id
    >>> str(id)
    'ether-2'

    Now find the top matches for 2 probes.  We’ll get one ether and one acid:
    >>> fps = []
    >>> fps.append(Chem.RDKFingerprint(Chem.MolFromSmiles('COC')))
    >>> fps.append(Chem.RDKFingerprint(Chem.MolFromSmiles('CC(=O)O')))
    >>> picker = SpreadPicker(numToPick=3,probeFps=fps,dataSet=probefps)
    >>> len(picker)
    3
    >>> fp,score = picker[0]
    >>> id = fp._id
    >>> str(id)
    'ether-1'
    >>> score
    1.0
    >>> fp,score = picker[1]
    >>> id = fp._id
    >>> str(id)
    'acid-1'
    >>> score
    1.0
    >>> fp,score = picker[2]
    >>> id = fp._id
    >>> str(id)
    'ether-2'

    dataSet should be a sequence of BitVectors or, if expectPickles
    is False, a set of strings that can be converted to bit vectors"""

    numToPick: Incomplete
    probes: Incomplete
    data: Incomplete
    simMetric: Incomplete
    expectPickles: Incomplete
    onlyNames: Incomplete

    def __init__(
        self,
        numToPick: int = ...,
        probeFps: Incomplete | None = ...,
        dataSet: Incomplete | None = ...,
        simMetric=...,
        expectPickles: bool = ...,
        onlyNames: bool = ...,
    ) -> None: ...
    def MakePicks(self, force=False, silent=False): ...

class TopNOverallPicker(GenericPicker):
    """
    A class for picking the top N overall best matches across a library
    Connect to a database and build molecules:
    >>> from rdkit import Chem
    >>> from rdkit import RDConfig
    >>> import os.path
    >>> from rdkit.Dbase.DbConnection import DbConnect
    >>> dbName = RDConfig.RDTestDatabase
    >>> conn = DbConnect(dbName,'simple_mols1')
    >>> [x.upper() for x in conn.GetColumnNames()]
    ['SMILES', 'ID']
    >>> mols = []
    >>> for smi,id in conn.GetData():
    ...   mol = Chem.MolFromSmiles(str(smi))
    ...   mol.SetProp('_Name',str(id))
    ...   mols.append(mol)
    >>> len(mols)
    12

    Calculate fingerprints:
    >>> probefps = []
    >>> for mol in mols:
    ...   fp = Chem.RDKFingerprint(mol)
    ...   fp._id = mol.GetProp('_Name')
    ...   probefps.append(fp)

    Start by finding the top matches for a single probe.  This ether should pull
    other ethers from the db:
    >>> mol = Chem.MolFromSmiles('COC')
    >>> probeFp = Chem.RDKFingerprint(mol)
    >>> picker = TopNOverallPicker(numToPick=2,probeFps=[probeFp],dataSet=probefps)
    >>> len(picker)
    2
    >>> fp,score = picker[0]
    >>> id = fp._id
    >>> str(id)
    'ether-1'
    >>> score
    1.0

    The results come back in order:
    >>> fp,score = picker[1]
    >>> id = fp._id
    >>> str(id)
    'ether-2'

    Now find the top matches for 2 probes.  We’ll get one ether and one acid:
    >>> fps = []
    >>> fps.append(Chem.RDKFingerprint(Chem.MolFromSmiles('COC')))
    >>> fps.append(Chem.RDKFingerprint(Chem.MolFromSmiles('CC(=O)O')))
    >>> picker = TopNOverallPicker(numToPick=3,probeFps=fps,dataSet=probefps)
    >>> len(picker)
    3
    >>> fp,score = picker[0]
    >>> id = fp._id
    >>> str(id)
    'acid-1'
    >>> fp,score = picker[1]
    >>> id = fp._id
    >>> str(id)
    'ether-1'
    >>> score
    1.0
    >>> fp,score = picker[2]
    >>> id = fp._id
    >>> str(id)
    'acid-2'

    dataSet should be a sequence of BitVectors"""

    numToPick: Incomplete
    probes: Incomplete
    data: Incomplete
    simMetric: Incomplete

    def __init__(
        self,
        numToPick: int = ...,
        probeFps: Incomplete | None = ...,
        dataSet: Incomplete | None = ...,
        simMetric=...,
    ) -> None: ...
    def MakePicks(self, force=False): ...

class SpreadPicker(GenericPicker):
    """
    A class for picking the best matches across a library
    Connect to a database:
    >>> from rdkit import Chem
    >>> from rdkit import RDConfig
    >>> import os.path
    >>> from rdkit.Dbase.DbConnection import DbConnect
    >>> dbName = RDConfig.RDTestDatabase
    >>> conn = DbConnect(dbName,'simple_mols1')
    >>> [x.upper() for x in conn.GetColumnNames()]
    ['SMILES', 'ID']
    >>> mols = []
    >>> for smi,id in conn.GetData():
    ...   mol = Chem.MolFromSmiles(str(smi))
    ...   mol.SetProp('_Name',str(id))
    ...   mols.append(mol)
    >>> len(mols)
    12

    Calculate fingerprints:
    >>> probefps = []
    >>> for mol in mols:
    ...   fp = Chem.RDKFingerprint(mol)
    ...   fp._id = mol.GetProp('_Name')
    ...   probefps.append(fp)

    Start by finding the top matches for a single probe.  This ether should pull
    other ethers from the db:
    >>> mol = Chem.MolFromSmiles('COC')
    >>> probeFp = Chem.RDKFingerprint(mol)
    >>> picker = SpreadPicker(numToPick=2,probeFps=[probeFp],dataSet=probefps)
    >>> len(picker)
    2
    >>> fp,score = picker[0]
    >>> id = fp._id
    >>> str(id)
    'ether-1'
    >>> score
    1.0

    The results come back in order:
    >>> fp,score = picker[1]
    >>> id = fp._id
    >>> str(id)
    'ether-2'

    Now find the top matches for 2 probes.  We’ll get one ether and one acid:
    >>> fps = []
    >>> fps.append(Chem.RDKFingerprint(Chem.MolFromSmiles('COC')))
    >>> fps.append(Chem.RDKFingerprint(Chem.MolFromSmiles('CC(=O)O')))
    >>> picker = SpreadPicker(numToPick=3,probeFps=fps,dataSet=probefps)
    >>> len(picker)
    3
    >>> fp,score = picker[0]
    >>> id = fp._id
    >>> str(id)
    'ether-1'
    >>> score
    1.0
    >>> fp,score = picker[1]
    >>> id = fp._id
    >>> str(id)
    'acid-1'
    >>> score
    1.0
    >>> fp,score = picker[2]
    >>> id = fp._id
    >>> str(id)
    'ether-2'

    dataSet should be a sequence of BitVectors or, if expectPickles
    is False, a set of strings that can be converted to bit vectors"""

    numToPick: Incomplete
    probes: Incomplete
    data: Incomplete
    simMetric: Incomplete
    expectPickles: Incomplete
    onlyNames: Incomplete

    def __init__(
        self,
        numToPick: int = ...,
        probeFps: Incomplete | None = ...,
        dataSet: Incomplete | None = ...,
        simMetric=...,
        expectPickles: bool = ...,
        onlyNames: bool = ...,
    ) -> None: ...
    def MakePicks(self, force=False, silent=False): ...
