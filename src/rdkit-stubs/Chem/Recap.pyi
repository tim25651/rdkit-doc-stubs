"""
rdkit.Chem.Recap module¶
Implementation of the RECAP algorithm from Lewell et al. JCICS 38 511-522 (1998)
The published algorithm is implemented more or less without
modification. The results are returned as a hierarchy of nodes instead
of just as a set of fragments. The hope is that this will allow a bit
more flexibility in working with the results.
For example:
>>> from rdkit import Chem
>>> from rdkit.Chem import Recap
>>> m = Chem.MolFromSmiles(‘C1CC1Oc1ccccc1-c1ncc(OC)cc1’)
>>> res = Recap.RecapDecompose(m)
>>> res
<…Chem.Recap.RecapHierarchyNode object at …>
>>> sorted(res.children.keys())
[’C1CC1’, ‘*c1ccc(OC)cn1’, ‘*c1ccccc1-c1ccc(OC)cn1’, ‘*c1ccccc1OC1CC1’]
>>> sorted(res.GetAllChildren().keys())
[‘*C1CC1’, ‘*c1ccc(OC)cn1’, ‘*c1ccccc1’, ‘*c1ccccc1-c1ccc(OC)cn1’, ‘*c1ccccc1OC1CC1’]
To get the standard set of RECAP results, use GetLeaves():
>>> leaves=res.GetLeaves()
>>> sorted(leaves.keys())
[’C1CC1’, ‘*c1ccc(OC)cn1’, ‘*c1ccccc1’]
>>> leaf = leaves[’*C1CC1’]
>>> leaf.mol
<…Chem.rdchem.Mol object at …>
"""
from _typeshed import Incomplete
from rdkit import Chem as Chem

reactionDefs: Incomplete
reactions: Incomplete

class RecapHierarchyNode(object):
    """
    This class is used to hold the Recap hiearchy"""

    children: None = ...
    mol: None = ...
    children: None = ...
    parents: None = ...
    smiles: None = ...

    def __init__(self, mol) -> None: ...
    def GetAllChildren(self):
        """
        returns a dictionary, keyed by SMILES, of children"""
        ...
    def GetLeaves(self):
        """
        returns a dictionary, keyed by SMILES, of leaf (terminal) nodes"""
        ...
    def getUltimateParents(self):
        """
        returns all the nodes in the hierarchy tree that contain this
        node as a child"""
        ...
    def __del__(self) -> None: ...

def RecapDecompose(self, mol, allNodes=None, minFragmentSize=0, onlyUseReactions=None):
    """
    returns the recap decomposition for a molecule"""
    ...
