"""
rdkit.Chem.SaltRemover module¶
"""
from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit import RDConfig as RDConfig
from rdkit.Chem.rdmolfiles import SDMolSupplier as SDMolSupplier
from rdkit.Chem.rdmolfiles import SmilesMolSupplier as SmilesMolSupplier

class InputFormat(object):
    MOL: str = ...
    SMARTS: str = ...
    MOL: str = ...
    SMILES: str = ...

class SaltRemover(object):
    defnFilename: str = ...
    defnData: Incomplete
    salts: Incomplete
    defnFormat: Incomplete

    def __init__(
        self,
        defnFilename: Incomplete | None = ...,
        defnData: Incomplete | None = ...,
        defnFormat=...,
    ) -> None: ...
    def StripMol(self, mol, dontRemoveEverything=False, sanitize=True):
        """
        >>> remover = SaltRemover(defnData="[Cl,Br]")
        >>> len(remover.salts)
        1

        >>> mol = Chem.MolFromSmiles('CN(C)C.Cl')
        >>> res = remover.StripMol(mol)
        >>> res is not None
        True
        >>> res.GetNumAtoms()
        4

        Notice that all salts are removed:
        >>> mol = Chem.MolFromSmiles('CN(C)C.Cl.Cl.Br')
        >>> res = remover.StripMol(mol)
        >>> res.GetNumAtoms()
        4

        Matching (e.g. “salt-like”) atoms in the molecule are unchanged:
        >>> mol = Chem.MolFromSmiles('CN(Br)Cl')
        >>> res = remover.StripMol(mol)
        >>> res.GetNumAtoms()
        4

        >>> mol = Chem.MolFromSmiles('CN(Br)Cl.Cl')
        >>> res = remover.StripMol(mol)
        >>> res.GetNumAtoms()
        4

        Charged salts are handled reasonably:
        >>> mol = Chem.MolFromSmiles('C[NH+](C)(C).[Cl-]')
        >>> res = remover.StripMol(mol)
        >>> res.GetNumAtoms()
        4

        Watch out for this case (everything removed):
        >>> remover = SaltRemover()
        >>> len(remover.salts)>1
        True
        >>> mol = Chem.MolFromSmiles('CC(=O)O.[Na]')
        >>> res = remover.StripMol(mol)
        >>> res.GetNumAtoms()
        0

        dontRemoveEverything helps with this by leaving the last salt:
        >>> res = remover.StripMol(mol,dontRemoveEverything=True)
        >>> res.GetNumAtoms()
        4

        but in cases where the last salts are the same, it can’t choose
        between them, so it returns all of them:
        >>> mol = Chem.MolFromSmiles('Cl.Cl')
        >>> res = remover.StripMol(mol,dontRemoveEverything=True)
        >>> res.GetNumAtoms()
        2"""
        ...
    def StripMolWithDeleted(self, mol, dontRemoveEverything=False):
        """
        Strips given molecule and returns it, with the fragments which have been deleted.
        >>> remover = SaltRemover(defnData="[Cl,Br]")
        >>> len(remover.salts)
        1

        >>> mol = Chem.MolFromSmiles('CN(C)C.Cl.Br')
        >>> res, deleted = remover.StripMolWithDeleted(mol)
        >>> Chem.MolToSmiles(res)
        'CN(C)C'
        >>> [Chem.MolToSmarts(m) for m in deleted]
        ['[Cl,Br]']

        >>> mol = Chem.MolFromSmiles('CN(C)C.Cl')
        >>> res, deleted = remover.StripMolWithDeleted(mol)
        >>> res.GetNumAtoms()
        4
        >>> len(deleted)
        1
        >>> deleted[0].GetNumAtoms()
        1
        >>> Chem.MolToSmarts(deleted[0])
        '[Cl,Br]'

        Multiple occurrences of ‘Cl’ and without tuple destructuring
        >>> mol = Chem.MolFromSmiles('CN(C)C.Cl.Cl')
        >>> tup = remover.StripMolWithDeleted(mol)

        >>> tup.mol.GetNumAtoms()
        4
        >>> len(tup.deleted)
        1
        >>> tup.deleted[0].GetNumAtoms()
        1
        >>> Chem.MolToSmarts(deleted[0])
        '[Cl,Br]'"""
        ...
    def __call__(self, mol, dontRemoveEverything: bool = ...): ...
