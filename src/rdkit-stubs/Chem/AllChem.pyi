"""
rdkit.Chem.AllChem module¶
Import all RDKit chemistry modules
"""
from collections.abc import Generator
from typing import Callable

from _typeshed import Incomplete
from rdkit import DataStructs as DataStructs
from rdkit import ForceField as ForceField
from rdkit import RDConfig as RDConfig
from rdkit import rdBase as rdBase
from rdkit.Chem import *
from rdkit.Chem.ChemicalFeatures import *
from rdkit.Chem.EnumerateStereoisomers import (
    EnumerateStereoisomers as EnumerateStereoisomers,
)
from rdkit.Chem.EnumerateStereoisomers import (
    StereoEnumerationOptions as StereoEnumerationOptions,
)
from rdkit.Chem.rdChemReactions import *
from rdkit.Chem.rdDepictor import *
from rdkit.Chem.rdDistGeom import *
from rdkit.Chem.rdFingerprintGenerator import *
from rdkit.Chem.rdForceFieldHelpers import *
from rdkit.Chem.rdMolAlign import *
from rdkit.Chem.rdMolDescriptors import *
from rdkit.Chem.rdMolEnumerator import *
from rdkit.Chem.rdMolTransforms import *
from rdkit.Chem.rdPartialCharges import *
from rdkit.Chem.rdqueries import *
from rdkit.Chem.rdReducedGraphs import *
from rdkit.Chem.rdShapeHelpers import *
from rdkit.Chem.rdSLNParse import *
from rdkit.Geometry import rdGeometry as rdGeometry
from rdkit.RDLogger import logger as logger

def AssignBondOrdersFromTemplate(self, refmol, mol):
    """
    assigns bond orders to a molecule based on the
    bond orders in a template molecule

    Arguments
    refmol: the template molecule
    mol: the molecule to assign bond orders to

    An example, start by generating a template from a SMILES
    and read in the PDB structure of the molecule
    >>> import os
    >>> from rdkit.Chem import AllChem
    >>> template = AllChem.MolFromSmiles("CN1C(=NC(C1=O)(c2ccccc2)c3ccccc3)N")
    >>> mol = AllChem.MolFromPDBFile(os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', '4DJU_lig.pdb'))
    >>> len([1 for b in template.GetBonds() if b.GetBondTypeAsDouble() == 1.0])
    8
    >>> len([1 for b in mol.GetBonds() if b.GetBondTypeAsDouble() == 1.0])
    22

    Now assign the bond orders based on the template molecule
    >>> newMol = AllChem.AssignBondOrdersFromTemplate(template, mol)
    >>> len([1 for b in newMol.GetBonds() if b.GetBondTypeAsDouble() == 1.0])
    8

    Note that the template molecule should have no explicit hydrogens
    else the algorithm will fail.
    It also works if there are different formal charges (this was github issue 235):
    >>> template=AllChem.MolFromSmiles('CN(C)C(=O)Cc1ccc2c(c1)NC(=O)c3ccc(cc3N2)c4ccc(c(c4)OC)[N+](=O)[O-]')
    >>> mol = AllChem.MolFromMolFile(os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', '4FTR_lig.mol'))
    >>> AllChem.MolToSmiles(mol)
    'COC1CC(C2CCC3C(O)NC4CC(CC(O)N(C)C)CCC4NC3C2)CCC1N(O)O'
    >>> newMol = AllChem.AssignBondOrdersFromTemplate(template, mol)
    >>> AllChem.MolToSmiles(newMol)
    'COc1cc(-c2ccc3c(c2)Nc2ccc(CC(=O)N(C)C)cc2NC3=O)ccc1[N+](=O)[O-]'"""
    ...

def TransformMol(self, mol, tform, confId=-1, keepConfs=False):
    """
    Applies the transformation (usually a 4x4 double matrix) to a molecule
    if keepConfs is False then all but that conformer are removed"""
    ...

def ComputeMolShape(self, mol, confId=-1, boxDim=(20, 20, 20), spacing=0.5, **kwargs):
    """
    returns a grid representation of the molecule’s shape"""
    ...

def ComputeMolVolume(self, mol, confId=-1, gridSpacing=0.2, boxMargin=2.0):
    """
    Calculates the volume of a particular conformer of a molecule
    based on a grid-encoding of the molecular shape.
    A bit of demo as well as a test of github #1883:
    >>> from rdkit import Chem
    >>> from rdkit.Chem import AllChem
    >>> mol = Chem.AddHs(Chem.MolFromSmiles('C'))
    >>> AllChem.EmbedMolecule(mol)
    0
    >>> ComputeMolVolume(mol)
    28...
    >>> mol = Chem.AddHs(Chem.MolFromSmiles('O'))
    >>> AllChem.EmbedMolecule(mol)
    0
    >>> ComputeMolVolume(mol)
    20..."""
    ...

def ConstrainedEmbed(
    self,
    mol,
    core,
    useTethers=True,
    coreConfId=-1,
    randomseed=2342,
    getForceField: Callable = ...,
    **kwargs,
):
    """
    generates an embedding of a molecule where part of the molecule
    is constrained to have particular coordinates

    Arguments
    mol: the molecule to embed
    core: the molecule to use as a source of constraints

    useTethers: (optional) if True, the final conformation will beoptimized subject to a series of extra forces that pull the
    matching atoms to the positions of the core atoms. Otherwise
    simple distance constraints based on the core atoms will be
    used in the optimization.

    coreConfId: (optional) id of the core conformation to use
    randomSeed: (optional) seed for the random number generator

    An example, start by generating a template with a 3D structure:
    >>> from rdkit.Chem import AllChem
    >>> template = AllChem.MolFromSmiles("c1nn(Cc2ccccc2)cc1")
    >>> AllChem.EmbedMolecule(template)
    0
    >>> AllChem.UFFOptimizeMolecule(template)
    0

    Here’s a molecule:
    >>> mol = AllChem.MolFromSmiles("c1nn(Cc2ccccc2)cc1-c3ccccc3")

    Now do the constrained embedding
    >>> mol = AllChem.ConstrainedEmbed(mol, template)
    Demonstrate that the positions are nearly the same with template:
    >>> import math
    >>> molp = mol.GetConformer().GetAtomPosition(0)
    >>> templatep = template.GetConformer().GetAtomPosition(0)
    >>> all(math.isclose(v, 0.0, abs_tol=0.01) for v in molp-templatep)
    True
    >>> molp = mol.GetConformer().GetAtomPosition(1)
    >>> templatep = template.GetConformer().GetAtomPosition(1)
    >>> all(math.isclose(v, 0.0, abs_tol=0.01) for v in molp-templatep)
    True"""
    ...

def EnumerateLibraryFromReaction(self, reaction, sidechainSets, returnReactants=False):
    """
    Returns a generator for the virtual library defined by
    a reaction and a sequence of sidechain sets
    >>> from rdkit import Chem
    >>> from rdkit.Chem import AllChem
    >>> s1=[Chem.MolFromSmiles(x) for x in ('NC','NCC')]
    >>> s2=[Chem.MolFromSmiles(x) for x in ('OC=O','OC(=O)C')]
    >>> rxn = AllChem.ReactionFromSmarts('[O:2]=[C:1][OH].[N:3]>>[O:2]=[C:1][N:3]')
    >>> r = AllChem.EnumerateLibraryFromReaction(rxn,[s2,s1])
    >>> [Chem.MolToSmiles(x[0]) for x in list(r)]
    ['CNC=O', 'CCNC=O', 'CNC(C)=O', 'CCNC(C)=O']

    Note that this is all done in a lazy manner, so “infinitely” large libraries can
    be done without worrying about running out of memory. Your patience will run out first:
    Define a set of 10000 amines:
    >>> amines = (Chem.MolFromSmiles('N'+'C'*x) for x in range(10000))

    … a set of 10000 acids
    >>> acids = (Chem.MolFromSmiles('OC(=O)'+'C'*x) for x in range(10000))

    … now the virtual library (1e8 compounds in principle):
    >>> r = AllChem.EnumerateLibraryFromReaction(rxn,[acids,amines])

    … look at the first 4 compounds:
    >>> [Chem.MolToSmiles(next(r)[0]) for x in range(4)]
    ['NC=O', 'CNC=O', 'CCNC=O', 'CCCNC=O']"""
    ...

def GetConformerRMS(self, mol, confId1, confId2, atomIds=None, prealigned=False):
    """
    Returns the RMS between two conformations.
    By default, the conformers will be aligned to the first conformer
    before the RMS calculation and, as a side-effect, the second will be left
    in the aligned state.

    Arguments:
    mol:        the molecule
    confId1:    the id of the first conformer
    confId2:    the id of the second conformer

    atomIds:    (optional) list of atom ids to use a points foralingment - defaults to all atoms

    prealigned: (optional) by default the conformers are assumedbe unaligned and the second conformer be aligned
    to the first"""
    ...

def GetConformerRMSMatrix(self, mol, atomIds=None, prealigned=False):
    """
    Returns the RMS matrix of the conformers of a molecule.
    As a side-effect, the conformers will be aligned to the first
    conformer (i.e. the reference) and will left in the aligned state.

    Arguments:
    mol:     the molecule

    atomIds: (optional) list of atom ids to use a points foralingment - defaults to all atoms

    prealigned: (optional) by default the conformers are assumedbe unaligned and will therefore be aligned to the
    first conformer

    Note that the returned RMS matrix is symmetrical, i.e. it is the
    lower half of the matrix, e.g. for 5 conformers:
    rmsmatrix = [ a,
                  b, c,
                  d, e, f,
                  g, h, i, j]

    where a is the RMS between conformers 0 and 1, b is the RMS between
    conformers 0 and 2, etc.
    This way it can be directly used as distance matrix in e.g. Butina
    clustering."""
    ...

def TransformMol(self, mol, tform, confId=-1, keepConfs=False):
    """
    Applies the transformation (usually a 4x4 double matrix) to a molecule
    if keepConfs is False then all but that conformer are removed"""
    ...

def EnumerateLibraryFromReaction(self, reaction, sidechainSets, returnReactants=False):
    """
    Returns a generator for the virtual library defined by
    a reaction and a sequence of sidechain sets
    >>> from rdkit import Chem
    >>> from rdkit.Chem import AllChem
    >>> s1=[Chem.MolFromSmiles(x) for x in ('NC','NCC')]
    >>> s2=[Chem.MolFromSmiles(x) for x in ('OC=O','OC(=O)C')]
    >>> rxn = AllChem.ReactionFromSmarts('[O:2]=[C:1][OH].[N:3]>>[O:2]=[C:1][N:3]')
    >>> r = AllChem.EnumerateLibraryFromReaction(rxn,[s2,s1])
    >>> [Chem.MolToSmiles(x[0]) for x in list(r)]
    ['CNC=O', 'CCNC=O', 'CNC(C)=O', 'CCNC(C)=O']

    Note that this is all done in a lazy manner, so “infinitely” large libraries can
    be done without worrying about running out of memory. Your patience will run out first:
    Define a set of 10000 amines:
    >>> amines = (Chem.MolFromSmiles('N'+'C'*x) for x in range(10000))

    … a set of 10000 acids
    >>> acids = (Chem.MolFromSmiles('OC(=O)'+'C'*x) for x in range(10000))

    … now the virtual library (1e8 compounds in principle):
    >>> r = AllChem.EnumerateLibraryFromReaction(rxn,[acids,amines])

    … look at the first 4 compounds:
    >>> [Chem.MolToSmiles(next(r)[0]) for x in range(4)]
    ['NC=O', 'CNC=O', 'CCNC=O', 'CCCNC=O']"""
    ...

def ConstrainedEmbed(
    self,
    mol,
    core,
    useTethers=True,
    coreConfId=-1,
    randomseed=2342,
    getForceField: Callable = ...,
    **kwargs,
):
    """
    generates an embedding of a molecule where part of the molecule
    is constrained to have particular coordinates

    Arguments
    mol: the molecule to embed
    core: the molecule to use as a source of constraints

    useTethers: (optional) if True, the final conformation will beoptimized subject to a series of extra forces that pull the
    matching atoms to the positions of the core atoms. Otherwise
    simple distance constraints based on the core atoms will be
    used in the optimization.

    coreConfId: (optional) id of the core conformation to use
    randomSeed: (optional) seed for the random number generator

    An example, start by generating a template with a 3D structure:
    >>> from rdkit.Chem import AllChem
    >>> template = AllChem.MolFromSmiles("c1nn(Cc2ccccc2)cc1")
    >>> AllChem.EmbedMolecule(template)
    0
    >>> AllChem.UFFOptimizeMolecule(template)
    0

    Here’s a molecule:
    >>> mol = AllChem.MolFromSmiles("c1nn(Cc2ccccc2)cc1-c3ccccc3")

    Now do the constrained embedding
    >>> mol = AllChem.ConstrainedEmbed(mol, template)
    Demonstrate that the positions are nearly the same with template:
    >>> import math
    >>> molp = mol.GetConformer().GetAtomPosition(0)
    >>> templatep = template.GetConformer().GetAtomPosition(0)
    >>> all(math.isclose(v, 0.0, abs_tol=0.01) for v in molp-templatep)
    True
    >>> molp = mol.GetConformer().GetAtomPosition(1)
    >>> templatep = template.GetConformer().GetAtomPosition(1)
    >>> all(math.isclose(v, 0.0, abs_tol=0.01) for v in molp-templatep)
    True"""
    ...

def AssignBondOrdersFromTemplate(self, refmol, mol):
    """
    assigns bond orders to a molecule based on the
    bond orders in a template molecule

    Arguments
    refmol: the template molecule
    mol: the molecule to assign bond orders to

    An example, start by generating a template from a SMILES
    and read in the PDB structure of the molecule
    >>> import os
    >>> from rdkit.Chem import AllChem
    >>> template = AllChem.MolFromSmiles("CN1C(=NC(C1=O)(c2ccccc2)c3ccccc3)N")
    >>> mol = AllChem.MolFromPDBFile(os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', '4DJU_lig.pdb'))
    >>> len([1 for b in template.GetBonds() if b.GetBondTypeAsDouble() == 1.0])
    8
    >>> len([1 for b in mol.GetBonds() if b.GetBondTypeAsDouble() == 1.0])
    22

    Now assign the bond orders based on the template molecule
    >>> newMol = AllChem.AssignBondOrdersFromTemplate(template, mol)
    >>> len([1 for b in newMol.GetBonds() if b.GetBondTypeAsDouble() == 1.0])
    8

    Note that the template molecule should have no explicit hydrogens
    else the algorithm will fail.
    It also works if there are different formal charges (this was github issue 235):
    >>> template=AllChem.MolFromSmiles('CN(C)C(=O)Cc1ccc2c(c1)NC(=O)c3ccc(cc3N2)c4ccc(c(c4)OC)[N+](=O)[O-]')
    >>> mol = AllChem.MolFromMolFile(os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', '4FTR_lig.mol'))
    >>> AllChem.MolToSmiles(mol)
    'COC1CC(C2CCC3C(O)NC4CC(CC(O)N(C)C)CCC4NC3C2)CCC1N(O)O'
    >>> newMol = AllChem.AssignBondOrdersFromTemplate(template, mol)
    >>> AllChem.MolToSmiles(newMol)
    'COc1cc(-c2ccc3c(c2)Nc2ccc(CC(=O)N(C)C)cc2NC3=O)ccc1[N+](=O)[O-]'"""
    ...
