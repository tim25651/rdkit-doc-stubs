"""
rdkit.Chem.rdMolAlign module¶
Module containing functions to align a molecule to a second molecule
"""
from typing import Any

import Boost.Python
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdMolDescriptors import AtomPairsParameters

class O3A(Boost.Python.instance):
    """
    Open3DALIGN object
    Raises an exception
    This class cannot be instantiated from Python"""

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def Align(self: O3A) -> float:
        """
        aligns probe molecule onto reference molecule

        C++ signature :double Align(RDKit::MolAlign::PyO3A {lvalue})"""
        ...
    def Matches(self: O3A) -> list:
        """
        returns the AtomMap as found by Open3DALIGN

        C++ signature :boost::python::list Matches(RDKit::MolAlign::PyO3A {lvalue})"""
        ...
    def Score(self: O3A) -> float:
        """
        returns the O3AScore of the alignment

        C++ signature :double Score(RDKit::MolAlign::PyO3A {lvalue})"""
        ...
    def Trans(self: O3A) -> object:
        """
        returns the transformation which aligns probe molecule onto reference molecule

        C++ signature :_object* Trans(RDKit::MolAlign::PyO3A {lvalue})"""
        ...
    def Weights(self: O3A) -> list:
        """
        returns the weight vector as found by Open3DALIGN

        C++ signature :boost::python::list Weights(RDKit::MolAlign::PyO3A {lvalue})"""
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

def AlignMol(
    self,
    prbMol: Mol,
    refMol: Mol,
    prbCid: int = -1,
    refCid: int = -1,
    atomMap: AtomPairsParameters = [],
    weights: AtomPairsParameters = [],
    reflect: bool = False,
    maxIters: int = 50,
) -> float:
    """
    Optimally (minimum RMSD) align a molecule to another molecule

    The 3D transformation required to align the specied conformation in the probe molecule
    to a specified conformation in the reference molecule is computed so that the root mean
    squared distance between a specified set of atoms is minimized.
    This transform is then applied to the specified conformation in the probe molecule

    ARGUMENTS
    prbMol    molecule that is to be aligned
    refMol    molecule used as the reference for the alignment

    prbCid    ID of the conformation in the probe to be used for the alignment (defaults to first conformation)

    refCid    ID of the conformation in the ref molecule to which the alignment is computed (defaults to first conformation)

    atomMap   a vector of pairs of atom IDs (probe AtomId, ref AtomId)used to compute the alignments. If this mapping is
    not specified an attempt is made to generate one by
    substructure matching

    weights   Optionally specify weights for each of the atom pairs
    reflect   if true reflect the conformation of the probe molecule
    maxIters  maximum number of iterations used in minimizing the RMSD

    RETURNS
    RMSD value

    C++ signature :double AlignMol(RDKit::ROMol {lvalue},RDKit::ROMol [,int=-1 [,int=-1 [,boost::python::api::object=[] [,boost::python::api::object=[] [,bool=False [,unsigned int=50]]]]]])
    """
    ...

def AlignMolConformers(
    self,
    mol: Mol,
    atomIds: AtomPairsParameters = [],
    confIds: AtomPairsParameters = [],
    weights: AtomPairsParameters = [],
    reflect: bool = False,
    maxIters: int = 50,
    RMSlist: AtomPairsParameters = None,
) -> None:
    """
    Align conformations in a molecule to each other

    The first conformation in the molecule is used as the reference

    ARGUMENTS
    mol          molecule of interest
    atomIds      List of atom ids to use a points for alingment - defaults to all atoms
    confIds      Ids of conformations to align - defaults to all conformers
    weights      Optionally specify weights for each of the atom pairs
    reflect      if true reflect the conformation of the probe molecule
    maxIters     maximum number of iterations used in minimizing the RMSD

    RMSlist      if provided, fills in the RMS values between the referenceconformation and the other aligned conformations

    C++ signature :void AlignMolConformers(RDKit::ROMol {lvalue} [,boost::python::api::object=[] [,boost::python::api::object=[] [,boost::python::api::object=[] [,bool=False [,unsigned int=50 [,boost::python::api::object=None]]]]]])
    """
    ...

def CalcRMS(
    self,
    prbMol: Mol,
    refMol: Mol,
    prbId: int = -1,
    refId: int = -1,
    map: AtomPairsParameters = None,
    maxMatches: int = 1000000,
    symmetrizeConjugatedTerminalGroups: bool = True,
    weights: AtomPairsParameters = [],
) -> float:
    """
    Returns the RMS between two molecules, taking symmetry into account.
    In contrast to getBestRMS, the RMS is computed ‘in place’, i.e.
    probe molecules are not aligned to the reference ahead of the
    RMS calculation. This is useful, for example, to compute
    the RMSD between docking poses and the co-crystallized ligand.
    Note:
    This function will attempt to match all permutations of matching atom
    orders in both molecules, for some molecules it will lead to
    ‘combinatorial explosion’ especially if hydrogens are present.

    ARGUMENTS
    prbMol:      the molecule to be aligned to the reference
    refMol:      the reference molecule
    prbCId:      (optional) probe conformation to use
    refCId:      (optional) reference conformation to use

    map:         (optional) a list of lists of (probeAtomId, refAtomId)tuples with the atom-atom mappings of the two
    molecules. If not provided, these will be generated
    using a substructure search.

    maxMatches:  (optional) if map isn’t specified, this will bethe max number of matches found in a SubstructMatch()

    symmetrizeConjugatedTerminalGroups:  (optional) if set, conjugatedterminal functional groups (like nitro or carboxylate)
    will be considered symmetrically

    weights:     (optional) weights for mapping

    RETURNS
    The best RMSD found

    C++ signature :double CalcRMS(RDKit::ROMol {lvalue},RDKit::ROMol {lvalue} [,int=-1 [,int=-1 [,boost::python::api::object=None [,int=1000000 [,bool=True [,boost::python::api::object=[]]]]]]])
    """
    ...

def GetAlignmentTransform(
    self,
    prbMol: Mol,
    refMol: Mol,
    prbCid: int = -1,
    refCid: int = -1,
    atomMap: AtomPairsParameters = [],
    weights: AtomPairsParameters = [],
    reflect: bool = False,
    maxIters: int = 50,
) -> object:
    """
    Compute the transformation required to align a molecule

    The 3D transformation required to align the specied conformation in the probe molecule
    to a specified conformation in the reference molecule is computed so that the root mean
    squared distance between a specified set of atoms is minimized

    ARGUMENTS
    prbMol    molecule that is to be aligned
    refMol    molecule used as the reference for the alignment

    prbCid    ID of the conformation in the probe to be used for the alignment (defaults to first conformation)

    refCid    ID of the conformation in the ref molecule to which the alignment is computed (defaults to first conformation)

    atomMap   a vector of pairs of atom IDs (probe AtomId, ref AtomId)used to compute the alignments. If this mapping is
    not specified an attempt is made to generate one by
    substructure matching

    weights   Optionally specify weights for each of the atom pairs
    reflect   if true reflect the conformation of the probe molecule
    maxIters  maximum number of iterations used in minimizing the RMSD

    RETURNS
    a tuple of (RMSD value, transform matrix)

    C++ signature :_object* GetAlignmentTransform(RDKit::ROMol,RDKit::ROMol [,int=-1 [,int=-1 [,boost::python::api::object=[] [,boost::python::api::object=[] [,bool=False [,unsigned int=50]]]]]])
    """
    ...

def GetBestAlignmentTransform(
    self,
    prbMol: Mol,
    refMol: Mol,
    prbCid: int = -1,
    refCid: int = -1,
    map: AtomPairsParameters = [],
    maxMatches: int = 1000000,
    symmetrizeConjugatedTerminalGroups: bool = True,
    weights: AtomPairsParameters = [],
    reflect: bool = False,
    maxIters: int = 50,
) -> object:
    """
    Compute the optimal RMS, transformation and atom map for aligning
    two molecules, taking symmetry into account. Molecule coordinates
    are left unaltered.
    This function will attempt to align all permutations of matching atom
    orders in both molecules, for some molecules it will lead to ‘combinatorial
    explosion’ especially if hydrogens are present.
    Use ‘GetAlignmentTransform’ to align molecules without changing the atom order.

    ARGUMENTS
    prbMol      molecule that is to be aligned
    refMol      molecule used as the reference for the alignment

    prbCid      ID of the conformation in the probe to be used for the alignment (defaults to first conformation)

    refCid      ID of the conformation in the ref molecule to which the alignment is computed (defaults to first conformation)

    map:        (optional) a list of lists of (probeAtomId, refAtomId)tuples with the atom-atom mappings of the two
    molecules. If not provided, these will be generated
    using a substructure search.

    maxMatches  (optional) if atomMap is empty, this will be the max number ofmatches found in a SubstructMatch().

    symmetrizeConjugatedTerminalGroups (optional) if set, conjugatedterminal functional groups (like nitro or carboxylate)
    will be considered symmetrically.

    weights     Optionally specify weights for each of the atom pairs
    reflect     if true reflect the conformation of the probe molecule
    maxIters    maximum number of iterations used in minimizing the RMSD

    RETURNS
    a tuple of (RMSD value, best transform matrix, best atom map)

    C++ signature :_object* GetBestAlignmentTransform(RDKit::ROMol,RDKit::ROMol [,int=-1 [,int=-1 [,boost::python::api::object=[] [,int=1000000 [,bool=True [,boost::python::api::object=[] [,bool=False [,unsigned int=50]]]]]]]])
    """
    ...

def GetBestRMS(
    self,
    prbMol: Mol,
    refMol: Mol,
    prbId: int = -1,
    refId: int = -1,
    map: AtomPairsParameters = None,
    maxMatches: int = 1000000,
    symmetrizeConjugatedTerminalGroups: bool = True,
    weights: AtomPairsParameters = [],
) -> float:
    """
    Returns the optimal RMS for aligning two molecules, taking
    symmetry into account. As a side-effect, the probe molecule is
    left in the aligned state.
    Note:
    This function will attempt to align all permutations of matching atom
    orders in both molecules, for some molecules it will lead to
    ‘combinatorial explosion’ especially if hydrogens are present.
    Use ‘rdkit.Chem.AllChem.AlignMol’ to align molecules without changing
    the atom order.

    ARGUMENTS
    prbMol:      the molecule to be aligned to the reference
    refMol:      the reference molecule
    prbId:       (optional) probe conformation to use
    refId:       (optional) reference conformation to use

    map:         (optional) a list of lists of (probeAtomId,refAtomId)tuples with the atom-atom mappings of the two
    molecules. If not provided, these will be generated
    using a substructure search.

    maxMatches:  (optional) if map isn’t specified, this will bethe max number of matches found in a SubstructMatch()

    symmetrizeConjugatedTerminalGroups:  (optional) if set, conjugatedterminal functional groups (like nitro or carboxylate)
    will be considered symmetrically

    weights:     (optional) weights for mapping

    RETURNS
    The best RMSD found

    C++ signature :double GetBestRMS(RDKit::ROMol {lvalue},RDKit::ROMol {lvalue} [,int=-1 [,int=-1 [,boost::python::api::object=None [,int=1000000 [,bool=True [,boost::python::api::object=[]]]]]]])
    """
    ...

def GetCrippenO3A(
    self,
    prbMol: Mol,
    refMol: Mol,
    prbCrippenContribs: list = [],
    refCrippenContribs: list = [],
    prbCid: int = -1,
    refCid: int = -1,
    reflect: bool = False,
    maxIters: int = 50,
    options: int = 0,
    constraintMap: list = [],
    constraintWeights: list = [],
) -> O3A:
    """
    Get an O3A object with atomMap and weights vectors to overlay
    the probe molecule onto the reference molecule based on
    Crippen logP atom contributions

    ARGUMENTS
    prbMol                   molecule that is to be aligned
    refMol                   molecule used as the reference for the alignment

    prbCrippenContribs       Crippen atom contributions for the probe moleculeas a list of (logp, mr) tuples, as returned
    by _CalcCrippenContribs()

    refCrippenContribs       Crippen atom contributions for the reference moleculeas a list of (logp, mr) tuples, as returned
    by _CalcCrippenContribs()

    prbCid                   ID of the conformation in the probe to be used for the alignment (defaults to first conformation)

    refCid                   ID of the conformation in the ref molecule to which the alignment is computed (defaults to first conformation)

    reflect                  if true reflect the conformation of the probe molecule(defaults to false)

    maxIters                 maximum number of iterations used in minimizing the RMSD(defaults to 50)

    options                  least 2 significant bits encode accuracy(0: maximum, 3: minimum; defaults to 0)
    bit 3 triggers local optimization of the alignment
    (no computation of the cost matrix; defaults: off)

    constraintMap            a vector of pairs of atom IDs (probe AtomId, ref AtomId)which shall be used for the alignment (defaults to [])

    constraintWeights        optionally specify weights for each of the constraints(weights default to 100.0)

    RETURNS
    The O3A object

    C++ signature :RDKit::MolAlign::PyO3A* GetCrippenO3A(RDKit::ROMol {lvalue},RDKit::ROMol {lvalue} [,boost::python::list=[] [,boost::python::list=[] [,int=-1 [,int=-1 [,bool=False [,unsigned int=50 [,unsigned int=0 [,boost::python::list=[] [,boost::python::list=[]]]]]]]]]])
    """
    ...

def GetCrippenO3AForProbeConfs(
    self,
    prbMol: Mol,
    refMol: Mol,
    numThreads: int = 1,
    prbCrippenContribs: list = [],
    refCrippenContribs: list = [],
    refCid: int = -1,
    reflect: bool = False,
    maxIters: int = 50,
    options: int = 0,
    constraintMap: list = [],
    constraintWeights: list = [],
) -> tuple:
    """
    Get a vector of O3A objects for the overlay of all
    the probe molecule’s conformations onto the reference molecule based on
    MMFF atom types and charges

    ARGUMENTS
    prbMol                   molecule that is to be aligned
    refMol                   molecule used as the reference for the alignment

    numThreadsthe number of threads to use, only has an effect ifthe RDKit was built with thread support (defaults to 1)

    prbCrippenContribs       Crippen atom contributions for the probe moleculeas a list of (logp, mr) tuples, as returned
    by _CalcCrippenContribs()

    refCrippenContribs       Crippen atom contributions for the reference moleculeas a list of (logp, mr) tuples, as returned
    by _CalcCrippenContribs()

    refCid                   ID of the conformation in the ref molecule to which the alignment is computed (defaults to first conformation)

    reflect                  if true reflect the conformation of the probe molecule(defaults to false)

    maxIters                 maximum number of iterations used in minimizing the RMSD(defaults to 50)

    options                  least 2 significant bits encode accuracy(0: maximum, 3: minimum; defaults to 0)
    bit 3 triggers local optimization of the alignment
    (no computation of the cost matrix; defaults: off)

    constraintMap            a vector of pairs of atom IDs (probe AtomId, ref AtomId)which shall be used for the alignment (defaults to [])

    constraintWeights        optionally specify weights for each of the constraints(weights default to 100.0)

    RETURNS
    A vector of O3A objects

    C++ signature :boost::python::tuple GetCrippenO3AForProbeConfs(RDKit::ROMol {lvalue},RDKit::ROMol {lvalue} [,int=1 [,boost::python::list=[] [,boost::python::list=[] [,int=-1 [,bool=False [,unsigned int=50 [,unsigned int=0 [,boost::python::list=[] [,boost::python::list=[]]]]]]]]]])
    """
    ...

def GetO3A(
    self,
    prbMol: Mol,
    refMol: Mol,
    prbPyMMFFMolProperties: AtomPairsParameters = None,
    refPyMMFFMolProperties: AtomPairsParameters = None,
    prbCid: int = -1,
    refCid: int = -1,
    reflect: bool = False,
    maxIters: int = 50,
    options: int = 0,
    constraintMap: list = [],
    constraintWeights: list = [],
) -> O3A:
    """
    Get an O3A object with atomMap and weights vectors to overlay
    the probe molecule onto the reference molecule based on
    MMFF atom types and charges

    ARGUMENTS
    prbMol                   molecule that is to be aligned
    refMol                   molecule used as the reference for the alignment

    prbPyMMFFMolProperties   PyMMFFMolProperties object for the probe molecule as returnedby SetupMMFFForceField()

    refPyMMFFMolProperties   PyMMFFMolProperties object for the reference molecule as returnedby SetupMMFFForceField()

    prbCid                   ID of the conformation in the probe to be used for the alignment (defaults to first conformation)

    refCid                   ID of the conformation in the ref molecule to which the alignment is computed (defaults to first conformation)

    reflect                  if true reflect the conformation of the probe molecule(defaults to false)

    maxIters                 maximum number of iterations used in minimizing the RMSD(defaults to 50)

    options                  least 2 significant bits encode accuracy(0: maximum, 3: minimum; defaults to 0)
    bit 3 triggers local optimization of the alignment
    (no computation of the cost matrix; defaults: off)

    constraintMap            a vector of pairs of atom IDs (probe AtomId, ref AtomId)which shall be used for the alignment (defaults to [])

    constraintWeights        optionally specify weights for each of the constraints(weights default to 100.0)

    RETURNS
    The O3A object

    C++ signature :RDKit::MolAlign::PyO3A* GetO3A(RDKit::ROMol {lvalue},RDKit::ROMol {lvalue} [,boost::python::api::object=None [,boost::python::api::object=None [,int=-1 [,int=-1 [,bool=False [,unsigned int=50 [,unsigned int=0 [,boost::python::list=[] [,boost::python::list=[]]]]]]]]]])
    """
    ...

def GetO3AForProbeConfs(
    self,
    prbMol: Mol,
    refMol: Mol,
    numThreads: int = 1,
    prbPyMMFFMolProperties: AtomPairsParameters = None,
    refPyMMFFMolProperties: AtomPairsParameters = None,
    refCid: int = -1,
    reflect: bool = False,
    maxIters: int = 50,
    options: int = 0,
    constraintMap: list = [],
    constraintWeights: list = [],
) -> tuple:
    """
    Get a vector of O3A objects for the overlay of all
    the probe molecule’s conformations onto the reference molecule based on
    MMFF atom types and charges

    ARGUMENTS
    prbMol                   molecule that is to be aligned
    refMol                   molecule used as the reference for the alignment

    numThreadsthe number of threads to use, only has an effect ifthe RDKit was built with thread support (defaults to 1)
    If set to zero, the max supported by the system will be used.

    prbPyMMFFMolProperties   PyMMFFMolProperties object for the probe molecule as returnedby SetupMMFFForceField()

    refPyMMFFMolProperties   PyMMFFMolProperties object for the reference molecule as returnedby SetupMMFFForceField()

    refCid                   ID of the conformation in the ref molecule to which the alignment is computed (defaults to first conformation)

    reflect                  if true reflect the conformation of the probe molecule(defaults to false)

    maxIters                 maximum number of iterations used in minimizing the RMSD(defaults to 50)

    options                  least 2 significant bits encode accuracy(0: maximum, 3: minimum; defaults to 0)
    bit 3 triggers local optimization of the alignment
    (no computation of the cost matrix; defaults: off)

    constraintMap            a vector of pairs of atom IDs (probe AtomId, ref AtomId)which shall be used for the alignment (defaults to [])

    constraintWeights        optionally specify weights for each of the constraints(weights default to 100.0)

    RETURNS
    A vector of O3A objects

    C++ signature :boost::python::tuple GetO3AForProbeConfs(RDKit::ROMol {lvalue},RDKit::ROMol {lvalue} [,int=1 [,boost::python::api::object=None [,boost::python::api::object=None [,int=-1 [,bool=False [,unsigned int=50 [,unsigned int=0 [,boost::python::list=[] [,boost::python::list=[]]]]]]]]]])
    """
    ...

def RandomTransform(self, mol: Mol, cid: int = -1, seed: int = -1) -> None:
    """
    Perform a random transformation on a molecule

    ARGUMENTS
    mol    molecule that is to be transformed

    cid    ID of the conformation in the mol to be transformed(defaults to first conformation)

    seed   seed used to initialize the random generator(defaults to -1, that is no seeding)

    C++ signature :void RandomTransform(RDKit::ROMol {lvalue} [,int=-1 [,int=-1]])"""
    ...
