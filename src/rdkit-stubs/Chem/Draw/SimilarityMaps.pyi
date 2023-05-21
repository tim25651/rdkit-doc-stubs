"""
rdkit.Chem.Draw.SimilarityMaps module¶
"""
from typing import Callable

from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit import DataStructs as DataStructs
from rdkit import Geometry as Geometry
from rdkit.Chem import Draw as Draw
from rdkit.Chem import rdDepictor as rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D as rdMolDraw2D

def GetAPFingerprint(
    self,
    mol,
    atomId=-1,
    fpType="normal",
    nBits=2048,
    minLength=1,
    maxLength=30,
    nBitsPerEntry=4,
    **kwargs,
):
    """
    Calculates the atom pairs fingerprint with the torsions of atomId removed.

    Parameters:mol – the molecule of interest
    atomId – the atom to remove the pairs for (if -1, no pair is removed)
    fpType – the type of AP fingerprint (‘normal’, ‘hashed’, ‘bv’)
    nBits – the size of the bit vector (only for fpType=’bv’)
    minLength – the minimum path length for an atom pair
    maxLength – the maxmimum path length for an atom pair
    nBitsPerEntry – the number of bits available for each pair"""
    ...

def GetAtomicWeightsForFingerprint(
    self, refMol, probeMol, fpFunction, metric: Callable = ...
):
    """
    Calculates the atomic weights for the probe molecule
    based on a fingerprint function and a metric.

    Parameters:refMol – the reference molecule
    probeMol – the probe molecule
    fpFunction – the fingerprint function
    metric – the similarity metric

    Note:If fpFunction needs additional parameters, use a lambda construct"""
    ...

def GetAtomicWeightsForModel(self, probeMol, fpFunction, predictionFunction):
    """
    Calculates the atomic weights for the probe molecule based on
    a fingerprint function and the prediction function of a ML model.

    Parameters:probeMol – the probe molecule
    fpFunction – the fingerprint function
    predictionFunction – the prediction function of the ML model"""
    ...

def GetMorganFingerprint(
    self, mol, atomId=-1, radius=2, fpType="bv", nBits=2048, useFeatures=False, **kwargs
):
    """
    Calculates the Morgan fingerprint with the environments of atomId removed.

    Parameters:mol – the molecule of interest
    radius – the maximum radius
    fpType – the type of Morgan fingerprint: ‘count’ or ‘bv’
    atomId – the atom to remove the environments for (if -1, no environments is removed)
    nBits – the size of the bit vector (only for fpType = ‘bv’)
    useFeatures – if false: ConnectivityMorgan, if true: FeatureMorgan

    any additional keyword arguments will be passed to the fingerprinting function."""
    ...

def GetRDKFingerprint(
    self,
    mol,
    atomId=-1,
    fpType="bv",
    nBits=2048,
    minPath=1,
    maxPath=5,
    nBitsPerHash=2,
    **kwargs,
):
    """
    Calculates the RDKit fingerprint with the paths of atomId removed.

    Parameters:mol – the molecule of interest
    atomId – the atom to remove the paths for (if -1, no path is removed)
    fpType – the type of RDKit fingerprint: ‘bv’
    nBits – the size of the bit vector
    minPath – minimum path length
    maxPath – maximum path length
    nBitsPerHash – number of to set per path"""
    ...

def GetStandardizedWeights(self, weights):
    """
    Normalizes the weights,
    such that the absolute maximum weight equals 1.0.

    Parameters:weights – the list with the atomic weights"""
    ...

def GetSimilarityMapFromWeights(
    self,
    mol,
    weights,
    colorMap=None,
    scale=-1,
    size=(250, 250),
    sigma=None,
    coordScale=1.5,
    step=0.01,
    colors="k",
    contourLines=10,
    alpha=0.5,
    draw2d=None,
    **kwargs,
):
    """
    Generates the similarity map for a molecule given the atomic weights.

    Parameters:mol – the molecule of interest
    colorMap – the matplotlib color map scheme, default is custom PiWG color map
    scale – the scaling: scale < 0 -> the absolute maximum weight is used as maximum scale

    scale = double -> this is the maximum scale

    size – the size of the figure
    sigma – the sigma for the Gaussians
    coordScale – scaling factor for the coordinates
    step – the step for calcAtomGaussian
    colors – color of the contour lines
    contourLines – if integer number N: N contour lines are drawn

    if list(numbers): contour lines at these numbers are drawn

    alpha – the alpha blending value for the contour lines
    kwargs – additional arguments for drawing"""
    ...

def GetSimilarityMapForFingerprint(
    self, refMol, probeMol, fpFunction, metric: Callable = ..., **kwargs
):
    """
    Generates the similarity map for a given reference and probe molecule,
    fingerprint function and similarity metric.

    Parameters:refMol – the reference molecule
    probeMol – the probe molecule
    fpFunction – the fingerprint function
    metric – the similarity metric.
    kwargs – additional arguments for drawing"""
    ...

def GetSimilarityMapForModel(self, probeMol, fpFunction, predictionFunction, **kwargs):
    """
    Generates the similarity map for a given ML model and probe molecule,
    and fingerprint function.

    Parameters:probeMol – the probe molecule
    fpFunction – the fingerprint function
    predictionFunction – the prediction function of the ML model
    kwargs – additional arguments for drawing"""
    ...

def GetSimilarityMapFromWeights(
    self,
    mol,
    weights,
    colorMap=None,
    scale=-1,
    size=(250, 250),
    sigma=None,
    coordScale=1.5,
    step=0.01,
    colors="k",
    contourLines=10,
    alpha=0.5,
    draw2d=None,
    **kwargs,
):
    """
    Generates the similarity map for a molecule given the atomic weights.

    Parameters:mol – the molecule of interest
    colorMap – the matplotlib color map scheme, default is custom PiWG color map
    scale – the scaling: scale < 0 -> the absolute maximum weight is used as maximum scale

    scale = double -> this is the maximum scale

    size – the size of the figure
    sigma – the sigma for the Gaussians
    coordScale – scaling factor for the coordinates
    step – the step for calcAtomGaussian
    colors – color of the contour lines
    contourLines – if integer number N: N contour lines are drawn

    if list(numbers): contour lines at these numbers are drawn

    alpha – the alpha blending value for the contour lines
    kwargs – additional arguments for drawing"""
    ...

def GetStandardizedWeights(self, weights):
    """
    Normalizes the weights,
    such that the absolute maximum weight equals 1.0.

    Parameters:weights – the list with the atomic weights"""
    ...

apDict: Incomplete

def GetAPFingerprint(
    self,
    mol,
    atomId=-1,
    fpType="normal",
    nBits=2048,
    minLength=1,
    maxLength=30,
    nBitsPerEntry=4,
    **kwargs,
):
    """
    Calculates the atom pairs fingerprint with the torsions of atomId removed.

    Parameters:mol – the molecule of interest
    atomId – the atom to remove the pairs for (if -1, no pair is removed)
    fpType – the type of AP fingerprint (‘normal’, ‘hashed’, ‘bv’)
    nBits – the size of the bit vector (only for fpType=’bv’)
    minLength – the minimum path length for an atom pair
    maxLength – the maxmimum path length for an atom pair
    nBitsPerEntry – the number of bits available for each pair"""
    ...

ttDict: Incomplete

def GetTTFingerprint(
    self,
    mol,
    atomId=-1,
    fpType="normal",
    nBits=2048,
    targetSize=4,
    nBitsPerEntry=4,
    **kwargs,
):
    """
    Calculates the topological torsion fingerprint with the pairs of atomId removed.

    Parameters:mol – the molecule of interest
    atomId – the atom to remove the torsions for (if -1, no torsion is removed)
    fpType – the type of TT fingerprint (‘normal’, ‘hashed’, ‘bv’)
    nBits – the size of the bit vector (only for fpType=’bv’)
    minLength – the minimum path length for an atom pair
    maxLength – the maxmimum path length for an atom pair
    nBitsPerEntry – the number of bits available for each torsion

    any additional keyword arguments will be passed to the fingerprinting function."""
    ...

def GetMorganFingerprint(
    self, mol, atomId=-1, radius=2, fpType="bv", nBits=2048, useFeatures=False, **kwargs
):
    """
    Calculates the Morgan fingerprint with the environments of atomId removed.

    Parameters:mol – the molecule of interest
    radius – the maximum radius
    fpType – the type of Morgan fingerprint: ‘count’ or ‘bv’
    atomId – the atom to remove the environments for (if -1, no environments is removed)
    nBits – the size of the bit vector (only for fpType = ‘bv’)
    useFeatures – if false: ConnectivityMorgan, if true: FeatureMorgan

    any additional keyword arguments will be passed to the fingerprinting function."""
    ...

def GetRDKFingerprint(
    self,
    mol,
    atomId=-1,
    fpType="bv",
    nBits=2048,
    minPath=1,
    maxPath=5,
    nBitsPerHash=2,
    **kwargs,
):
    """
    Calculates the RDKit fingerprint with the paths of atomId removed.

    Parameters:mol – the molecule of interest
    atomId – the atom to remove the paths for (if -1, no path is removed)
    fpType – the type of RDKit fingerprint: ‘bv’
    nBits – the size of the bit vector
    minPath – minimum path length
    maxPath – maximum path length
    nBitsPerHash – number of to set per path"""
    ...
