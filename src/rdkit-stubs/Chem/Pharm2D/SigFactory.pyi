"""
rdkit.Chem.Pharm2D.SigFactory module¶
contains factory class for producing signatures
"""
from _typeshed import Incomplete
from rdkit.Chem.Pharm2D import Utils as Utils
from rdkit.Chem.Pharmacophores import cUtils as cUtils
from rdkit.DataStructs import IntSparseIntVect as IntSparseIntVect
from rdkit.DataStructs import LongSparseIntVect as LongSparseIntVect
from rdkit.DataStructs import SparseBitVect as SparseBitVect

class SigFactory(object):
    """
    SigFactory’s are used by creating one, setting the relevant
    parameters, then calling the GetSignature() method each time a
    signature is required."""

    featFactory: Incomplete
    useCounts: Incomplete
    minPointCount: Incomplete
    maxPointCount: Incomplete
    shortestPathsOnly: Incomplete
    includeBondOrder: Incomplete
    trianglePruneBins: Incomplete
    skipFeats: Incomplete
    sigKlass: Incomplete

    def __init__(
        self,
        featFactory,
        useCounts: bool = ...,
        minPointCount: int = ...,
        maxPointCount: int = ...,
        shortestPathsOnly: bool = ...,
        includeBondOrder: bool = ...,
        skipFeats: Incomplete | None = ...,
        trianglePruneBins: bool = ...,
    ) -> None: ...
    def SetBins(self, bins):
        """
        bins should be a list of 2-tuples"""
        ...
    def GetBins(self): ...
    def GetBitDescription(self, bitIdx):
        """
        returns a text description of the bit
        Arguments

        bitIdx: an integer bit index

        Returns

        a string"""
        ...
    def GetNumBins(self): ...
    def GetSignature(self): ...
    def GetBitDescriptionAsText(self, bitIdx, includeBins=0, fullPage=1):
        """
        returns text with a description of the bit
        Arguments

        bitIdx: an integer bit index
        includeBins: (optional) if nonzero, information about the bins will be
        included as well
        fullPage: (optional) if nonzero, html headers and footers will
        be included (so as to make the output a complete page)

        Returns

        a string with the HTML"""
        ...
    def GetBitIdx(self, featIndices, dists, sortIndices=True):
        """
        returns the index for a pharmacophore described using a set offeature indices and distances

        Arguments*

        featIndices: a sequence of feature indices
        dists: a sequence of distance between the features, only the
        unique distances should be included, and they should be in the
        order defined in Utils.
        sortIndices : sort the indices

        Returns

        the integer bit index"""
        ...
    def GetBitInfo(self, idx):
        """
        returns information about the given bit
        Arguments

        idx: the bit index to be considered

        Returns

        a 3-tuple:

        the number of points in the pharmacophore
        the proto-pharmacophore (tuple of pattern indices)
        the scaffold (tuple of distance indices)"""
        ...
    def GetBitDescription(self, bitIdx):
        """
        returns a text description of the bit
        Arguments

        bitIdx: an integer bit index

        Returns

        a string"""
        ...
    def GetFeatFamilies(self): ...
    def GetMolFeats(self, mol): ...
    def GetNumBins(self): ...
    def GetSigSize(self): ...
    def GetSignature(self): ...
    def GetBitIdx(self, featIndices, dists, sortIndices=True):
        """
        returns the index for a pharmacophore described using a set offeature indices and distances

        Arguments*

        featIndices: a sequence of feature indices
        dists: a sequence of distance between the features, only the
        unique distances should be included, and they should be in the
        order defined in Utils.
        sortIndices : sort the indices

        Returns

        the integer bit index"""
        ...
    def GetBitInfo(self, idx):
        """
        returns information about the given bit
        Arguments

        idx: the bit index to be considered

        Returns

        a 3-tuple:

        the number of points in the pharmacophore
        the proto-pharmacophore (tuple of pattern indices)
        the scaffold (tuple of distance indices)"""
        ...
    def Init(self):
        """
        Initializes internal parameters.  This must be called after
        making any changes to the signature parameters"""
        ...
    def SetBins(self, bins):
        """
        bins should be a list of 2-tuples"""
        ...
    def GetSigSize(self): ...
