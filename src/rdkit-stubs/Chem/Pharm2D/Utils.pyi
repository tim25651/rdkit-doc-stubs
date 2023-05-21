"""
rdkit.Chem.Pharm2D.Utils module¶
utility functionality for the 2D pharmacophores code
See Docs/Chem/Pharm2D.triangles.jpg for an illustration of the way
pharmacophores are broken into triangles and labelled.
See Docs/Chem/Pharm2D.signatures.jpg for an illustration of bit
numbering
"""
from _typeshed import Incomplete

nPointDistDict: Incomplete
nDistPointDict: Incomplete

def GetTriangles(self, nPts):
    """
    returns a tuple with the distance indices for
    triangles composing an nPts-pharmacophore"""
    ...

def BinsTriangleInequality(self, d1, d2, d3):
    """
    checks the triangle inequality for combinations
    of distance bins.

    the general triangle inequality is:d1 + d2 >= d3

    the conservative binned form of this is:d1(upper) + d2(upper) >= d3(lower)"""
    ...

def ScaffoldPasses(self, combo, bins=None):
    """
    checks the scaffold passed in to see if all
    contributing triangles can satisfy the triangle inequality
    the scaffold itself (encoded in combo) is a list of binned distances"""
    ...

def NumCombinations(self, nItems, nSlots):
    """
    returns the number of ways to fit nItems into nSlots
    We assume that (x, y) and (y, x) are equivalent, and
    (x, x) is allowed.

    General formula is, for N items and S slots:res = (N+S-1)! / ( (N-1)! * S! )"""
    ...

def CountUpTo(self, nItems, nSlots, vs, idx=0, startAt=0):
    """
    Figures out where a given combination of indices wouldoccur in the combinatorial explosion generated by _GetIndexCombinations_
    Arguments

    nItems: the number of items to distribute
    nSlots: the number of slots in which to distribute them
    vs: a sequence containing the values to find
    idx: used in the recursion
    startAt: used in the recursion

    Returns

    an integer"""
    ...

def GetAllCombinations(self, choices, noDups=1, which=0):
    """
    Does the combinatorial explosion of the possible combinations
    of the elements of _choices_.
    Arguments

    choices: sequence of sequences with the elements to be enumerated
    noDups: (optional) if this is nonzero, results with duplicates,
    e.g. (1,1,0), will not be generated
    which: used in recursion

    Returns

    a list of lists

    >>> GetAllCombinations([(0, ), (1, ), (2, )])
    [[0, 1, 2]]
    >>> GetAllCombinations([(0, ), (1, 3), (2, )])
    [[0, 1, 2], [0, 3, 2]]

    >>> GetAllCombinations([(0, 1), (1, 3), (2, )])
    [[0, 1, 2], [0, 3, 2], [1, 3, 2]]"""
    ...

def GetIndexCombinations(self, nItems, nSlots, slot=0, lastItemVal=0):
    """
    Generates all combinations of nItems in nSlots without includingduplicates

    Arguments

    nItems: the number of items to distribute
    nSlots: the number of slots in which to distribute them
    slot: used in recursion
    lastItemVal: used in recursion

    Returns

    a list of lists"""
    ...

def GetPossibleScaffolds(self, nPts, bins, useTriangleInequality=True):
    """
    gets all realizable scaffolds (passing the triangle inequality) with the
    given number of points and returns them as a list of tuples"""
    ...

def GetTriangles(self, nPts):
    """
    returns a tuple with the distance indices for
    triangles composing an nPts-pharmacophore"""
    ...

def GetAllCombinations(self, choices, noDups=1, which=0):
    """
    Does the combinatorial explosion of the possible combinations
    of the elements of _choices_.
    Arguments

    choices: sequence of sequences with the elements to be enumerated
    noDups: (optional) if this is nonzero, results with duplicates,
    e.g. (1,1,0), will not be generated
    which: used in recursion

    Returns

    a list of lists

    >>> GetAllCombinations([(0, ), (1, ), (2, )])
    [[0, 1, 2]]
    >>> GetAllCombinations([(0, ), (1, 3), (2, )])
    [[0, 1, 2], [0, 3, 2]]

    >>> GetAllCombinations([(0, 1), (1, 3), (2, )])
    [[0, 1, 2], [0, 3, 2], [1, 3, 2]]"""
    ...

def GetUniqueCombinations(self, choices, classes, which=0):
    """
    Does the combinatorial explosion of the possible combinations
    of the elements of _choices_."""
    ...

def GetUniqueCombinations_new(self, choices, classes, which=0):
    """
    Does the combinatorial explosion of the possible combinations
    of the elements of _choices_."""
    ...

def NumCombinations(self, nItems, nSlots):
    """
    returns the number of ways to fit nItems into nSlots
    We assume that (x, y) and (y, x) are equivalent, and
    (x, x) is allowed.

    General formula is, for N items and S slots:res = (N+S-1)! / ( (N-1)! * S! )"""
    ...

def OrderTriangle(self, featIndices, dists):
    """
    put the distances for a triangle into canonical order
    It’s easy if the features are all different:
    >>> OrderTriangle([0, 2, 4], [1, 2, 3])
    ([0, 2, 4], [1, 2, 3])

    It’s trickiest if they are all the same:
    >>> OrderTriangle([0, 0, 0], [1, 2, 3])
    ([0, 0, 0], [3, 2, 1])
    >>> OrderTriangle([0, 0, 0], [2, 1, 3])
    ([0, 0, 0], [3, 2, 1])
    >>> OrderTriangle([0, 0, 0], [1, 3, 2])
    ([0, 0, 0], [3, 2, 1])
    >>> OrderTriangle([0, 0, 0], [3, 1, 2])
    ([0, 0, 0], [3, 2, 1])
    >>> OrderTriangle([0, 0, 0], [3, 2, 1])
    ([0, 0, 0], [3, 2, 1])

    >>> OrderTriangle([0, 0, 1], [3, 2, 1])
    ([0, 0, 1], [3, 2, 1])
    >>> OrderTriangle([0, 0, 1], [1, 3, 2])
    ([0, 0, 1], [1, 3, 2])
    >>> OrderTriangle([0, 0, 1], [1, 2, 3])
    ([0, 0, 1], [1, 3, 2])"""
    ...

def ScaffoldPasses(self, combo, bins=None):
    """
    checks the scaffold passed in to see if all
    contributing triangles can satisfy the triangle inequality
    the scaffold itself (encoded in combo) is a list of binned distances"""
    ...

def UniquifyCombinations(self, combos):
    """
    uniquifies the combinations in the argument
    Arguments:

    combos: a sequence of sequences

    Returns

    a list of tuples containing the unique combos"""
    ...

def GetPossibleScaffolds(self, nPts, bins, useTriangleInequality=True):
    """
    gets all realizable scaffolds (passing the triangle inequality) with the
    given number of points and returns them as a list of tuples"""
    ...

def OrderTriangle(self, featIndices, dists):
    """
    put the distances for a triangle into canonical order
    It’s easy if the features are all different:
    >>> OrderTriangle([0, 2, 4], [1, 2, 3])
    ([0, 2, 4], [1, 2, 3])

    It’s trickiest if they are all the same:
    >>> OrderTriangle([0, 0, 0], [1, 2, 3])
    ([0, 0, 0], [3, 2, 1])
    >>> OrderTriangle([0, 0, 0], [2, 1, 3])
    ([0, 0, 0], [3, 2, 1])
    >>> OrderTriangle([0, 0, 0], [1, 3, 2])
    ([0, 0, 0], [3, 2, 1])
    >>> OrderTriangle([0, 0, 0], [3, 1, 2])
    ([0, 0, 0], [3, 2, 1])
    >>> OrderTriangle([0, 0, 0], [3, 2, 1])
    ([0, 0, 0], [3, 2, 1])

    >>> OrderTriangle([0, 0, 1], [3, 2, 1])
    ([0, 0, 1], [3, 2, 1])
    >>> OrderTriangle([0, 0, 1], [1, 3, 2])
    ([0, 0, 1], [1, 3, 2])
    >>> OrderTriangle([0, 0, 1], [1, 2, 3])
    ([0, 0, 1], [1, 3, 2])"""
    ...