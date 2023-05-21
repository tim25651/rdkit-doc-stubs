"""
rdkit.Numerics.rdAlignment module¶
Module containing functions to align pairs of points in 3D
"""
from typing import Any

from rdkit.Chem.rdMolDescriptors import AtomPairsParameters

def GetAlignmentTransform(
    self,
    refPoints: AtomPairsParameters,
    probePoints: AtomPairsParameters,
    weights: AtomPairsParameters = [],
    reflect: bool = False,
    maxIterations: int = 50,
) -> object:
    """
    Compute the optimal alignment (minimum RMSD) between two set of points

    ARGUMENTS:

    refPointsreference points specified as a N by 3 Numeric array or sequence of 3-sequences or sequence of Point3Ds

    probePointsprobe points to align to reference points - same format restrictions as reference points apply here

    weights : optional numeric vector or list of weights to associate to each pair of points
    reflect : reflect the probe points before attempting alignment
    maxIteration : maximum number of iterations to try to minimize RMSD

    RETURNS:

    a 2-tuple:
    SSD value for the alignment
    the 4x4 transform matrix, as a Numeric array

    C++ signature :_object* GetAlignmentTransform(boost::python::api::object,boost::python::api::object [,boost::python::api::object=[] [,bool=False [,unsigned int=50]]])
    """
    ...
