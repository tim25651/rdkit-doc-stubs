"""
rdkit.ML.Descriptors.MoleculeDescriptors module¶
Various bits and pieces for calculating Molecular descriptors
"""
from _typeshed import Incomplete
from rdkit.ML.Descriptors import Descriptors
from rdkit.ML.Descriptors import Descriptors as Descriptors
from rdkit.RDLogger import logger as logger

class MolecularDescriptorCalculator(Descriptors.DescriptorCalculator):
    """
    used for calculating descriptors for molecules
    Constructor
    Arguments

    simpleList: list of simple descriptors to be calculated(see below for format)

    Note

    format of simpleList:

    a list of strings which are functions in the rdkit.Chem.Descriptors module"""

    simpleList: Incomplete
    descriptorNames: Incomplete
    compoundList: Incomplete

    def __init__(self, simpleList, *args, **kwargs) -> None: ...
    def SaveState(self, fileName):
        """
        Writes this calculator off to a file so that it can be easily loaded later
        Arguments

        fileName: the name of the file to be written"""
        ...
    def CalcDescriptors(self, mol, *args, **kwargs):
        """
        calculates all descriptors for a given molecule
        Arguments

        mol: the molecule to be used

        Returnsa tuple of all descriptor values"""
        ...
    def GetDescriptorFuncs(self):
        """
        returns a tuple of the functions used to generate this calculator’s descriptors
        """
        ...
    def GetDescriptorNames(self):
        """
        returns a tuple of the names of the descriptors this calculator generates"""
        ...
    def GetDescriptorSummaries(self):
        """
        returns a tuple of summaries for the descriptors this calculator generates"""
        ...
    def GetDescriptorFuncs(self):
        """
        returns a tuple of the functions used to generate this calculator’s descriptors
        """
        ...
    def GetDescriptorVersions(self):
        """
        returns a tuple of the versions of the descriptor calculators"""
        ...
    def SaveState(self, fileName):
        """
        Writes this calculator off to a file so that it can be easily loaded later
        Arguments

        fileName: the name of the file to be written"""
        ...
