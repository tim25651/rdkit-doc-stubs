"""
rdkit.ML.Descriptors.Descriptors module¶
Various bits and pieces for calculating descriptors
"""
from _typeshed import Incomplete

class DescriptorCalculator(object):
    """
    abstract base class for descriptor calculators
    Constructor"""

    def CalcDescriptors(self, what, *args, **kwargs): ...
    simpleList: Incomplete
    descriptorNames: Incomplete
    compoundList: Incomplete

    def __init__(self, *args, **kwargs) -> None: ...
    def ShowDescriptors(self):
        """
        prints out a list of the descriptors"""
        ...
    def GetDescriptorNames(self):
        """
        returns a list of the names of the descriptors this calculator generates"""
        ...
    def SaveState(self, fileName):
        """
        Writes this calculator off to a file so that it can be easily loaded later
        Arguments

        fileName: the name of the file to be written"""
        ...
    def ShowDescriptors(self):
        """
        prints out a list of the descriptors"""
        ...
    def CalcDescriptors(self, what, *args, **kwargs): ...
