"""
rdkit.Chem.MolKey.InchiInfo module¶
"""
from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit.Chem import inchi as inchi

console: Incomplete
UPD_APP: Incomplete
version_re: Incomplete
reconnected_re: Incomplete
fixed_h_re: Incomplete
isotope_re: Incomplete
stereo_re: Incomplete
stereo_all_re: Incomplete
undef_stereo_re: Incomplete
all_stereo_re: Incomplete
defined_stereo_re: Incomplete
h_layer_re: Incomplete
mobile_h_group_re: Incomplete
mobile_h_atoms_re: Incomplete

class InchiInfo(object):
    def get_mobile_h(self):
        """
        retrieve mobile H (tautomer) information
        return a 2-item tuple containing
        1) Number of mobile hydrogen groups detected. If 0, next item = ‘’
        2) List of groups"""
        ...
    parsed_inchi: Incomplete

    def __init__(self, inchi_str) -> None: ...
    def get_sp3_stereo(self):
        """
        retrieve sp3 stereo information
        return a 4-item tuple containing
        1) Number of stereocenters detected. If 0, the remaining items of the tuple = None
        2) Number of undefined stereocenters. Must be smaller or equal to above
        3) True if the molecule is a meso form (with chiral centers and a plane of symmetry)
        4) Comma-separated list of internal atom numbers with sp3 stereochemistry"""
        ...
    def get_mobile_h(self):
        """
        retrieve mobile H (tautomer) information
        return a 2-item tuple containing
        1) Number of mobile hydrogen groups detected. If 0, next item = ‘’
        2) List of groups"""
        ...
