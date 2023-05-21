"""
rdkit.Chem.RegistrationHash module¶
Generate a unique hash code for a molecule based on chemistry. If two
molecules are chemically “the same”, they should have the same hash.
Using molhash adds value beyond using SMILES because it:

Ignores SMILES features that are not chemically meaningful

(e.g. atom map numbers)
* Canonicalizes enhanced stereochemistry groups. For example
C[C@H](O)CC |&1:1| and C[C@@H](O)CC |&1:1| have the same
molhash
* Canonicalizes S group data (for example, polymer data)
There are two hash schemes, the default, and one in which
tautomers are considered equivalent.
"""
import enum
from typing import Iterable, Optional

from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit.Chem import rdMolHash as rdMolHash

class EnhancedStereoUpdateMode(enum.Enum):
    """
    An enumeration."""

    ADD_WEIGHTS: int = ...
    REMOVE_WEIGHTS: int = ...

ATOM_PROP_MAP_NUMBER: str
logger: Incomplete
DEFAULT_CXFLAG: Incomplete
ENHANCED_STEREO_GROUP_REGEX: Incomplete
ENHANCED_STEREO_GROUP_WEIGHTS: Incomplete
EMPTY_MOL_TAUTOMER_HASH: str

class HashLayer(enum.Enum):
    """
    Variables

    CANONICAL_SMILES – RDKit canonical SMILES (excluding enhanced stereo)
    ESCAPE – arbitrary other information to be incorporated
    FORMULA – a simple molecular formula for the molecule
    NO_STEREO_SMILES – RDKit canonical SMILES with all stereo removed
    SGROUP_DATA – canonicalization of all SGroups data present
    TAUTOMER_HASH – SMILES-like representation for a generic tautomer form
    NO_STEREO_TAUTOMER_HASH – the above tautomer hash lacking all stereo"""

    CANONICAL_SMILES: int = ...
    ESCAPE: int = ...
    FORMULA: int = ...
    NO_STEREO_SMILES: int = ...
    NO_STEREO_TAUTOMER_HASH: int = ...
    SGROUP_DATA: int = ...
    TAUTOMER_HASH: int = ...

class HashScheme(enum.Enum):
    """
    Which hash layers to use to when deduplicating molecules
    Typically the “ALL_LAYERS” scheme is used, but some users may want
    the “TAUTOMER_INSENSITIVE_LAYERS” scheme.

    Variables

    ALL_LAYERS – most strict hash scheme utilizing all layers
    STEREO_INSENSITIVE_LAYERS – excludes stereo sensitive layers
    TAUTOMER_INSENSITIVE_LAYERS – excludes tautomer sensitive layers"""

    ALL_LAYERS: tuple[HashLayer] = ...
    STEREO_INSENSITIVE_LAYERS: tuple[HashLayer] = ...
    TAUTOMER_INSENSITIVE_LAYERS: tuple[HashLayer] = ...

def GetMolHash(
    self, all_layers, hash_scheme: HashScheme = HashScheme.ALL_LAYERS
) -> str:
    """
    Generate a molecular hash using a specified set of layers.

    Parameters

    all_layers – a dictionary of layers
    hash_scheme – enum encoding information layers for the hash

    Returns
    hash for the given scheme constructed from the input layers"""
    ...

def GetMolLayers(self):
    """
    Generate layers of data about that could be used to identify a molecule

    Parameters

    original_molecule – molecule to obtain canonicalization layers from
    data_field_names – optional sequence of names of SGroup DAT fields which
    will be included in the hash.
    escape – optional field which can contain arbitrary information
    enable_tautomer_hash_v2 – use v2 of the tautomer hash

    Returns
    dictionary of HashLayer enum to calculated hash"""
    ...

def GetNoStereoLayers(self, mol, enable_tautomer_hash_v2=False): ...
def GetStereoTautomerHash(self, molecule, cxflag=65, enable_tautomer_hash_v2=False): ...
def GetNoStereoLayers(self, mol, enable_tautomer_hash_v2=False): ...

class EnhancedStereoUpdateMode(enum.Enum):
    """
    An enumeration."""

    ADD_WEIGHTS: int = ...
    REMOVE_WEIGHTS: int = ...
