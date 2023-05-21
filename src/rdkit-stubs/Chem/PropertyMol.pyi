"""
rdkit.Chem.PropertyMol module¶
"""
from rdkit import Chem as Chem
from rdkit.Chem import rdchem
from rdkit.Chem.rdchem import Mol

class PropertyMol(rdchem.Mol):
    """
    allows rdkit molecules to be pickled with their properties saved.
    >>> import os
    >>> import pickle
    >>> from rdkit import RDConfig
    >>> m = Chem.MolFromMolFile(os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data/benzene.mol'))
    >>> m.GetProp('_Name')
    'benzene.mol'

    by default pickling removes properties:
    >>> m2 = pickle.loads(pickle.dumps(m))
    >>> m2.HasProp('_Name')
    0

    Property mols solve this:
    >>> pm = PropertyMol(m)
    >>> pm.GetProp('_Name')
    'benzene.mol'
    >>> pm.SetProp('MyProp','foo')
    >>> pm.HasProp('MyProp')
    1

    >>> pm2 = pickle.loads(pickle.dumps(pm))
    >>> Chem.MolToSmiles(pm2)
    'c1ccccc1'
    >>> pm2.GetProp('_Name')
    'benzene.mol'
    >>> pm2.HasProp('MyProp')
    1
    >>> pm2.GetProp('MyProp')
    'foo'
    >>> pm2.HasProp('MissingProp')
    0

    Property mols are a bit more permissive about the types
    of property values:
    >>> pm.SetProp('IntVal',1)

    That wouldn’t work with a standard mol
    but the Property mols still convert all values to strings before storing:
    >>> pm.GetProp('IntVal')
    '1'

    This is a test for sf.net issue 2880943: make sure properties end up in SD files:
    >>> import tempfile, os
    >>> fn = tempfile.NamedTemporaryFile(suffix='.sdf', delete=False).name
    >>> w = Chem.SDWriter(fn)
    >>> w.write(pm)
    >>> w=None
    >>> with open(fn,'r') as inf:
    ...   txt = inf.read()
    >>> '<IntVal>' in txt
    True
    >>> try:
    ...   os.unlink(fn)
    ... except Exception:
    ...   pass

    The next level of that bug: does writing a depickled propertymol
    to an SD file include properties:
    >>> fn = tempfile.NamedTemporaryFile(suffix='.sdf', delete=False).name
    >>> w = Chem.SDWriter(fn)
    >>> pm = pickle.loads(pickle.dumps(pm))
    >>> w.write(pm)
    >>> w=None
    >>> with open(fn,'r') as inf:
    ...   txt = inf.read()
    >>> '<IntVal>' in txt
    True
    >>> try:
    ...   os.unlink(fn)
    ... except Exception:
    ...   pass

    Constructor, takes no arguments

    C++ signature :void __init__(_object*)

    __init__( (object)arg1, (str)pklString) -> None :

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)

    __init__( (object)arg1, (str)pklString, (int)propertyFlags) -> None :

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,unsigned int)

    __init__( (object)arg1, (Mol)mol [, (bool)quickCopy=False [, (int)confId=-1]]) -> None :

    C++ signature :void __init__(_object*,RDKit::ROMol [,bool=False [,int=-1]])"""

    __getstate_manages_dict__: bool

    def __init__(self, mol) -> None: ...
    def SetProp(self: Mol, key: str, val: str, computed: bool = False) -> None:
        """
        Sets a molecular property

        ARGUMENTS:
        key: the name of the property to be set (a string).
        value: the property value (a string).

        computed: (optional) marks the property as being computed.Defaults to False.

        C++ signature :void SetProp(RDKit::ROMol,char const*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])
        """
        ...
