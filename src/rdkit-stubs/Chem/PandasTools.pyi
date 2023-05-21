"""
rdkit.Chem.PandasTools module¶
Importing pandasTools enables several features that allow for using RDKit molecules as columns of a
Pandas dataframe.
If the dataframe is containing a molecule format in a column (e.g. smiles), like in this example:
>>> from rdkit.Chem import PandasTools
>>> import pandas as pd
>>> import os
>>> from rdkit import RDConfig
>>> antibiotics = pd.DataFrame(columns=['Name','Smiles'])
>>> antibiotics = pd.concat([antibiotics, pd.DataFrame.from_records([{'Smiles':'CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C',
...   'Name':'Penicilline G'}])], ignore_index=True) #Penicilline G
>>> antibiotics = pd.concat([antibiotics,pd.DataFrame.from_records([{
...   'Smiles':'CC1(C2CC3C(C(=O)C(=C(C3(C(=O)C2=C(C4=C1C=CC=C4O)O)O)O)C(=O)N)N(C)C)O',
...   'Name':'Tetracycline'}])], ignore_index=True) #Tetracycline
>>> antibiotics = pd.concat([antibiotics,pd.DataFrame.from_records([{
...   'Smiles':'CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=CC=C3)N)C(=O)O)C',
...   'Name':'Ampicilline'}])], ignore_index=True) #Ampicilline
>>> print([str(x) for x in  antibiotics.columns])
['Name', 'Smiles']
>>> print(antibiotics)
            Name                                             Smiles
0  Penicilline G    CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C
1   Tetracycline  CC1(C2CC3C(C(=O)C(=C(C3(C(=O)C2=C(C4=C1C=CC=C4...
2  Ampicilline  CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=CC=C3)N)C(=O...

a new column can be created holding the respective RDKit molecule objects. The fingerprint can be
included to accelerate substructure searches on the dataframe.
>>> PandasTools.AddMoleculeColumnToFrame(antibiotics,'Smiles','Molecule',includeFingerprints=True)
>>> print([str(x) for x in  antibiotics.columns])
['Name', 'Smiles', 'Molecule']

A substructure filter can be applied on the dataframe using the RDKit molecule column,
because the “>=” operator has been modified to work as a substructure check.
Such the antibiotics containing the beta-lactam ring “C1C(=O)NC1” can be obtained by
>>> beta_lactam = Chem.MolFromSmiles('C1C(=O)NC1')
>>> beta_lactam_antibiotics = antibiotics[antibiotics['Molecule'] >= beta_lactam]
>>> print(beta_lactam_antibiotics[['Name','Smiles']])
            Name                                             Smiles
0  Penicilline G    CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C
2  Ampicilline  CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=CC=C3)N)C(=O...

It is also possible to load an SDF file can be load into a dataframe.
>>> sdfFile = os.path.join(RDConfig.RDDataDir,'NCI/first_200.props.sdf')
>>> frame = PandasTools.LoadSDF(sdfFile,smilesName='SMILES',molColName='Molecule',
...            includeFingerprints=True)
>>> frame.info 
<bound method DataFrame.info of <class 'pandas.core.frame.DataFrame'>
Int64Index: 200 entries, 0 to 199
Data columns:
AMW                       200  non-null values
CLOGP                     200  non-null values
CP                        200  non-null values
CR                        200  non-null values
DAYLIGHT.FPG              200  non-null values
DAYLIGHT_CLOGP            200  non-null values
FP                        200  non-null values
ID                        200  non-null values
ISM                       200  non-null values
LIPINSKI_VIOLATIONS       200  non-null values
NUM_HACCEPTORS            200  non-null values
NUM_HDONORS               200  non-null values
NUM_HETEROATOMS           200  non-null values
NUM_LIPINSKIHACCEPTORS    200  non-null values
NUM_LIPINSKIHDONORS       200  non-null values
NUM_RINGS                 200  non-null values
NUM_ROTATABLEBONDS        200  non-null values
P1                        30  non-null values
SMILES                    200  non-null values
Molecule                  200  non-null values
dtypes: object(20)>

The standard ForwardSDMolSupplier keywords are also available:
>>> sdfFile = os.path.join(RDConfig.RDDataDir,'NCI/first_200.props.sdf')
>>> frame = PandasTools.LoadSDF(sdfFile, smilesName='SMILES', molColName='Molecule',
...            includeFingerprints=True, removeHs=False, strictParsing=True)

Conversion to html is quite easy:
>>> htm = frame.to_html() # doctest:
...
>>> str(htm[:36])
'<table border="1" class="dataframe">'

In order to support rendering the molecules as images in the HTML export of the
dataframe, we use a custom formatter for columns containing RDKit molecules,
and also disable escaping of HTML where needed.
"""
from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit import DataStructs as DataStructs
from rdkit.Chem import AllChem as AllChem
from rdkit.Chem import Draw as Draw
from rdkit.Chem import PandasPatcher as PandasPatcher
from rdkit.Chem import SDWriter as SDWriter
from rdkit.Chem import rdchem as rdchem
from rdkit.Chem.Draw.IPythonConsole import InteractiveRenderer as InteractiveRenderer
from rdkit.Chem.Draw.IPythonConsole import drawOptions as drawOptions
from rdkit.Chem.Scaffolds import MurckoScaffold as MurckoScaffold

def AddMoleculeColumnToFrame(
    self, frame, smilesCol="Smiles", molCol="ROMol", includeFingerprints=False
):
    """
    Converts the molecules contains in “smilesCol” to RDKit molecules and appends them to the
    dataframe “frame” using the specified column name.
    If desired, a fingerprint can be computed and stored with the molecule objects to accelerate
    substructure matching"""
    ...

def AddMurckoToFrame(
    self, frame, molCol="ROMol", MurckoCol="Murcko_SMILES", Generic=False
):
    """
    Adds column with SMILES of Murcko scaffolds to pandas DataFrame.
    Generic set to true results in SMILES of generic framework."""
    ...

def AlignMol(self, mol, scaffold):
    """
    Aligns mol (RDKit mol object) to scaffold (SMILES string)"""
    ...

def AlignToScaffold(self, frame, molCol="ROMol", scaffoldCol="Murcko_SMILES"):
    """
    Aligns molecules in molCol to scaffolds in scaffoldCol"""
    ...

def ChangeMoleculeRendering(self, frame=None, renderer="image"):
    """
    Allows to change the rendering of the molecules between image and string
    representations.
    This serves two purposes: First it allows to avoid the generation of images if this is
    not desired and, secondly, it allows to enable image rendering for newly created dataframe
    that already contains molecules, without having to rerun the time-consuming
    AddMoleculeColumnToFrame. Note: this behaviour is, because some pandas methods, e.g. head()
    returns a new dataframe instance that uses the default pandas rendering (thus not drawing
    images for molecules) instead of the monkey-patched one."""
    ...

def FrameToGridImage(self, frame, column="ROMol", legendsCol=None, **kwargs):
    """
    Draw grid image of mols in pandas DataFrame."""
    ...

def InstallPandasTools(self):
    """
    Monkey patch an RDKit method of Chem.Mol and pandas"""
    ...

def LoadSDF(
    self,
    filename,
    idName="ID",
    molColName="ROMol",
    includeFingerprints=False,
    isomericSmiles=True,
    smilesName=None,
    embedProps=False,
    removeHs=True,
    strictParsing=True,
):
    """
    Read file in SDF format and return as Pandas data frame.
    If embedProps=True all properties also get embedded in Mol objects in the molecule column.
    If molColName=None molecules would not be present in resulting DataFrame (only properties
    would be read)."""
    ...

InteractiveRenderer: Incomplete
drawOptions: Incomplete
log: Incomplete
highlightSubstructures: bool
molRepresentation: str
molSize: Incomplete
molJustify: str

def PrintAsImageString(self, x):
    """
    Returns the molecules as base64 encoded PNG image or as SVG"""
    ...

def RenderImagesInAllDataFrames(self, images=True):
    """
    Changes the default dataframe rendering to not escape HTML characters, thus allowing
    rendered images in all dataframes.
    IMPORTANT: THIS IS A GLOBAL CHANGE THAT WILL AFFECT TO COMPLETE PYTHON SESSION. If you want
    to change the rendering only for a single dataframe use the “ChangeMoleculeRendering” method
    instead."""
    ...

def LoadSDF(
    self,
    filename,
    idName="ID",
    molColName="ROMol",
    includeFingerprints=False,
    isomericSmiles=True,
    smilesName=None,
    embedProps=False,
    removeHs=True,
    strictParsing=True,
):
    """
    Read file in SDF format and return as Pandas data frame.
    If embedProps=True all properties also get embedded in Mol objects in the molecule column.
    If molColName=None molecules would not be present in resulting DataFrame (only properties
    would be read)."""
    ...

def RGroupDecompositionToFrame(
    self, groups, mols, include_core=False, redraw_sidechains=False
):
    """
    returns a dataframe with the results of R-Group Decomposition
    >>> from rdkit import Chem
    >>> from rdkit.Chem import rdRGroupDecomposition
    >>> from rdkit.Chem import PandasTools
    >>> import pandas as pd
    >>> scaffold = Chem.MolFromSmiles('c1ccccn1')
    >>> mols = [Chem.MolFromSmiles(smi) for smi in 'c1c(F)cccn1 c1c(Cl)c(C)ccn1 c1c(O)cccn1 c1c(F)c(C)ccn1 c1cc(Cl)c(F)cn1'.split()]
    >>> groups,_ = rdRGroupDecomposition.RGroupDecompose([scaffold],mols,asSmiles=False,asRows=False)
    >>> df = PandasTools.RGroupDecompositionToFrame(groups,mols,include_core=False)
    >>> list(df.columns)
    ['Mol', 'R1', 'R2']
    >>> df = PandasTools.RGroupDecompositionToFrame(groups,mols,include_core=True)
    >>> list(df.columns)
    ['Mol', 'Core', 'R1', 'R2']
    >>> len(df)
    5
    >>> df.columns()
    <class 'pandas*...*DataFrame'>
    RangeIndex: 5 entries, 0 to 4
    Data columns (total 4 columns):
    Mol     5 non-null object
    Core    5 non-null object
    R1      5 non-null object
    R2      5 non-null object
    dtypes: object(4)
    memory usage: *...*"""
    ...

def AddMoleculeColumnToFrame(
    self, frame, smilesCol="Smiles", molCol="ROMol", includeFingerprints=False
):
    """
    Converts the molecules contains in “smilesCol” to RDKit molecules and appends them to the
    dataframe “frame” using the specified column name.
    If desired, a fingerprint can be computed and stored with the molecule objects to accelerate
    substructure matching"""
    ...

def ChangeMoleculeRendering(self, frame=None, renderer="image"):
    """
    Allows to change the rendering of the molecules between image and string
    representations.
    This serves two purposes: First it allows to avoid the generation of images if this is
    not desired and, secondly, it allows to enable image rendering for newly created dataframe
    that already contains molecules, without having to rerun the time-consuming
    AddMoleculeColumnToFrame. Note: this behaviour is, because some pandas methods, e.g. head()
    returns a new dataframe instance that uses the default pandas rendering (thus not drawing
    images for molecules) instead of the monkey-patched one."""
    ...

def WriteSDF(
    self,
    df,
    out,
    molColName="ROMol",
    idName=None,
    properties=None,
    allNumeric=False,
    forceV3000=False,
):
    """
    Write an SD file for the molecules in the dataframe. Dataframe columns can be exported as
    SDF tags if specified in the “properties” list. “properties=list(df.columns)” would export
    all columns.
    The “allNumeric” flag allows to automatically include all numeric columns in the output.
    User has to make sure that correct data type is assigned to column.
    “idName” can be used to select a column to serve as molecule title. It can be set to
    “RowID” to use the dataframe row key as title."""
    ...

def RemoveSaltsFromFrame(self, frame, molCol="ROMol"):
    """
    Removes salts from mols in pandas DataFrame’s ROMol column"""
    ...

def RenderImagesInAllDataFrames(self, images=True):
    """
    Changes the default dataframe rendering to not escape HTML characters, thus allowing
    rendered images in all dataframes.
    IMPORTANT: THIS IS A GLOBAL CHANGE THAT WILL AFFECT TO COMPLETE PYTHON SESSION. If you want
    to change the rendering only for a single dataframe use the “ChangeMoleculeRendering” method
    instead."""
    ...

def SaveSMILESFromFrame(
    self, frame, outFile, molCol="ROMol", NamesCol="", isomericSmiles=False
):
    """
    Saves smi file. SMILES are generated from column with RDKit molecules. Column
    with names is optional."""
    ...

def SaveXlsxFromFrame(
    self, frame, outFile, molCol="ROMol", size=(300, 300), formats=None
):
    """
    Saves pandas DataFrame as a xlsx file with embedded images.
    molCol can be either a single column label or a list of column labels.
    It maps numpy data types to excel cell types:
    int, float -> number
    datetime -> datetime
    object -> string (limited to 32k character - xlsx limitations)
    The formats parameter can be optionally set to a dict of XlsxWriter
    formats (https://xlsxwriter.readthedocs.io/format.html#format), e.g.:
    {

    ‘write_string’:  {‘text_wrap’: True}

    }
    Currently supported keys for the formats dict are:
    ‘write_string’, ‘write_number’, ‘write_datetime’.
    Cells with compound images are a bit larger than images due to excel.
    Column width weirdness explained (from xlsxwriter docs):
    The width corresponds to the column width value that is specified in Excel.
    It is approximately equal to the length of a string in the default font of Calibri 11.
    Unfortunately, there is no way to specify “AutoFit” for a column in the Excel file format.
    This feature is only available at runtime from within Excel."""
    ...

def FrameToGridImage(self, frame, column="ROMol", legendsCol=None, **kwargs):
    """
    Draw grid image of mols in pandas DataFrame."""
    ...

def AddMurckoToFrame(
    self, frame, molCol="ROMol", MurckoCol="Murcko_SMILES", Generic=False
):
    """
    Adds column with SMILES of Murcko scaffolds to pandas DataFrame.
    Generic set to true results in SMILES of generic framework."""
    ...

def AlignMol(self, mol, scaffold):
    """
    Aligns mol (RDKit mol object) to scaffold (SMILES string)"""
    ...

def AlignToScaffold(self, frame, molCol="ROMol", scaffoldCol="Murcko_SMILES"):
    """
    Aligns molecules in molCol to scaffolds in scaffoldCol"""
    ...

def InstallPandasTools(self):
    """
    Monkey patch an RDKit method of Chem.Mol and pandas"""
    ...

def UninstallPandasTools(self):
    """
    Unpatch an RDKit method of Chem.Mol and pandas"""
    ...

def WriteSDF(
    self,
    df,
    out,
    molColName="ROMol",
    idName=None,
    properties=None,
    allNumeric=False,
    forceV3000=False,
):
    """
    Write an SD file for the molecules in the dataframe. Dataframe columns can be exported as
    SDF tags if specified in the “properties” list. “properties=list(df.columns)” would export
    all columns.
    The “allNumeric” flag allows to automatically include all numeric columns in the output.
    User has to make sure that correct data type is assigned to column.
    “idName” can be used to select a column to serve as molecule title. It can be set to
    “RowID” to use the dataframe row key as title."""
    ...
