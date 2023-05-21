"""
rdkit.Chem.Fragments module¶
functions to match a bunch of fragment descriptors from a file
No user-servicable parts inside.  ;-)
"""
from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit import RDConfig as RDConfig
from rdkit.Chem import rdchem

def fr_Al_COO(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of aliphatic carboxylic acids"""
    ...

def fr_Al_OH(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of aliphatic hydroxyl groups"""
    ...

def fr_Al_OH_noTert(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of aliphatic hydroxyl groups excluding tert-OH"""
    ...

def fr_ArN(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of N functional groups attached to aromatics"""
    ...

def fr_Ar_COO(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of Aromatic carboxylic acide"""
    ...

def fr_Ar_N(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of aromatic nitrogens"""
    ...

def fr_Ar_NH(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of aromatic amines"""
    ...

def fr_Ar_OH(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of aromatic hydroxyl groups"""
    ...

def fr_COO(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of carboxylic acids"""
    ...

def fr_COO2(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of carboxylic acids"""
    ...

def fr_C_O(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of carbonyl O"""
    ...

def fr_C_O_noCOO(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of carbonyl O, excluding COOH"""
    ...

def fr_C_S(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of thiocarbonyl"""
    ...

def fr_HOCCN(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of C(OH)CCN-Ctert-alkyl or  C(OH)CCNcyclic"""
    ...

def fr_Imine(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of Imines"""
    ...

def fr_NH0(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of Tertiary amines"""
    ...

def fr_NH1(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of Secondary amines"""
    ...

def fr_NH2(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of Primary amines"""
    ...

def fr_N_O(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of hydroxylamine groups"""
    ...

def fr_Ndealkylation1(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of XCCNR groups"""
    ...

def fr_Ndealkylation2(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of tert-alicyclic amines (no heteroatoms, not quinine-like bridged N)"""
    ...

def fr_Nhpyrrole(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of H-pyrrole nitrogens"""
    ...

def fr_SH(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of thiol groups"""
    ...

def fr_aldehyde(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of aldehydes"""
    ...

def fr_alkyl_carbamate(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of alkyl carbamates (subject to hydrolysis)"""
    ...

def fr_alkyl_halide(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of alkyl halides"""
    ...

def fr_allylic_oxid(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of allylic oxidation sites excluding steroid dienone"""
    ...

def fr_amide(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of amides"""
    ...

def fr_amidine(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of amidine groups"""
    ...

def fr_aniline(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of anilines"""
    ...

def fr_aryl_methyl(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of aryl methyl sites for hydroxylation"""
    ...

def fr_azide(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of azide groups"""
    ...

def fr_azo(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of azo groups"""
    ...

def fr_barbitur(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of barbiturate groups"""
    ...

def fr_benzene(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of benzene rings"""
    ...

def fr_benzodiazepine(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of benzodiazepines with no additional fused rings"""
    ...

def fr_bicyclic(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Bicyclic"""
    ...

def fr_diazo(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of diazo groups"""
    ...

def fr_dihydropyridine(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of dihydropyridines"""
    ...

def fr_epoxide(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of epoxide rings"""
    ...

def fr_ester(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of esters"""
    ...

def fr_ether(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of ether oxygens (including phenoxy)"""
    ...

def fr_furan(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of furan rings"""
    ...

def fr_guanido(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of guanidine groups"""
    ...

def fr_halogen(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of halogens"""
    ...

def fr_hdrzine(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of hydrazine groups"""
    ...

def fr_hdrzone(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of hydrazone groups"""
    ...

def fr_imidazole(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of imidazole rings"""
    ...

def fr_imide(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of imide groups"""
    ...

def fr_isocyan(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of isocyanates"""
    ...

def fr_isothiocyan(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of isothiocyanates"""
    ...

def fr_ketone(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of ketones"""
    ...

def fr_ketone_Topliss(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of ketones excluding diaryl, a,b-unsat. dienones, heteroatom on Calpha"""
    ...

def fr_lactam(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of beta lactams"""
    ...

def fr_lactone(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of cyclic esters (lactones)"""
    ...

def fr_methoxy(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of methoxy groups -OCH3"""
    ...

def fr_morpholine(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of morpholine rings"""
    ...

def fr_nitrile(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of nitriles"""
    ...

def fr_nitro(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of nitro groups"""
    ...

def fr_nitro_arom(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of nitro benzene ring substituents"""
    ...

def fr_nitro_arom_nonortho(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of non-ortho nitro benzene ring substituents"""
    ...

def fr_nitroso(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of nitroso groups, excluding NO2"""
    ...

def fr_oxazole(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of oxazole rings"""
    ...

def fr_oxime(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of oxime groups"""
    ...

def fr_para_hydroxylation(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of para-hydroxylation sites"""
    ...

def fr_phenol(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of phenols"""
    ...

def fr_phenol_noOrthoHbond(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of phenolic OH excluding ortho intramolecular Hbond substituents"""
    ...

def fr_phos_acid(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of phosphoric acid groups"""
    ...

def fr_phos_ester(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of phosphoric ester groups"""
    ...

def fr_piperdine(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of piperdine rings"""
    ...

def fr_piperzine(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of piperzine rings"""
    ...

def fr_priamide(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of primary amides"""
    ...

def fr_prisulfonamd(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of primary sulfonamides"""
    ...

def fr_pyridine(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of pyridine rings"""
    ...

def fr_quatN(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of quaternary nitrogens"""
    ...

def fr_sulfide(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of thioether"""
    ...

def fr_sulfonamd(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of sulfonamides"""
    ...

def fr_sulfone(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of sulfone groups"""
    ...

def fr_term_acetylene(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of terminal acetylenes"""
    ...

def fr_tetrazole(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of tetrazole rings"""
    ...

def fr_thiazole(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of thiazole rings"""
    ...

def fr_thiocyan(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of thiocyanates"""
    ...

def fr_thiophene(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of thiophene rings"""
    ...

def fr_unbrch_alkane(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of unbranched alkanes  of at least 4 members (excludes halogenated alkanes)
    """
    ...

def fr_urea(self, mol, countUnique=True, pattern: rdchem.Mol = ...):
    """
    Number of urea groups"""
    ...

defaultPatternFileName: Incomplete
fns: Incomplete
fn: Incomplete
