from fr3d.modified_parent_mapping import modified_nucleotides

RNAbaseheavyatoms = {}
RNAbasehydrogens = {}
RNAconnections = {}

NAbaseheavyatoms = {}
NAbasehydrogens = {}
NAconnections = {}

Ribophos_connect = {}
tilt_cutoff = {}
aa_connections = {}
aa_backconnect = {}
aa_hydrogen_connections = {}
nt_phosphate = {}
nt_sugar = {}
nt_backbone = {}
HB_donors = {}
HB_acceptors = {}
HB_weak_donors = {}

#List of amino acids for perpendicular interactions

Perp_aa = set(['ARG','LYS','GLU','ASP','GLN','ASN','HIS','TYR','TRP','PHE'])

#Dictionaries for normal to plane calculations
planar_atoms = {}
planar_sugar = {}
sugar_plane = ["C2'","C1'","C3'"]

planar_sugar['A'] = sugar_plane
planar_sugar['U'] = sugar_plane
planar_sugar['C'] = sugar_plane
planar_sugar['G'] = sugar_plane

planar_atoms['A'] = ['C4','C5','N3']
planar_atoms['U'] = ['C2','N1','N3']
planar_atoms['C'] = ['C6','N1','C5']
planar_atoms['G'] = ['C4','C5','N3']

planar_sugar['DA'] = sugar_plane
planar_sugar['DT'] = sugar_plane
planar_sugar['DC'] = sugar_plane
planar_sugar['DG'] = sugar_plane

planar_atoms['DA'] = ['C4','C5','N3']
planar_atoms['U'] = ['C2','N1','N3']    # need DNA atoms for DT
planar_atoms['DC'] = ['C6','N1','C5']
planar_atoms['DG'] = ['C4','C5','N3']

# three atoms in certain amino acids used to calculate the plane of the amino acid
planar_atoms['ARG'] = ['CZ','NH1','NH2']
planar_atoms['ASN'] = ['CG','OD1','ND2']
planar_atoms['ASP'] = ['CG','OD1','OD2']
planar_atoms['CYS'] = ['CB','SG','CB']
planar_atoms['GLN'] = ['CD','OE1','NE2']
planar_atoms['GLU'] = ['CD','OE1','OE2']
planar_atoms['HIS'] = ['NE2','CD2','CE1']
planar_atoms['PHE'] = ['CG','CD1','CZ']
#planar_atoms['PRO'] = ['CA','N','CG']
planar_atoms['TRP'] = ['CD2','CG','CE2']
planar_atoms['TYR'] = ['CE2','CZ','CE1']
planar_atoms['LYS'] = ['CE','CD','NZ']         # need to verify
planar_atoms['SER'] = ['CB','OG','CA']         # need to verify
planar_atoms['THR'] = ['CB','OG1','CG2']       # need to verify

"""
planar_atoms['ILE'] =['CG1','CB','CD1']
planar_atoms['LEU'] =['CG','CD1','CD2']
planar_atoms['LYS'] = ['CE','CD','NZ']
planar_atoms['MET'] =['SD','CG','CE']
planar_atoms['VAL'] =['CB','CG1','CG2']"""

#Hydrogen bond donor and acceptor atoms for each RNA base

HB_donors['A'] = ["N6", "C2", "C8", "O2'"]
HB_donors['G'] = ["N1", "N2", "C8", "O2'"]
HB_donors['C'] = ["N4", "C5", "C6", "O2'"]
HB_donors['U'] = ["N3", "C5", "C6", "O2'"]

HB_acceptors['A'] = ["N1", "N3", "N7", "O2'"]
HB_acceptors['G'] = ["N3", "O6", "N7", "O2'"]
HB_acceptors['C'] = ["O2", "N3", "O2'"]
HB_acceptors['U'] = ["O2", "O4", "O2'"]

#Hydrogen bond donor and acceptor atoms for each NA backbone

HB_backbone_acceptors = ["O4'", "O3'", "O5'", "OP1", "OP2"]

#Hydrogen bond donor and acceptor atoms for each DNA base

HB_donors['DA'] = ["N6", "C2", "C8", "O2'"]
HB_donors['DG'] = ["N1", "N2", "C8", "O2'"]
HB_donors['DC'] = ["N4", "C5", "C6", "O2'"]
HB_donors['DT'] = ["N3", "C6", "O2'"]

HB_acceptors ['DA'] = ["N1", "N3", "N7", "O2'"]
HB_acceptors ['DG'] = ["N3", "O6", "N7", "O2'"]
HB_acceptors ['DC'] = ["O2", "N3", "O2'"]
HB_acceptors ['DT'] = ["O2", "O4", "O2'"]

# Hydrogen bond donor and acceptor atoms for each amino acid

HB_donors['ARG'] = ['NE', 'NH1', 'NH2']
HB_acceptors['ARG'] = []
HB_weak_donors['ARG'] = ['CA','CD']

HB_donors['ASN'] = ['ND2']
HB_acceptors['ASN'] = ['OD1']

#HB_donors['ASP'] = ['OD1', 'OD2']        # could have OH, exactly one of the oxygens could be donor
HB_donors['ASP'] = []                     # 99.9% of the time OD1, OD2 do not act as donors
HB_acceptors['ASP'] = ['OD1', 'OD2']

HB_donors['GLN'] = ['NE2']
HB_acceptors['GLN'] = ['OE1']

#HB_donors['GLU'] = ['OE1', 'OE2']        # could have OH, exactly one of the oxygens could be donor
HB_donors['GLU'] = []                     # 99.9% of the time OE1, OE2 do not act as donors
HB_acceptors['GLU'] = ['OE1', 'OE2']

HB_donors['HIS'] = ['ND1', 'NE2']
HB_weak_donors['HIS'] = ['CE1','CD2']
HB_acceptors['HIS'] = ['ND1', 'NE2']     # exactly one can be an acceptor

HB_donors['LYS'] = ['NZ']
HB_acceptors['LYS'] = []

HB_donors['PHE'] = ['CD1','CD2','CE1','CZ','CE2']
HB_acceptors['PHE'] =[]         # can we list the functional group center as a "base"?

HB_donors ['SER'] = ['OG']
HB_acceptors ['SER'] = ['OG']

HB_donors ['THR'] = ['OG1']
HB_acceptors ['THR'] = ['OG1']

HB_donors ['TRP'] = ['CD1','NE1','CE3','CZ3','CH2','CZ2']
HB_acceptors ['TRP'] = []       # can we list the functional group center as a "base"?

HB_donors ['TYR'] = ['OH', 'CD1','CD2', 'CE1', 'CE2']
HB_acceptors ['TYR'] = ['OH']   # can we list the functional group center as a "base"?
#Creating dictionaries for detecting edges of nts

# define base edges

WC_1 = {}
WC_2 = {}
Hoogsteen_1 = {}
Hoogsteen_2 = {}
Sugar_1 = {}
Sugar_2 = {}

# define base edge limits for RNA

WC_1['A']= ["C2"]
WC_2['A']= ["C6", "N6"]
Hoogsteen_1['A']= ["C6", "N6"]
Hoogsteen_2['A']= ["C8"]
Sugar_1['A']= ["N1","C1'"]
Sugar_2['A']= ["C2"]

WC_1['C']= ["C2", "O2"]
WC_2['C']= ["C4", "N4"]
Hoogsteen_1['C']= ["C4", "N4"]
Hoogsteen_2['C']= ["C5"]         # why not go all the way to C6?
Sugar_1['C']= ["N1","C1'"]
Sugar_2['C']= ["C2", "O2"]

WC_1['G']= ["C2", "N2"]
WC_2['G']= ["C6", "O6"]
Hoogsteen_1['G']= ["C6", "O6"]
Hoogsteen_2['G']= ["C8"]
Sugar_1['G']= ["N1","C1'"]
Sugar_2['G']= ["C2", "N2"]

WC_1['U']= ["C2", "O2"]
WC_2['U']= ["C4", "O4"]
Hoogsteen_1['U']= ["C4", "O4"]
Hoogsteen_2['U']= ["C5"]
Sugar_1['U']= ["N1","C1'"]
Sugar_2['U']= ["C2", "O2"]

# define base edge limits for DNA

WC_1['DA']= ["C2"]
WC_2['DA']= ["C6", "N6"]
Hoogsteen_1['DA']= ["C6", "N6"]
Hoogsteen_2['DA']= ["C8"]
Sugar_1['DA']= ["N1","C1'"]
Sugar_2['DA']= ["C2"]

WC_1['DC']= ["C2", "O2"]
WC_2['DC']= ["C4", "N4"]
Hoogsteen_1['DC']= ["C4", "N4"]
Hoogsteen_2['DC']= ["C5"]
Sugar_1['DC']= ["N1","C1'"]
Sugar_2['DC']= ["C2", "O2"]

WC_1['DG']= ["C2", "N2"]
WC_2['DG']= ["C6", "O6"]
Hoogsteen_1['DG']= ["C6", "O6"]
Hoogsteen_2['DG']= ["C8"]
Sugar_1['DG']= ["N1","C1'"]
Sugar_2['DG']= ["C2", "N2"]

WC_1['DT']= ["C2", "O2"]
WC_2['DT']= ["C4", "O4"]
Hoogsteen_1['DT']= ["C4", "O4"]
Hoogsteen_2['DT']= ["C5", "C7"]
Sugar_1['DT']= ["N1","C1'"]
Sugar_2['DT']= ["C2", "O2"]

#RNA computation definitions

"""Defining the parts of RNA nucleotides that we use to compute RNA-amino acid interactions"""

# This "RNA" variable is here for backward compatibility
RNAbaseheavyatoms['A'] = ['N9','C4','N3','N1','C6','N6','C8','C5','C2','N7']
RNAbasehydrogens['A'] = ['H2','H8','H9','H61','H62']
RNAbaseheavyatoms['C'] = ['N1','C2','O2','N3','C4','N4','C6','C5']
RNAbasehydrogens['C'] = ['H1','H6','H5','H41','H42']
RNAbaseheavyatoms['G'] = ['N9','C4','N3','N1','C6','O6','C8','C5','C2','N7','N2']
RNAbasehydrogens['G'] = ['H1','H8','H9','H22','H22']
RNAbaseheavyatoms['U'] = ['N1','C2','O2','N3','C4','O4','C6','C5']
RNAbasehydrogens['U'] = ['H5','H1','H3','H6']

# This should be the new standard, for all nucleic acids
NAbaseheavyatoms['A'] = ['N9','C4','N3','N1','C6','N6','C8','C5','C2','N7']
NAbaseheavyatoms['C'] = ['N1','C2','O2','N3','C4','N4','C6','C5']
NAbaseheavyatoms['G'] = ['N9','C4','N3','N1','C6','O6','C8','C5','C2','N7','N2']
NAbaseheavyatoms['U'] = ['N1','C2','O2','N3','C4','O4','C6','C5']

# RNA backbone
ribose = ["C1'","C2'","O2'","C3'","O3'","C4'","O4'","C5'"]
phosphate = ["O5'","P","OP1","OP2"]

nt_phosphate['A'] = phosphate
nt_phosphate['C'] = phosphate
nt_phosphate['G'] = phosphate
nt_phosphate['U'] = phosphate

nt_sugar['A'] = ribose
nt_sugar['C'] = ribose
nt_sugar['G'] = ribose
nt_sugar['U'] = ribose

nt_backbone['A'] = ribose + phosphate
nt_backbone['C'] = ribose + phosphate
nt_backbone['G'] = ribose + phosphate
nt_backbone['U'] = ribose + phosphate

NAbasehydrogens['A'] = ['H2','H8','H9','H61','H62']
NAbasehydrogens['C'] = ['H1','H6','H5','H41','H42']
NAbasehydrogens['G'] = ['H1','H8','H9','H21','H22']
NAbasehydrogens['U'] = ['H5','H1','H3','H6']

"""Defining the parts of DNA nucleotides"""

DNAribose = ["C1'","C2'","C3'","O3'","C4'","O4'","C5'"]

nt_sugar['DA'] = DNAribose
nt_sugar['DC'] = DNAribose
nt_sugar['DG'] = DNAribose
nt_sugar['DT'] = DNAribose

nt_phosphate['DA'] = phosphate
nt_phosphate['DC'] = phosphate
nt_phosphate['DG'] = phosphate
nt_phosphate['DT'] = phosphate

nt_backbone['DA'] = DNAribose + phosphate
nt_backbone['DC'] = DNAribose + phosphate
nt_backbone['DG'] = DNAribose + phosphate
nt_backbone['DT'] = DNAribose + phosphate

NAbaseheavyatoms['DA'] = ['N9','C4','N3','N1','C6','N6','C8','C5','C2','N7']
NAbaseheavyatoms['DC'] = ['N1','C2','O2','N3','C4','N4','C6','C5']
NAbaseheavyatoms['DG'] = ['N9','C4','N3','N1','C6','O6','C8','C5','C2','N7','N2']
NAbaseheavyatoms['DT'] = ['N1','C2','O2','N3','C4','O4','C6','C5','C7']

NAbasehydrogens['DA'] = ['H2','H8','H9','H61','H62']
NAbasehydrogens['DC'] = ['H1','H6','H5','H41','H42']
NAbasehydrogens['DG'] = ['H1','H8','H9','H21','H22']
NAbasehydrogens['DT'] = ['H3','H71','H72','H73','H6']

#Amino acid computation definitions

#Defining center-to-center and tilt cutoffs for stacking
tilt_cutoff= {'ALA': 2,'VAL': 0.7,'ILE': 1.9,'LEU': 2.1,'ARG': 1.5,'LYS': 1.5,'HIS': 1.2,'ASP': 1.5,'GLU': 1.5,'ASN': 1.4,'GLN': 1.4,'THR': 0.5,'SER': 0.5,'TYR': 2.1,'TRP': 2.1,'PHE': 1.5,'PRO': 3.1,'CYS': 1.0, 'MET': 1.5}

"""Defining the functional groups of sidechains of amino acids that interact
with nts by Hydrogenbonding. aa_fg refers to the functional group of the
sidechain. aa_backbone refers to the peptide backbone. aa_linker is the carbon chain
that links the fg with the peptide backbone"""

backbone = ['N','CA','C','O']

aa_backbone = {}
aa_linker = {}
aa_fg = {}

aa_backbone['ARG'] = backbone
aa_linker['ARG'] = ['CB','CG','CD']
aa_fg['ARG'] = ['NE','CZ','NH1','NH2']

aa_backbone['LYS'] = backbone
aa_linker['LYS'] = ['CB','CG','CD','CE']
aa_fg['LYS'] = ['NZ']

aa_backbone['HIS'] = backbone
aa_linker['HIS'] = ['CB']
aa_fg['HIS'] = ['CG','CD2','NE2','CE1','ND1']

aa_backbone['GLN'] = backbone
aa_linker['GLN'] = ['CB','CG']
aa_fg['GLN'] = ['CD','OE1','NE2']

aa_backbone['ASN'] = backbone
aa_linker['ASN'] = ['CB']
aa_fg['ASN'] = ['CG','OD1','ND2']

aa_backbone['GLU'] = backbone
aa_linker['GLU'] = ['CB','CG']
aa_fg['GLU'] = ['CD','OE1','OE2']

aa_backbone['ASP'] = backbone
aa_linker['ASP'] = ['CB']
aa_fg['ASP'] = ['CG','OD1','OD2']

aa_backbone['TRP'] = backbone
aa_linker['TRP'] = ['CB']
aa_fg['TRP'] = ['CG','CD1','NE1','CE2','CD2','CE3','CZ3','CH2','CZ2']

aa_backbone['TYR'] = backbone
aa_linker['TYR'] = ['CB']
aa_fg['TYR'] = ['CG','CD1','CE1','CZ','OH','CE2','CD2']

aa_backbone['PHE'] = backbone
aa_linker['PHE'] = ['CB']
aa_fg['PHE'] = ['CG','CD1','CD2','CE1','CZ','CE2']

aa_backbone['PRO'] = backbone
aa_linker['PRO'] = []
aa_fg['PRO'] = ['CB','CG','CD']

aa_backbone['MET'] = backbone
aa_linker['MET'] = []
aa_fg['MET'] = ['CB','CG','SD','CE']

aa_backbone['ILE'] = backbone
aa_linker['ILE'] = []
aa_fg['ILE'] = ['CB','CG1','CG2','CD1']

aa_backbone['LEU'] = backbone
aa_linker['LEU'] = []
aa_fg['LEU'] = ['CB','CG','CD2','CD1']

aa_backbone['VAL'] = backbone
aa_linker['VAL'] = []
aa_fg['VAL'] = ['CB','CG1','CG2']

aa_backbone['ALA'] = backbone
aa_linker['ALA'] = []
aa_fg['ALA'] = ['CB']

aa_backbone['GLY'] = backbone
aa_linker['GLY'] = []
aa_fg['GLY'] = []

aa_backbone['SER'] = backbone
aa_linker['SER'] = []
aa_fg['SER'] = ['CB','OG']

aa_backbone['THR'] = backbone
aa_linker['THR'] = ['CB']
aa_fg['THR'] = ['OG1','CG2']

aa_backbone['CYS'] = backbone
aa_linker['CYS'] = []
aa_fg['CYS'] = ['CB','SG']

# Quantum mechanics optimized base atom locations, by Jiri Sponer
# OK to change the listing of nucleotides to a numpy array, if that helps
# The numbers are generated by zStandardBases.m; the formatting can be changed there
# Base centers are placed at the origin (0,0,0)
RNAbasecoordinates = {}
RNAbasecoordinates['A'] = {}
RNAbasecoordinates['A'][ 'N9'] = [ -1.110515,  -1.823319,   0.000000]
RNAbasecoordinates['A'][ 'C4'] = [  0.007975,  -1.020192,   0.000000]
RNAbasecoordinates['A'][ 'N3'] = [  1.298514,  -1.383716,   0.000000]
RNAbasecoordinates['A'][ 'N1'] = [  1.754329,   1.001072,   0.000000]
RNAbasecoordinates['A'][ 'C6'] = [  0.453129,   1.320188,   0.000000]
RNAbasecoordinates['A'][ 'N6'] = [  0.092589,   2.623715,   0.000000]
RNAbasecoordinates['A'][ 'C8'] = [ -2.200426,  -0.989732,   0.000000]
RNAbasecoordinates['A'][ 'C5'] = [ -0.504186,   0.282479,   0.000000]
RNAbasecoordinates['A'][ 'C2'] = [  2.092069,  -0.307357,   0.000000]
RNAbasecoordinates['A'][ 'N7'] = [ -1.883477,   0.296861,   0.000000]
RNAbasecoordinates['A'][ 'H2'] = [  3.160932,  -0.505395,   0.000000]
RNAbasecoordinates['A'][ 'H8'] = [ -3.210978,  -1.376174,   0.000000]
RNAbasecoordinates['A'][ 'H9'] = [ -1.110515,  -2.832638,   0.000000]
RNAbasecoordinates['A']['1H6'] = [  0.810429,   3.327414,   0.000000]
RNAbasecoordinates['A']['2H6'] = [ -0.880080,   2.876400,   0.000000]
RNAbasecoordinates['C'] = {}
RNAbasecoordinates['C'][ 'N1'] = [ -0.380579,  -1.484583,   0.000000]
RNAbasecoordinates['C'][ 'C2'] = [  0.908756,  -0.885330,   0.000000]
RNAbasecoordinates['C'][ 'O2'] = [  1.888871,  -1.605366,   0.000000]
RNAbasecoordinates['C'][ 'N3'] = [  0.931558,   0.493192,   0.000000]
RNAbasecoordinates['C'][ 'C4'] = [ -0.197209,   1.170043,   0.000000]
RNAbasecoordinates['C'][ 'N4'] = [ -0.099223,   2.524619,   0.000000]
RNAbasecoordinates['C'][ 'C6'] = [ -1.542632,  -0.786228,   0.000000]
RNAbasecoordinates['C'][ 'C5'] = [ -1.509542,   0.573654,   0.000000]
RNAbasecoordinates['C'][ 'H1'] = [ -0.380579,  -2.494842,   0.000000]
RNAbasecoordinates['C'][ 'H6'] = [ -2.463312,  -1.360951,   0.000000]
RNAbasecoordinates['C'][ 'H5'] = [ -2.417351,   1.162512,   0.000000]
RNAbasecoordinates['C']['1H4'] = [ -0.908549,   3.117619,   0.000000]
RNAbasecoordinates['C']['2H4'] = [  0.823255,   2.927413,   0.000000]
RNAbasecoordinates['G'] = {}
RNAbasecoordinates['G'][ 'N9'] = [ -1.456680,  -1.711888,   0.000000]
RNAbasecoordinates['G'][ 'C4'] = [ -0.339093,  -0.921183,   0.000000]
RNAbasecoordinates['G'][ 'N3'] = [  0.947138,  -1.370892,   0.000000]
RNAbasecoordinates['G'][ 'N1'] = [  1.440727,   0.946157,   0.000000]
RNAbasecoordinates['G'][ 'C6'] = [  0.110442,   1.479180,   0.000000]
RNAbasecoordinates['G'][ 'O6'] = [ -0.059694,   2.682325,   0.000000]
RNAbasecoordinates['G'][ 'C8'] = [ -2.545320,  -0.870005,   0.000000]
RNAbasecoordinates['G'][ 'C5'] = [ -0.831466,   0.385014,   0.000000]
RNAbasecoordinates['G'][ 'C2'] = [  1.802732,  -0.378447,   0.000000]
RNAbasecoordinates['G'][ 'N7'] = [ -2.208072,   0.408244,   0.000000]
RNAbasecoordinates['G'][ 'N2'] = [  3.139287,  -0.648506,   0.000000]
RNAbasecoordinates['G'][ 'H1'] = [  2.156655,   1.662659,   0.000000]
RNAbasecoordinates['G'][ 'H8'] = [ -3.558253,  -1.248304,   0.000000]
RNAbasecoordinates['G'][ 'H9'] = [ -1.456680,  -2.721417,   0.000000]
RNAbasecoordinates['G']['1H2'] = [  3.415836,  -1.614107,   0.000000]
RNAbasecoordinates['G']['2H2'] = [  3.838495,   0.070395,   0.000000]
RNAbasecoordinates['U'] = {}
RNAbasecoordinates['U'][ 'N1'] = [ -0.326420,  -1.514422,   0.000000]
RNAbasecoordinates['U'][ 'C2'] = [  0.933866,  -0.925171,   0.000000]
RNAbasecoordinates['U'][ 'O2'] = [  1.966070,  -1.562816,   0.000000]
RNAbasecoordinates['U'][ 'N3'] = [  0.868896,   0.459470,   0.000000]
RNAbasecoordinates['U'][ 'C4'] = [ -0.266417,   1.293935,   0.000000]
RNAbasecoordinates['U'][ 'O4'] = [ -0.149264,   2.505638,   0.000000]
RNAbasecoordinates['U'][ 'C6'] = [ -1.504181,  -0.804727,   0.000000]
RNAbasecoordinates['U'][ 'C5'] = [ -1.522550,   0.548093,   0.000000]
RNAbasecoordinates['U'][ 'H5'] = [ -2.450997,   1.102115,   0.000000]
RNAbasecoordinates['U'][ 'H1'] = [ -0.326420,  -2.523369,   0.000000]
RNAbasecoordinates['U'][ 'H3'] = [  1.765732,   0.930757,   0.000000]
RNAbasecoordinates['U'][ 'H6'] = [ -2.409200,  -1.402586,   0.000000]

# RNA base coordinates in the xy plane, with heavy atoms centered at (0,0,0)
NAbasecoordinates = {}
NAbasecoordinates['A'] = {}
NAbasecoordinates['A'][ 'N9'] = [ -1.110515,  -1.823319,   0.000000]
NAbasecoordinates['A'][ "C1'"] = [ -1.110515,  (-1.823319 - 1.48),   0.000000] #should be roughly 1.48 angstroms lower than N9 
NAbasecoordinates['A'][ 'C4'] = [  0.007975,  -1.020192,   0.000000]
NAbasecoordinates['A'][ 'N3'] = [  1.298514,  -1.383716,   0.000000]
NAbasecoordinates['A'][ 'N1'] = [  1.754329,   1.001072,   0.000000]
NAbasecoordinates['A'][ 'C6'] = [  0.453129,   1.320188,   0.000000]
NAbasecoordinates['A'][ 'N6'] = [  0.092589,   2.623715,   0.000000]
NAbasecoordinates['A'][ 'C8'] = [ -2.200426,  -0.989732,   0.000000]
NAbasecoordinates['A'][ 'C5'] = [ -0.504186,   0.282479,   0.000000]
NAbasecoordinates['A'][ 'C2'] = [  2.092069,  -0.307357,   0.000000]
NAbasecoordinates['A'][ 'N7'] = [ -1.883477,   0.296861,   0.000000]
NAbasecoordinates['A'][ 'H2'] = [  3.160932,  -0.505395,   0.000000]
NAbasecoordinates['A'][ 'H8'] = [ -3.210978,  -1.376174,   0.000000]
NAbasecoordinates['A'][ 'H9'] = [ -1.110515,  -2.832638,   0.000000]
NAbasecoordinates['A']['H61'] = [  0.810429,   3.327414,   0.000000]
NAbasecoordinates['A']['H62'] = [ -0.880080,   2.876400,   0.000000]

NAbasecoordinates['C'] = {}
NAbasecoordinates['C'][ 'N1'] = [ -0.380579,  -1.484583,   0.000000]
NAbasecoordinates['C'][ "C1'"] = [ -0.380579,  (-1.484583-1.48),   0.000000] #should be roughly 1.48 angstroms lower than N1 
NAbasecoordinates['C'][ 'C2'] = [  0.908756,  -0.885330,   0.000000]
NAbasecoordinates['C'][ 'O2'] = [  1.888871,  -1.605366,   0.000000]
NAbasecoordinates['C'][ 'N3'] = [  0.931558,   0.493192,   0.000000]
NAbasecoordinates['C'][ 'C4'] = [ -0.197209,   1.170043,   0.000000]
NAbasecoordinates['C'][ 'N4'] = [ -0.099223,   2.524619,   0.000000]
NAbasecoordinates['C'][ 'C6'] = [ -1.542632,  -0.786228,   0.000000]
NAbasecoordinates['C'][ 'C5'] = [ -1.509542,   0.573654,   0.000000]
NAbasecoordinates['C'][ 'H1'] = [ -0.380579,  -2.494842,   0.000000]
NAbasecoordinates['C'][ 'H6'] = [ -2.463312,  -1.360951,   0.000000]
NAbasecoordinates['C'][ 'H5'] = [ -2.417351,   1.162512,   0.000000]
NAbasecoordinates['C']['H41'] = [ -0.908549,   3.117619,   0.000000]
NAbasecoordinates['C']['H42'] = [  0.823255,   2.927413,   0.000000]

NAbasecoordinates['G'] = {}
NAbasecoordinates['G'][ 'N9'] = [ -1.456680,  -1.711888,   0.000000]
NAbasecoordinates['G'][ "C1'"] = [ -1.456680,  (-1.711888-1.48),   0.000000] #should be roughly 1.48 angstroms lower than N9 
NAbasecoordinates['G'][ 'C4'] = [ -0.339093,  -0.921183,   0.000000]
NAbasecoordinates['G'][ 'N3'] = [  0.947138,  -1.370892,   0.000000]
NAbasecoordinates['G'][ 'N1'] = [  1.440727,   0.946157,   0.000000]
NAbasecoordinates['G'][ 'C6'] = [  0.110442,   1.479180,   0.000000]
NAbasecoordinates['G'][ 'O6'] = [ -0.059694,   2.682325,   0.000000]
NAbasecoordinates['G'][ 'C8'] = [ -2.545320,  -0.870005,   0.000000]
NAbasecoordinates['G'][ 'C5'] = [ -0.831466,   0.385014,   0.000000]
NAbasecoordinates['G'][ 'C2'] = [  1.802732,  -0.378447,   0.000000]
NAbasecoordinates['G'][ 'N7'] = [ -2.208072,   0.408244,   0.000000]
NAbasecoordinates['G'][ 'N2'] = [  3.139287,  -0.648506,   0.000000]
NAbasecoordinates['G'][ 'H1'] = [  2.156655,   1.662659,   0.000000]
NAbasecoordinates['G'][ 'H8'] = [ -3.558253,  -1.248304,   0.000000]
NAbasecoordinates['G'][ 'H9'] = [ -1.456680,  -2.721417,   0.000000]  # if no C1' atom
NAbasecoordinates['G']['H21'] = [  3.415836,  -1.614107,   0.000000]
NAbasecoordinates['G']['H22'] = [  3.838495,   0.070395,   0.000000]

NAbasecoordinates['U'] = {}
NAbasecoordinates['U'][ 'N1'] = [ -0.326420,  -1.514422,   0.000000]
NAbasecoordinates['U'][ "C1'"] = [ -0.326420,  (-1.514422-1.48),   0.000000] #should be roughly 1.48 angstroms lower than N1
NAbasecoordinates['U'][ 'C2'] = [  0.933866,  -0.925171,   0.000000]
NAbasecoordinates['U'][ 'O2'] = [  1.966070,  -1.562816,   0.000000]
NAbasecoordinates['U'][ 'N3'] = [  0.868896,   0.459470,   0.000000]
NAbasecoordinates['U'][ 'C4'] = [ -0.266417,   1.293935,   0.000000]
NAbasecoordinates['U'][ 'O4'] = [ -0.149264,   2.505638,   0.000000]
NAbasecoordinates['U'][ 'C6'] = [ -1.504181,  -0.804727,   0.000000]
NAbasecoordinates['U'][ 'C5'] = [ -1.522550,   0.548093,   0.000000]
NAbasecoordinates['U'][ 'H5'] = [ -2.450997,   1.102115,   0.000000]
NAbasecoordinates['U'][ 'H1'] = [ -0.326420,  -2.523369,   0.000000]
NAbasecoordinates['U'][ 'H3'] = [  1.765732,   0.930757,   0.000000]
NAbasecoordinates['U'][ 'H6'] = [ -2.409200,  -1.402586,   0.000000]

NAbasecoordinates['DA'] = {}
NAbasecoordinates['DA'][ 'N9'] = [ -1.110515,  -1.823319,   0.000000]
NAbasecoordinates['DA'][ "C1'"] = [ -1.110515,  (-1.823319-1.48),   0.000000] #should be roughly 1.48 angstroms lower than N9
NAbasecoordinates['DA'][ 'C4'] = [  0.007975,  -1.020192,   0.000000]
NAbasecoordinates['DA'][ 'N3'] = [  1.298514,  -1.383716,   0.000000]
NAbasecoordinates['DA'][ 'N1'] = [  1.754329,   1.001072,   0.000000]
NAbasecoordinates['DA'][ 'C6'] = [  0.453129,   1.320188,   0.000000]
NAbasecoordinates['DA'][ 'N6'] = [  0.092589,   2.623715,   0.000000]
NAbasecoordinates['DA'][ 'C8'] = [ -2.200426,  -0.989732,   0.000000]
NAbasecoordinates['DA'][ 'C5'] = [ -0.504186,   0.282479,   0.000000]
NAbasecoordinates['DA'][ 'C2'] = [  2.092069,  -0.307357,   0.000000]
NAbasecoordinates['DA'][ 'N7'] = [ -1.883477,   0.296861,   0.000000]
NAbasecoordinates['DA'][ 'H2'] = [  3.160932,  -0.505395,   0.000000]
NAbasecoordinates['DA'][ 'H8'] = [ -3.210978,  -1.376174,   0.000000]
NAbasecoordinates['DA'][ 'H9'] = [ -1.110515,  -2.832638,   0.000000]
NAbasecoordinates['DA']['H61'] = [  0.810429,   3.327414,   0.000000]
NAbasecoordinates['DA']['H62'] = [ -0.880080,   2.876400,   0.000000]
NAbasecoordinates['DC'] = {}
NAbasecoordinates['DC'][ 'N1'] = [ -0.380579,  -1.484583,   0.000000]
NAbasecoordinates['DC'][ "C1'"] = [ -0.380579,  (-1.484583-1.48),   0.000000] #should be roughly 1.48 angstroms lower than N1
NAbasecoordinates['DC'][ 'C2'] = [  0.908756,  -0.885330,   0.000000]
NAbasecoordinates['DC'][ 'O2'] = [  1.888871,  -1.605366,   0.000000]
NAbasecoordinates['DC'][ 'N3'] = [  0.931558,   0.493192,   0.000000]
NAbasecoordinates['DC'][ 'C4'] = [ -0.197209,   1.170043,   0.000000]
NAbasecoordinates['DC'][ 'N4'] = [ -0.099223,   2.524619,   0.000000]
NAbasecoordinates['DC'][ 'C6'] = [ -1.542632,  -0.786228,   0.000000]
NAbasecoordinates['DC'][ 'C5'] = [ -1.509542,   0.573654,   0.000000]
NAbasecoordinates['DC'][ 'H1'] = [ -0.380579,  -2.494842,   0.000000]
NAbasecoordinates['DC'][ 'H6'] = [ -2.463312,  -1.360951,   0.000000]
NAbasecoordinates['DC'][ 'H5'] = [ -2.417351,   1.162512,   0.000000]
NAbasecoordinates['DC']['H41'] = [ -0.908549,   3.117619,   0.000000]
NAbasecoordinates['DC']['H42'] = [  0.823255,   2.927413,   0.000000]
NAbasecoordinates['DG'] = {}
NAbasecoordinates['DG'][ 'N9'] = [ -1.456680,  -1.711888,   0.000000]
NAbasecoordinates['DG'][ "C1'"] = [ -1.456680,  (-1.711888-1.48),   0.000000] #should be roughly 1.48 angstroms lower than N9
NAbasecoordinates['DG'][ 'C4'] = [ -0.339093,  -0.921183,   0.000000]
NAbasecoordinates['DG'][ 'N3'] = [  0.947138,  -1.370892,   0.000000]
NAbasecoordinates['DG'][ 'N1'] = [  1.440727,   0.946157,   0.000000]
NAbasecoordinates['DG'][ 'C6'] = [  0.110442,   1.479180,   0.000000]
NAbasecoordinates['DG'][ 'O6'] = [ -0.059694,   2.682325,   0.000000]
NAbasecoordinates['DG'][ 'C8'] = [ -2.545320,  -0.870005,   0.000000]
NAbasecoordinates['DG'][ 'C5'] = [ -0.831466,   0.385014,   0.000000]
NAbasecoordinates['DG'][ 'C2'] = [  1.802732,  -0.378447,   0.000000]
NAbasecoordinates['DG'][ 'N7'] = [ -2.208072,   0.408244,   0.000000]
NAbasecoordinates['DG'][ 'N2'] = [  3.139287,  -0.648506,   0.000000]
NAbasecoordinates['DG'][ 'H1'] = [  2.156655,   1.662659,   0.000000]
NAbasecoordinates['DG'][ 'H8'] = [ -3.558253,  -1.248304,   0.000000]
NAbasecoordinates['DG'][ 'H9'] = [ -1.456680,  -2.721417,   0.000000]
NAbasecoordinates['DG']['H21'] = [  3.415836,  -1.614107,   0.000000]
NAbasecoordinates['DG']['H22'] = [  3.838495,   0.070395,   0.000000]

# improvised in fall 2020, should be improved in the near future
NAbasecoordinates['DT'] = {}
NAbasecoordinates['DT']['N1'] = [ 0.077892,-1.669793, 0.000000]
NAbasecoordinates['DT']["C1'"] = [ 0.077892,(-1.669793-1.48), 0.000000] #should be roughly 1.48 angstroms lower than N1
NAbasecoordinates['DT']['C2'] = [ 1.276471,-0.994656, 0.000000]
NAbasecoordinates['DT']['O2'] = [ 2.359169,-1.553674, 0.000000]
NAbasecoordinates['DT']['N3'] = [ 1.166654, 0.370792, 0.000000]
NAbasecoordinates['DT']['C4'] = [ 0.001313, 1.113698, 0.000000]
NAbasecoordinates['DT']['O4'] = [ 0.005353, 2.340664, 0.000000]
NAbasecoordinates['DT']['C5'] = [-1.219974, 0.345949, 0.000000]
NAbasecoordinates['DT']['C7'] = [-2.539534, 1.039263, 0.000000]
NAbasecoordinates['DT']['C6'] = [-1.127344,-0.992242, 0.000000]
NAbasecoordinates['DT']['H3'] = [ 1.903148, 0.814431, 0.000000]
NAbasecoordinates['DT']['H6'] = [-1.913009,-1.488609, 0.000000]
NAbasecoordinates['DT']['H71'] = [-2.422619, 1.940262, 0.339933]
NAbasecoordinates['DT']['H72'] = [-2.882094, 1.068812,-0.906410]
NAbasecoordinates['DT']['H73'] = [-3.153579, 0.548710, 0.569770]

# Coordinates of atoms as defined by the matlab FR3D code in zStackingOverlap.m
# NTs are all translated to try to match the C1' x,y coordinates to work with Draw Base Method from draw_residues 
# These coordinates are used to create checks in generate_location_checks which are then used in checkconvexhull.py
# The M is stands for Matlab and the second letter stands for the RNA base
# The indeces can probably be renamed to atoms 
NAbasecoordinates['MA'] = {}
NAbasecoordinates['MA']['0'] = [(-2.100463 - 1.110515), (0.447145 - 2.293999) ,0.000000] 
NAbasecoordinates['MA']['1'] = [(0.000000 - 1.110515), (-1.009320 - 2.293999) , 0.000000]
NAbasecoordinates['MA']['2'] = [(4.271447 - 1.110515), (1.317924 - 2.293999) ,0.000000]
NAbasecoordinates['MA']['3'] = [(1.920945 - 1.110515), (5.150733 - 2.293999) ,0.000000]
NAbasecoordinates['MA']['4'] = [(0.230436 - 1.110515), (4.699718 - 2.293999) ,0.000000]
NAbasecoordinates['MA']['5'] = [(-2.100463 - 1.110515), (0.447145 - 2.293999) ,0.00000]

NAbasecoordinates['MC'] = {}
NAbasecoordinates['MC']['0'] = [(-2.082733 - 0.380579), (0.123632 - 1.954324), 0.000000] 
NAbasecoordinates['MC']['1'] = [(0.000000 - 0.380579), (-1.010259 - 1.954324), 0.000000] 
NAbasecoordinates['MC']['2'] = [(2.269450 - 0.380579), (-0.120783 - 1.954324), 0.000000]
NAbasecoordinates['MC']['3'] = [(1.203833 - 0.380579), (4.41199 - 1.954324) , 0.000000]
NAbasecoordinates['MC']['4'] = [(-0.527970 - 0.380579), (4.602202 - 1.954324) , 0.000000]
NAbasecoordinates['MC']['5'] = [(-2.036772 - 0.380579), (2.647095 - 1.954324) , 0.00000]
NAbasecoordinates['MC']['6'] = [(-2.082733 - 0.380579), (0.123632 - 1.954324) , 0.000000]

NAbasecoordinates['MG'] = {}
NAbasecoordinates['MG']['0'] = [(-2.101572 - 1.456680), (0.463584 - 2.182359), 0.000000] 
NAbasecoordinates['MG']['1'] = [(0.000000 - 1.456680), (-1.009529 - 2.182359), 0.000000]
NAbasecoordinates['MG']['2'] = [(4.872516 - 1.456680), (0.0977814 - 2.182359), 0.000000]
NAbasecoordinates['MG']['3'] = [(5.295175 - 1.456680), (1.782283 - 2.182359), 0.000000]
NAbasecoordinates['MG']['4'] = [(3.613335 - 1.456680), (3.374547 - 2.182359), 0.000000]
NAbasecoordinates['MG']['5'] = [(1.396986 - 1.456680), (4.394213 - 2.182359), 0.00000]
NAbasecoordinates['MG']['6'] = [(-0.751391 - 1.456680), (2.120132 - 2.182359), 0.00000]
NAbasecoordinates['MG']['7'] = [(-2.101572 - 1.456680), (0.463584 - 2.182359), 0.000000]

NAbasecoordinates['MU'] = {}
NAbasecoordinates['MU']['0'] = [(-2.082780 - 0.326420), (0.111836 - 1.985475), 0.000000] 
NAbasecoordinates['MU']['1'] = [(0.000000 - 0.326420), (-1.008947 - 1.985475), 0.000000]
NAbasecoordinates['MU']['2'] = [(2.292490 - 0.326420), (-0.048394 - 1.985475), 0.000000]
NAbasecoordinates['MU']['3'] = [(2.092152 - 0.326420), (2.445179 - 1.985475), 0.000000]
NAbasecoordinates['MU']['4'] = [(0.177156 - 0.326420), (4.020060 - 1.985475), 0.000000]
NAbasecoordinates['MU']['5'] = [(-2.124577 - 0.326420), (2.616537 - 1.985475), 0.00000]
NAbasecoordinates['MU']['6'] = [(-2.082780 - 0.326420), (0.111836 - 1.985475), 0.00000]



#Definitions for drawing nucleotides
NAbasecolor = {}
NAbasecolor['A'] = [1,0,0]   # red
NAbasecolor['C'] = [1,214.0/255,0]   # yellowish
NAbasecolor['G'] = [0.5,1,0]   # green
NAbasecolor['U'] = [0,1,1]   # cyan
NAbasecolor['DA'] = [1,0,0]   # red
NAbasecolor['DC'] = [1,214.0/255,0]   # yellowish
NAbasecolor['DG'] = [0.5,1,0]   # green
NAbasecolor['DT'] = [0,0,1]   # blue


RNAconnections['A'] =['N1','C6','C6','N6','C6','C5','C5','C4','C5','N7','N7','C8','C8','N9','N9','C4','C4','C5','C4','N3','N3','C2','C2','N1']
RNAconnections['U'] =['N1','C2','C2','O2','C2','N3','N3','C4','C4','O4','C4','C5','C5','C6','C6','N1']
RNAconnections['G'] =['N1','C6','C6','O6','C6','C5','C5','C4','C5','N7','N7','C8','C8','N9','N9','C4','C4','C5','C4','N3','N3','C2','C2','N2','C2','N1']
RNAconnections['C'] =['N1','C2','C2','O2','C2','N3','N3','C4','C4','N4','C4','C5','C5','C6','C6','N1']

Ribophos_connect['A'] = ["N9","C1'","C1'","C2'","C2'","O2'","C2'","C3'","C3'","O3'","C3'","C4'","C4'","O4'","O4'","C1'","O4'","C4'","C4'","C5'","C5'","O5'","O5'","P","P","OP1","P","OP2"]
Ribophos_connect['U'] = ["N1","C1'","C1'","C2'","C2'","O2'","C2'","C3'","C3'","O3'","C3'","C4'","C4'","O4'","O4'","C1'","O4'","C4'","C4'","C5'","C5'","O5'","O5'","P","P","OP1","P","OP2"]
Ribophos_connect['G'] = ["N9","C1'","C1'","C2'","C2'","O2'","C2'","C3'","C3'","O3'","C3'","C4'","C4'","O4'","O4'","C1'","O4'","C4'","C4'","C5'","C5'","O5'","O5'","P","P","OP1","P","OP2"]
Ribophos_connect['C'] = ["N1","C1'","C1'","C2'","C2'","O2'","C2'","C3'","C3'","O3'","C3'","C4'","C4'","O4'","O4'","C1'","O4'","C4'","C4'","C5'","C5'","O5'","O5'","P","P","OP1","P","OP2"]

NAconnections[ 'A'] =['N1','C6','C6','N6','C6','C5','C5','C4','C5','N7','N7','C8','C8','N9','N9','C4','C4','C5','C4','N3','N3','C2','C2','N1']
NAconnections[ 'C'] =['N1','C2','C2','O2','C2','N3','N3','C4','C4','N4','C4','C5','C5','C6','C6','N1']
NAconnections[ 'G'] =['N1','C6','C6','O6','C6','C5','C5','C4','C5','N7','N7','C8','C8','N9','N9','C4','C4','C5','C4','N3','N3','C2','C2','N2','C2','N1']
NAconnections[ 'U'] =['N1','C2','C2','O2','C2','N3','N3','C4','C4','O4','C4','C5','C5','C6','C6','N1']
NAconnections['DA'] =['N1','C6','C6','N6','C6','C5','C5','C4','C5','N7','N7','C8','C8','N9','N9','C4','C4','C5','C4','N3','N3','C2','C2','N1']
NAconnections['DC'] =['N1','C2','C2','O2','C2','N3','N3','C4','C4','N4','C4','C5','C5','C6','C6','N1']
NAconnections['DG'] =['N1','C6','C6','O6','C6','C5','C5','C4','C5','N7','N7','C8','C8','N9','N9','C4','C4','C5','C4','N3','N3','C2','C2','N2','C2','N1']
NAconnections['DT'] =['N1','C2','C2','N3','N3','C4','C3','O4','C4','C5','C5','C7','C5','C6','C6','N1']

#NAhydrogenconnections = {}
#NAhydrogenconnections[ 'A'] =[]
#NAhydrogenconnections[ 'C'] =[]
#NAhydrogenconnections[ 'G'] =[]
#NAhydrogenconnections[ 'U'] =[]

# Definitions for drawing the amino acid backbones.  Each pair gets a line between them.
aa_backconnect['ARG']=['N','CA','CA','C','C','O','C','CA']
aa_backconnect['LYS']=['N','CA','CA','C','C','O','C','CA']
aa_backconnect['HIS']=['N','CA','CA','C','C','O','C','CA']
aa_backconnect['GLN']=['N','CA','CA','C','C','O','C','CA']
aa_backconnect['ASN']=['N','CA','CA','C','C','O','C','CA']
aa_backconnect['ASP']=['N','CA','CA','C','C','O','C','CA']
aa_backconnect['GLU']=['N','CA','CA','C','C','O','C','CA']
aa_backconnect['TRP']=['N','CA','CA','C','C','O','C','CA']
aa_backconnect['TYR']=['N','CA','CA','C','C','O','C','CA']
aa_backconnect['PHE']=['N','CA','CA','C','C','O','C','CA']
aa_backconnect['PRO']=['N','CA','CA','C','C','O','C','CA']
aa_backconnect['MET']=['N','CA','CA','C','C','O','C','CA']
aa_backconnect['ILE']=['N','CA','CA','C','C','O','C','CA']
aa_backconnect['LEU']=['N','CA','CA','C','C','O','C','CA']
aa_backconnect['ALA']=['N','CA','CA','C','C','O','C','CA']
aa_backconnect['VAL']=['N','CA','CA','C','C','O','C','CA']
aa_backconnect['GLY']=['N','CA','CA','C','C','O','C','CA']
aa_backconnect['SER']=['N','CA','CA','C','C','O','C','CA']
aa_backconnect['THR']=['N','CA','CA','C','C','O','C','CA']
aa_backconnect['CYS']=['N','CA','CA','C','C','O','C','CA']

# Definitions for drawing the amino acid sidechains.  Draw a red line from atom to atom to atom.
aa_connections['ARG'] =['CA','CB','CB','CG','CG','CD','CD','NE','NE','CZ','CZ','NH1','CZ','NH2']
aa_connections['LYS'] =['CA','CB','CB','CG','CG','CD','CD','CE','CE','NZ']
aa_connections['HIS'] =['CA','CB','CB','CG','CG','CD2','CD2','NE2','NE2','CE1','CE1','ND1','ND1','CG']
aa_connections['GLN'] =['CA','CB','CB','CG','CG','CD','CD','OE1','CD','NE2']
aa_connections['ASN'] =['CA','CB','CB','CG','CG','OD1','CG','ND2']
aa_connections['GLU'] =['CA','CB','CB','CG','CG','CD','CD','OE1','CD','OE2']
aa_connections['ASP'] =['CA','CB','CB','CG','CG','OD2']
aa_connections['TRP'] =['CA','CB','CB','CG','CG','CD1','CD1','NE1','NE1','CE2','CE2','CD2','CD2','CG','CD2','CE3','CE3','CZ3','CZ3','CH2','CH2','CZ2','CZ2','CE2']
aa_connections['TYR'] =['CA','CB','CB','CG','CG','CD1','CD1','CE1','CE1','CZ','CZ','OH','CZ','CE2','CE2','CD2','CG']
aa_connections['PHE'] =['CA','CB','CB','CG','CG','CD1','CD1','CE1','CE1','CZ','CZ','CE2','CE2','CD2','CD2','CG']
aa_connections['PRO'] =['CA','CB','CB','CG','CG','CD','CD','N']
aa_connections['MET'] =['CA','CB','CB','CG','CG','SD','SD','CE']
aa_connections['ILE'] =['CA','CB','CG2','CB','CB','CG1','CG1','CD1'] #Switch CB-CG2 instead of CG1-CG2
aa_connections['LEU'] =['CA','CB','CB','CG','CG','CD2','CG','CD1'] #Switch CD2 with CG because of unnessary bond between CD2 and CD1
aa_connections['ALA'] =['CA','CB']
aa_connections['VAL'] =['CA','CB','CB','CG1','CB','CG2']
aa_connections['GLY'] =['N','CA','CA','C']
aa_connections['SER'] =['CA','CB','CB','OG']
aa_connections['THR'] =['CA','CB','OG1','CB','CB','CG2']
aa_connections['CYS'] =['CA','CB','CB','SG']

# Definitions for drawing the amino acid hydrogens.  Each tuple defines a line to be drawn.
aa_hydrogen_connections['ALA']=[('CA','HA'),('CB','HB1'),('CB','HB2'),('CB','HB3')]
aa_hydrogen_connections['ARG']=[('N','H'),('N','H2'),('CA','HA'),('CB','HB2'),('CB','HB3'),('CG','HG2'),('CG','HG3'),('CD','HD2'),('CD','HD3'),('NE','HE'),('NH1','HH11'),('NH1','HH12'),('NH2','HH21'),('NH2','HH22')]
aa_hydrogen_connections['ASN']=[('CA','HA'),('CB','HB3'),('CB','HB2'),('ND2','HD21'),('ND2','HD22')]
aa_hydrogen_connections['ASP']=[('CA','HA'),('CB','HB2'),('CB','HB3'),('HD2','OD2')]
aa_hydrogen_connections['CYS']=[('CA','HA'),('CB','HB2'),('CB','HB3'),('SG','HG')]
aa_hydrogen_connections['GLU']=[('CA','HA'),('CB','HB2'),('CB','HB3'),('CG','HG2'),('CG','HG3')]
aa_hydrogen_connections['GLY']=[('CA','HA3'),('CA','HA2')]
aa_hydrogen_connections['HIS']=[('CA','HA'),('CB','HB2'),('CB','HB3'),('CD2','HD2'),('NE2','HE2'),('CE1','HE1'),('ND1','HD1')]
aa_hydrogen_connections['ILE']=[('CA','HA'),('CB','HB'),('CG1','HG12'),('CG1','HG13'),('CG2','HG21'),('CG2','HG22'),('CG2','HG23'),('CD1','HD11'),('CD1','HD12'),('CD1','HD13')]
aa_hydrogen_connections['LEU']=[('CA','HA'),('CB','HB2'),('CB','HB3'), ('CG','HG'),('CD2','HD21'),('CD2','HD22'),('CD2','HD23'),('CD1','HD11'),('CD1','HD12'),('CD1','HD13')]
aa_hydrogen_connections['LYS']=[('NZ','HZ1'),('NZ','HZ2'),('NZ','HZ3'),('CE','HE2'),('CE','HE3'),('HD2','CD'),('CD','HD3'),('CG','HG2'),('CG','HG3'),('CB','HB2'),('CB','HB3'),('CA','HA'),('N','H'),('OXT','HXT')]
aa_hydrogen_connections['MET']=[('CA','HA'),('CB','HB2'),('CB','HB3'),('CG','HG2'),('CG','HG3'),('CE','HE1'),('CE','HE2'),('CE','HE3')]
aa_hydrogen_connections['PHE']=[('CA','HA'),('CB','HB2'),('CB','HB3'),('CD1','HD1'),('CD2','HD2'),('CE1','HE1'),('CE2','HE2'),('CZ','HZ')]
aa_hydrogen_connections['PRO']=[('CA','HA'), ('N','H'),('CD','HD2'),('CD','HD3'),('CG','HG2'),('CG','HG3'),('CB','HB2'),('CB','HB3')]
aa_hydrogen_connections['SER']=[('CA','HA'),('CB','HB2'),('CB','HB3'),('OG','HG')]
aa_hydrogen_connections['THR']=[('CA','HA'),('CB','HB'),('CG2','HG21'),('CG2','HG22'),('CG2','HG23'),('OG1','HG1')]
aa_hydrogen_connections['TRP']=[('CA','HA'),('CB','HB2'),('CB','HB3'),('CD1','HD1'),('NE1','HE1'),('CZ2','HZ2'),('CH2','HH2'),('CZ3','HZ3'),('CE3','HE3')]
aa_hydrogen_connections['TYR']=[('CA','HA'),('CB','HB2'),('CB','HB3'),('CD2','HD2'),('CE2','HE2'),('OH','HH'),('CE1','HE1'),('CD1','HD1')]
aa_hydrogen_connections['VAL']=[('CA','HA'),('CB','HB'),('CG2','HG21'),('CG2','HG22'),('CG2','HG23'),('CG1','HG11'),('CG1','HG12'),('CG1','HG13')]

