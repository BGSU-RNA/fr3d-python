RNAbaseheavyatoms = {}
RNAbasehydrogens = {}
RNAconnections = {}
Ribophos_connect = {}
tilt_cutoff = {}
Normal_residue = {}
aa_connections = {}
aa_backconnect = {}

#Hard-coding for Maryam's code
ChainNames = {}

ChainNames ['5AJ3'] = {'B' :'uS2m' , 'C': 'uS3m/ uS24m', 'E':'uS5m', 'F':'bS6m',
                       'G':'uS7m', 'I':'uS9m', 'J':'uS10m', 'K':'uS11m', 'L':'uS12m',
                       'N':'uS14m', 'O':'uS15m','P':'bS16m','Q':'uS17m',
                       'R':'uS18m/ uS18c', 'T':'bL19m', 'U':'bS21m','a':'mS22',
                       'b':'mS23', 'c':'mS25', 'd': 'mS26', 'e':'mS27', 'f':'bS1m/mS28',
                       'g':'mS29','h':'mS31','i':'mS33','j':'mS34','k':'mS35','m':'mS37',
                       'n':'mS38','o':'mS39','p':'mS40'}


ChainNames ['3J9M'] = {'B' :'uS2m' , 'C': 'uS3m/ uS24m', 'D':'uS5m',
                       'E':'bS6m', 'F':'uS7m', 'G':'uS9m', 'H':'uS10m',
                       'I':'uS11m', 'J':'uS12m', 'K':'uS14m', 'L':'uS15m',
                       'M':'bS16m', 'N':'uS17m', 'O':'bS18b/ mS40','P':'uS18m/ uS18c',
                       'Q':'bS21m', 'R':'mS22', 'S':'mS23', 'T':'mS25', 'U':'mS26',
                       'V': 'mS27', 'W':'bL27m','X':'mS29', 'Y':'mS31','Z':'mS33','0': 'bL32m',
                       '1': 'bL33m', '2': 'bL34m', '3': 'bL35m', '4': 'bL36m',
                       '5': 'mL37','6': 'mL38','7': 'mL39', '8': 'mL40','9': 'mL41',
                       'a':'mS34', 'b':'mS35', 'c':'mS37', 'd': 'mS38', 'e':'ms39', 'f':'mL48',
                       'g':'mL49','h':'mL50','i':'mL51','j':'mL52','k':'mL53','l':'mL54',
                       'm':'bL31m', 'o':'mL63','p':'mL62','q':'mL62', 'r':'mL66/bS18a',
                       's':'mL65/ mS30', 'AB':'uS2m','AC':'uS3m', 'AD':'uS5m',
                       'AE':'bS6m', 'AF':'uS7m','AG':'uS9m', 'AH':'uS10m','AI':'uS11m',
                       'AJ':'uS12m', 'AK':'uS14m', 'AL':'uS15m','AM':'bS16m','AN':'uS17m',
                       'AO':'mS40/ bS18b', 'AP':'bS18m/ bS18c','AQ':'bS21m','AR':'mS22',
                       'AS': 'mS23','AT': 'mS25','AU': 'mS26','AV': 'mS27','AW': 'bS21m',
                       'AX': 'mS29','AY': 'mS31','AZ': 'mS33','A0': 'mS34','A1': 'mS35',
                       'A2': 'mS37','A3': 'mS38','A4': 'mS39',}

ChainNames ['3I8G'] = {'B':'bS1', 'E' :'uS2', 'F' :'uS3', 'G' :'uS4', 'H' :'uS5', 'I' :'bS6',
                       'J' :'uS7', 'K' :'uS8', 'L' :'uS9', 'M' :'uS10', 'N' :'uS11',
                       'O' :'uS12', 'P' :'uS13', 'Q': 'uS14', 'R' :'uS15', 'S' :'bS16',
                       'T' :'bS17', 'U' :'uS18', 'V' :'uS19', 'W' :'bS20', 'X' :'bS21'}

ChainNames ['2AW7'] = {'B' :'uS2', 'C' :'uS3', 'D' :'uS4', 'E' :'uS5', 'F' :'bS6',
                       'G' :'uS7', 'H' :'uS8', 'I' :'uS9', 'J' :'uS10', 'K' :'uS11',
                       'L' :'uS12', 'M' :'uS13', 'N': 'uS14', 'O' :'uS15', 'P' :'bS16',
                       'Q' :'uS17', 'R' :'bS18', 'S' :'uS19', 'T' :'bS20', 'U' :'bS21'}
                       
#List of amino acids for perpendicular interactions

Perp_aa = set(['ARG','LYS','GLU','ASP','GLN','ASN','HIS','TYR','TRP','PHE'])

#Dictionaries for normal to plane calculations

Normal_residue['A'] = ['C4','C5','N3']
Normal_residue['U'] = ['C2','N1','N3']
Normal_residue['C'] = ['C6','N1','C5']
Normal_residue['G'] = ['C4','C5','N3']

Normal_residue['ARG'] =['CZ','NH1', 'NH2']
Normal_residue['LYS'] =['CE','CD','NZ']
#Normal_residue['HIS'] = ['ND1','CE1','CG']
Normal_residue['HIS'] = ['NE2','CD2','CE1']
Normal_residue['GLN'] =['CD','OE1','NE2']
Normal_residue['ASN'] =['CG','OD1','ND2']
Normal_residue['GLU'] =['CD','OE1','OE2']
Normal_residue['ASP'] =['CG','OD1','OD2']
Normal_residue['TRP'] =['CD2','CG','CE2']
Normal_residue['TYR'] =['CE2','CZ','CE1']
Normal_residue['PHE'] =['CG','CD1','CZ']
Normal_residue['PRO'] =['CG','CD','CB']
Normal_residue['MET'] =['SD','CG','CE']
Normal_residue['ILE'] =['CG1','CB','CD1']
Normal_residue['LEU'] =['CG','CD1','CD2']
#Normal_residue['ALA'] =[]
#Normal_residue['GLY'] =[]
Normal_residue['VAL'] =['CB','CG1','CG2']
Normal_residue['SER'] =['CB','OG','CA']
Normal_residue['THR'] =['CB','OG1','CG2']                       
Normal_residue['CYS'] =['CB','SG','CB']



#Creating dictionaries for detecting edges of nts

"""WC_1 = {}
WC_2 = {}
WC_3 = {}
Hoogsteen_1 = {}
Hoogsteen_2 = {}
Sugar = {}

WC_1['A']['HIS']= [6.000, 6.000, 5.500]
WC_2['A']['HIS']= [4.000, 4.000, 4.500]
WC_3['A']['HIS']= [6.000, 6.000, 6.500]

Hoogsteen_1['A']['HIS']= [-1.000, -1.000, -1.100]
Hoogsteen_2['A']['HIS']= [-3.000, -3.000, -3.000]

Sugar['A']['HIS']= [4.000, 2.000, 4.000]
"""
#Defining center-to-center and tilt cutoffs for stacking
#dist_cent_cutoff = {'TRP':7,'TYR': 7,'PHE':7, 'ALA':6,'ARG':6}
tilt_cutoff= {'ALA': 2,'VAL': 0.7,'ILE': 1.9,'LEU': 2.1,'ARG': 1.5,'LYS': 1.5,'HIS': 1.2,'ASP': 1.5,'GLU': 1.5,'ASN': 1.4,'GLN': 1.4,'THR': 0.5,'SER': 0.5,'TYR': 2.1,'TRP': 2.1,'PHE': 1.5,'PRO': 3.1,'CYS': 1.0, 'MET': 1.5}

#Definitions for drawing the amino acids

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
aa_connections['ILE'] =['CA','CB','CB','CG1','CG1','CG2','CG2','CD1']
aa_connections['LEU'] =['CA','CB','CB','CG','CG','CD2','CD2','CD1']
aa_connections['ALA'] =['CA','CB']
aa_connections['VAL'] =['CA','CB','CB','CG1','CB','CG2']
aa_connections['GLY'] =[]
aa_connections['SER'] =['CA','CB','CB','OG']
aa_connections['THR'] =['CA','CB','CB','OG1','CB','CG2']                       
aa_connections['CYS'] =['CA','CB','CB','SG']

#Definitions for drawing the RNA nucleotides

RNAconnections['A'] =['N1','C6','C6','N6','C6','C5','C5','C4','C5','N7','N7','C8','C8','N9','N9','C4','C4','C5','C4','N3','N3','C2','C2','N1']   
RNAconnections['U'] =['N1','C2','C2','O2','C2','N3','N3','C4','C4','O4','C4','C5','C5','C6','C6','N1']            
RNAconnections['G'] =['N1','C6','C6','O6','C6','C5','C5','C4','C5','N7','N7','C8','C8','N9','N9','C4','C4','C5','C4','N3','N3','C2','C2','N2','C2','N1']   
RNAconnections['C'] =['N1','C2','C2','O2','C2','N3','N3','C4','C4','N4','C4','C5','C5','C6','C6','N1']

Ribophos_connect['A'] = ["N9","C1'","C1'","C2'","C2'","O2'","C2'","C3'","C3'","O3'","C3'","C4'","C4'","O4'","O4'","C1'","O4'","C4'","C4'","C5'","C5'","O5'","O5'","P","P","OP1","P","OP2"]
Ribophos_connect['U'] = ["N1","C1'","C1'","C2'","C2'","O2'","C2'","C3'","C3'","O3'","C3'","C4'","C4'","O4'","O4'","C1'","O4'","C4'","C4'","C5'","C5'","O5'","O5'","P","P","OP1","P","OP2"]
Ribophos_connect['G'] = ["N9","C1'","C1'","C2'","C2'","O2'","C2'","C3'","C3'","O3'","C3'","C4'","C4'","O4'","O4'","C1'","O4'","C4'","C4'","C5'","C5'","O5'","O5'","P","P","OP1","P","OP2"]
Ribophos_connect['C'] = ["N1","C1'","C1'","C2'","C2'","O2'","C2'","C3'","C3'","O3'","C3'","C4'","C4'","O4'","O4'","C1'","O4'","C4'","C4'","C5'","C5'","O5'","O5'","P","P","OP1","P","OP2"]

#RNA computation definitions
"""Defining the parts of nt that we use to compute RNA-amino acid interactions"""

RNAbaseheavyatoms['A'] = ['N9','C4','N3','N1','C6','N6','C8','C5','C2','N7']
RNAbasehydrogens['A'] = ['H2','H8','H9','1H6','2H6']
RNAbaseheavyatoms['C'] = ['N1','C2','O2','N3','C4','N4','C6','C5']
RNAbasehydrogens['C'] = ['H1','H6','H5','1H4','2H4']
RNAbaseheavyatoms['G'] = ['N9','C4','N3','N1','C6','O6','C8','C5','C2','N7','N2']
RNAbasehydrogens['G'] = ['H1','H8','H9','1H2','2H2']
RNAbaseheavyatoms['U'] = ['N1','C2','O2','N3','C4','O4','C6','C5']
RNAbasehydrogens['U'] = ['H5','H1','H3','H6']
sugar_back = ["C1'","C2'","O2'","C3'","O3'","C4'","O4'","C5'","O5'","P","OP1","OP2"]

nt_backbone = {}

nt_backbone['A'] = sugar_back
nt_backbone['U'] = sugar_back
nt_backbone['C'] = sugar_back
nt_backbone['G'] = sugar_back

#Amino acid computation definitions

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
aa_linker['ASN'] = []
aa_fg['ASN'] = ['CB','CG','OD1','ND2']

aa_backbone['GLU'] = backbone
aa_linker['GLU'] = ['CB']
aa_fg['GLU'] = ['CG','CD','OE1','OE2']

aa_backbone['ASP'] = backbone
aa_linker['ASP'] = []
aa_fg['ASP'] = ['CB','CG','OD1','OD2']

aa_backbone['TRP'] = backbone
aa_linker['TRP'] = ['CB']
aa_fg['TRP'] = ['CG','CD1','NE1','CE2','CD2','CE3','CZ3','CH2','CZ2']

aa_backbone['TYR'] = backbone
aa_linker['TYR'] = ['CB']
aa_fg['TYR'] = ['CG','CD1','CE1','CZ','OH','CZ','CE2','CD2']

aa_backbone['PHE'] = backbone
aa_linker['PHE'] = ['CB']
aa_fg['PHE'] = ['CG','CD1','CD2','CE1','CZ','CE2']

aa_backbone['PRO'] = backbone
aa_linker['PRO'] = []
aa_fg['PRO'] = ['CA','CB','C','O','CG','CD']

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
aa_linker['THR'] = []
aa_fg['THR'] = ['CB','OG1','CG2']

aa_backbone['CYS'] = backbone
aa_linker['CYS'] = []
aa_fg['CYS'] = ['CB','SG']

# Quantum mechanics optimized base atom locations, by Jiri Sponer
# OK to change the listing of nucleotides to a Numpy array, if that helps
# The numbers are generated by zStandardBases.m; the formatting can be changed there
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
