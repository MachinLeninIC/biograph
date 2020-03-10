import os
import json
import configparser

def set_config(ini_path):
    cfg = configparser.ConfigParser()
    cfg.read(ini_path)
    global DATADIR
    global PROXY
    global PERSEUSPATH
    DATADIR = cfg["paths"]["DATADIR"]
    consurf = cfg["paths"]["CONSURF"]
    PROXY = cfg["proxy"]["PROXY"]
    PERSEUSPATH = cfg["paths"]["PERSEUSPATH"]

    with open(consurf, 'r') as f:
        global CONSURFCONFIG
        CONSURFCONFIG = json.load(f)

HOME = os.path.expanduser("~")
MODULEDIR = os.path.dirname(os.path.abspath(__file__))
# DATADIR = '/data/biocoso/out/'
# DATADIR = "/home/usuario/Documentos/protein_data/"
PERSEUSPATH = MODULEDIR + '/perseus/perseus'
PROXY = {'http': 'http://proxy.fcen.uba.ar:8080', 'https': 'http://proxy.fcen.uba.ar:8080'}
CONSURFURL = 'http://bental.tau.ac.il/new_ConSurfDB/DB/'
OUT = os.path.join(HOME, "bio", "data", "out")


SIDECHAIN = {'ALA': 'CA', 'ARG': 'CA NE', 'ASN': 'CA OD1', 'ASP': 'CA OD1', 'CYS': 'CA',
             'GLN': 'CA CD', 'GLU': 'CA CD', 'GLY': 'CA', 'HIS': 'CA CE1', 'ILE': 'CA CD1',
             'LEU': 'CA CD1', 'LYS': 'CA CE', 'MET': 'CA SD', 'PHE': 'CA CZ', 'PRO': 'CA',
             'SER': 'CA', 'THR': 'CA', 'TRP': 'CA CD2', 'TYR': 'CA CZ', 'VAL': 'CA'}



AMINOACIDS_3 = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
              'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

AMINOACIDS_3_TO_1 = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
    'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
    'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

FREQUENCIES = {'ALA':  0.0825, 'ARG':  0.0553, 'ASN':  0.0406, 'ASP':  0.0545, 'CYS':  0.0137, 'GLN':  0.0393,
              'GLU':  0.0675, 'GLY':  0.0707, 'HIS':  0.0227, 'ILE':  0.0596, 'LEU':  0.0966, 'LYS':  0.0584,
              'MET':  0.0242, 'PHE':  0.0386, 'PRO':  0.0470,'SER':  0.0656, 'THR':  0.0534, 'TRP':  0.0108,
              'TYR':  0.0292,  'VAL':  0.0687}

SURFACEFREQUENCIES = {'ALA': 0.0650510094915312,'ARG': 0.0603908302736197,'ASN': 0.0585967254645202,
                     'ASP': 0.0862229440008216,'CYS': 0.00621357229054869, 'GLN': 0.0493523249544252,
                     'GLU': 0.0997477341001874,'GLY': 0.0808620261723196,'HIS': 0.0244680805539151,
                     'ILE': 0.0303007933858834, 'LEU': 0.0524890663380149,'LYS': 0.0865631499216885,
                     'MET': 0.0114011776687978,'PHE': 0.0231350724488835,'PRO': 0.0560109464999444,
                     'SER': 0.0697582612267954,'THR': 0.0617858885151616,'TRP': 0.00979643275904862,
                     'TYR': 0.0263916347857345,'VAL': 0.0414623291481586}


HPIND = { 'ALA' : 1.8, 'CYS' :  2.5, 'ASP' : -3.5, 'GLU' : -3.5, 'PHE' :  2.8, 'GLY': -0.4, 'HIS':-3.2,
          'ILE' : 4.5, 'LYS' : -3.9, 'LEU' :  3.8, 'MET' : 1.9, 'ASN' : -3.5, 'PRO':-1.6, 'GLN':-3.5,
          'ARG' : -4.5, 'SER' : -0.8, 'THR' : -0.7, 'VAL' : 4.2, 'TRP' : -0.9, 'TYR': -1.3}
#          'UNK' : 0 }

# Helper functions for constants.

def amino_1code(amino_3code):
    return AMINOACIDS_3_TO_1[amino_3code]

def valid_amino_1(amino_1_code):
    return amino_1_code in AMINOACIDS_3_TO_1.keys()

def valid_amino_3(amino_3_code, unk_valid = False):
    return amino_3_code in AMINOACIDS_3 or (unk_valid and amino_3_code == "UNK")


