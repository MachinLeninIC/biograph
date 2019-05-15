from Bio.PDB import PDBParser
from pyprot.protein import Protein
from pyprot.structure import StructureModel
from pyprot.structure import Perseus
import numpy as np

parser = PDBParser()
pdb = parser.get_structure("1a3z","/home/usuario/Documentos/github/pyprot/tests/1A3Z.pdb")
prot = Protein(pdb)
prot



ca_atoms = np.array([atom.coord for atom in prot.get_atoms_() if atom.name == "CA" and atom.full_id[2] == "A"])
struc = StructureModel(ca_atoms)
struc
per = Perseus()
sutr = per.calculate_fatcore(prot)

sutr.structure.persistent_hom_params
