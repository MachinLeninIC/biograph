import numpy as np
from Bio import pairwise2
import operator
from pyprot.structure import StructureModel


class Protein:
    def __init__(self, pdb):
        self.pdb = pdb
        self.model = None
        self.aminoacids = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H',
                      'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q',
                      'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}
        self.atoms = self.get_atoms_()
        ca_atoms = np.array([atom.coord for atom in self.get_atoms_()
                             if atom.name == "CA" and atom.full_id[2] == "A"])
        self.structure = StructureModel(ca_atoms)


    def get_atoms_(self):
        """
        Get atoms from PDB
        :return: list of atoms
        """
        atoms = [i for i in self.pdb.get_atoms()]
        return atoms

    def get_residues_(self):
        """
        Get residues from PDB
        :return: list of residues
        """
        residues = [i for i in self.pdb.get_residues()]
        return residues

    def get_bfactor_by_atom_(self):
        """
        Get bfactor data
        :return: list of list of dicts, where each dict represents information of an atom and the list
        represents a residue
        """
        bfactor = [[{"enum": e[0], "bfactor": e[1].bfactor, "atom_full_id": e[1].full_id,
                     "res_full_id": e[1].parent.full_id} for e in enumerate(i.get_atoms())] for i in
                   self.pdb.get_residues()]
        return bfactor

    def get_avg_bfactor_by_residue_(self):
        """
        Get average bfactor by residue
        :return: list of dicts where each dict represents a residue
        """
        bfactor = [{"res_full_id": i.full_id, "avg_bfactor": np.mean([a.bfactor for a in i])} for i in
                   self.pdb.get_residues()]
        return bfactor

    def get_points_(self, dtype=np.float32):
        """
        Get coordinates of nodes (atoms)
        :param dtype: dtype of array to be returned, default is dtype32
        :return: np.array with coordinates of atoms
        """
        return np.array([atom.coord for atom in self.pdb.get_atoms() if atom["atom_full_id"][4][0] == "CA"], dtype=dtype)

    def read_conservation(self, path, chain_list):
        with open(path, 'r') as ifile:
            thelines = [x.rstrip('\n') for x in ifile.readlines()]
            thelines = [x.split('\t') for x in thelines]
            thelines = [x for x in thelines if len(x) == 14]
        if not thelines:
            return None
        seqSelf = ''.join([self.aminoacids[r.resname] for r in self.get_residues_() if r.resname in
                           self.aminoacids.keys() and r.parent.id in chain_list])
        seqCS = ''.join([y[1].lstrip() for y in thelines])
        alignment = max(pairwise2.align.globalxx(seqSelf, seqCS), key=operator.itemgetter(1))
        seqSelf, seqCS, _, _, _ = alignment
        enumres = enumerate([r for r in self.get_residues_() if r.resname in
                             self.aminoacids.keys() and r.parent.id in chain_list])
        thelines.reverse() # to pop() first element first
        prot_conservation = []
        count_good = 0
        count_bad = 0
        for x,y in zip(seqSelf, seqCS):
            res_conservation = {}
            i, res = next(enumres) if x != '-' else (-1, None)
            line = thelines.pop() if y != '-' else None
            if res and line:
                res_conservation["res_full_id"] = res.full_id
                res_conservation["score"] = float(line[3])
                res_conservation["color"] = line[5]
                res_conservation["score_confidence_interval"] = line[6]
                res_conservation["color_confidence_interval"] = line[9]
                res_conservation["residue_variety"] = line[-1]
                count_good += 1
                prot_conservation.append(res_conservation)
                if res.resname in self.aminoacids.values():
                    if self.aminoacids[res.resname] != x or x != y or y != line[1].lstrip():
                        print((x, y, i, res.resname, line[:3]))  # should never happen
                        count_bad += 1

        return prot_conservation, seqSelf, seqCS

    @staticmethod
    def distance_bet_res(r1, r2):
        atoms1 = np.array([(i.coord[0], i.coord[1], i.coord[2]) for i in r1.get_atoms()])
        atoms2 = np.array([(i.coord[0], i.coord[1], i.coord[2]) for i in r2.get_atoms()])
        v1 = np.repeat(atoms1.reshape(1,-1,3), atoms2.shape[0], axis=0)
        v2 = np.repeat(atoms2, atoms1.shape[0], axis=0).reshape(atoms2.shape[0],-1,3)
        return np.min(np.sum((v1 - v2)**2, axis=-1))
