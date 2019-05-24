import numpy as np
from Bio import pairwise2
import operator
from pyprot.structure import StructureModel
import Bio
import pandas as pd


class Protein:
    def __init__(self, pdb, generate_dataframe=True):
        self.pdb = pdb
        self.aminoacids = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
                           'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                           'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
                           'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
                           'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}
        #self.atoms = self._get_atoms()
        self.structure = None
        if generate_dataframe:
            self.df = self.generate_dataframe()
        else:
            self.df = None

    @property
    def pdb(self):
        return self.__pdb

    @pdb.setter
    def pdb(self, pdb):
        assert isinstance(pdb, Bio.PDB.Structure.Structure), """A Bio.PDB.Structure.Structure must be used"""
        self.__pdb = pdb

    @pdb.deleter
    def pdb(self):
        del self.__pdb

    def generate_structure(self, points):
        self.structure = StructureModel(points)

    def get_atoms(self, atom_as_dict=True, filter_atoms=lambda x: True,
                  filter_attr=lambda x: x):
        """Get atoms data from PDB.
        This method uses Bio.PDB.Structure.Structure.get_atoms() method to
        iterate through atoms. Atom can be represented both as a
        Bio.PDB.Atom.Atom object or as a dict. In addition, filters can be
        applied to choose what atoms and what atom's attributes retrieve.

        Parameters
        ----------
        atom_as_dict : bool
            Whether to represent an atom as a dict or as a Bio.PDB.Atom.Atom
            Default True
        filter_atoms : function
            Function applied to each atom, used to filter attributes.
        filter_attr : function
            Function applied to each atom, it must return True or False.
            If True the atom is returned, if False, the atom is filtered.

        Returns
        -------
        numpy.Array
            Array of atoms.

        Examples
        -------
        # Get CA atoms using dict representation
        ca_atoms = prot.get_atoms(filter_atoms=lambda x: x["name"] == "CA")

        # Get CA atoms using Bio representation
        ca_atoms = prot.get_atoms(
            filter_atoms=lambda x: x.name == "CA", atom_as_dict)

        # Just retreieve coordinates
        ca_atoms = prot.get_atoms(filter_attr=lambda x: x["coord"])
        """
        if atom_as_dict:
            atoms = [filter_attr(atom.__dict__) for atom in
                     self.pdb.get_atoms() if filter_atoms(atom.__dict__)]
        else:
            atoms = [filter_attr(atom) for atom in
                     self.pdb.get_atoms() if filter_atoms(atom)]
        return np.array(atoms)

    def get_residues(self, res_as_dict=True, filter_res=lambda x: True,
                     filter_attr=lambda x: x):
        """Get residues data from PDB.
        This method uses Bio.PDB.Structure.Structure.get_residues() method to
        iterate through residues. Residue can be represented both as a
        Bio.PDB.Residue.Residue object or as a dict. In addition, filters can
        be applied to choose what residue and what residue's attributes
        retrieve.

        ----------
        Parameters
        res_as_dict : bool
            Whether to represent a residue as a dict or as a
            Bio.PDB.Residue.Residue
            Default True
        filter_res : function
            Function applied to each residue, used to filter attributes.
        filter_attr : function
            Function applied to each residue, it must return True or False.
            If True the residue is returned, if False, the residue is filtered.

        Returns
        -------
        numpy.Array
            Array of residues.

        Examples
        -------
        # Get HOH using dict representation
        hoh = prot.get_residues(filter_res=lambda x: x["resname"] == "HOH")
        hoh[0]

        {'_id': ('W', 201, ' '),
         'child_dict': {'O': <Atom O>},
         'child_list': [<Atom O>],
         'disordered': 0,
         'full_id': ('1a3z', 0, 'A', ('W', 201, ' ')),
         'level': 'R',
         'parent': <Chain id=A>,
         'resname': 'HOH',
         'segid': '    ',
         'xtra': {}}

        # Get residues names
        resnames = prot.get_residues(filter_attr=lambda x: x["resname"])
        """
        if res_as_dict:
            atoms = [filter_attr(res.__dict__) for res in
                     self.pdb.get_residues() if filter_res(res.__dict__)]
        else:
            atoms = [filter_attr(res) for res in
                     self.pdb.get_residues() if filter_res(res)]
        return atoms

    def _get_bfactor_by_atom(self):
        """
        Get bfactor data
        :return: list of list of dicts, where each dict represents information of an atom and the list
        represents a residue
        """
        bfactor = [[{"enum": e[0], "bfactor": e[1].bfactor, "atom_full_id": e[1].full_id,
                     "res_full_id": e[1].parent.full_id} for e in enumerate(i.get_atoms())] for i in
                   self.pdb.get_residues()]
        return bfactor

    def _get_avg_bfactor_by_residue(self):
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
        seqSelf = ''.join([self.aminoacids[r.resname] for r in self._get_residues() if r.resname in
                           self.aminoacids.keys() and r.parent.id in chain_list])
        seqCS = ''.join([y[1].lstrip() for y in thelines])
        alignment = max(pairwise2.align.globalxx(seqSelf, seqCS), key=operator.itemgetter(1))
        seqSelf, seqCS, _, _, _ = alignment
        enumres = enumerate([r for r in self._get_residues() if r.resname in
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

    def head(self, n=10):
        """  Return the first `n` rows (atoms) of pandas.DataFrame representation
        of PDB.


        Parameters
        ----------
        n : int
            Number of rows (atoms) to select.

        Returns
        -------
        pandas.Dataframe
            n first rows of pandas.DataFrame representation of PDB

        """
        if self.df is not None:
            return self.df.head(n)
        else:
            self.df = self.generate_dataframe()
            return self.df.head(n)

    def generate_dataframe(self, columns=["bfactor", "chain", "coord",
                                          "disordered_flag", "element",
                                          "full_id", "mass", "resname",
                                          "occupancy"]):
        """ Generate a Pandas DataFrame from the PDB
        Parameters
        ----------
        columns : list
            list of column names to subset DataFrame

        Returns
        -------
        Pandas.DataFrame
            DataFrame which has an atom per row
        """
        full_atom = []
        for model in self.pdb:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        atom_i = atom.__dict__
                        atom_i["resname"] = residue.resname
                        atom_i["chain"] = chain.id
                        atom_i["model"] = str(model.id)
                        full_atom.append(atom_i)
                        del(atom_i)
        return pd.DataFrame(full_atom).loc[:, columns]
