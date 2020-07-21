import operator
import warnings
import numpy as np
import pandas as pd
import Bio
import Bio.SeqIO
from Bio import pairwise2
from Bio.PDB import PDBParser
from biograph.structure import StructureModel
from biograph.downloader import PdbDownloader
from biograph.constants import amino_1code, valid_amino_3, valid_amino_1
from biograph import alignment



class Protein:
    """Main BioGraph class for representing proteins through dataframes.
    Provides utilities for handling atoms, constructing structures or graphs
    and target variables."""
    # Bio often produces a couple of warnings when loading pdb files.
    # This can cloud logs when dealing with a large amount of pdbs.
    # We provide the option to suppress those warnings, and with
    # this module-level flag we make sure to only warn about it once.
    _warned_about_suppressing_bio = False
    @staticmethod
    def fetch(pdb_id, base_path=".", suppress_bio_warnings=True):
        """Fetch a PDB file and instantiate a Protein with it."""
        dw = PdbDownloader([pdb_id], base_path = base_path)
        filenames = dw.request_and_write()
        if filenames[0] is not None:
            return Protein(filenames[0], suppress_bio_warnings=suppress_bio_warnings)
        raise Exception("PDB could not be downloaded")

    def __init__(self, pdb, suppress_bio_warnings=True):
        self.suppress_bio_warnings = suppress_bio_warnings
        self.pdb = pdb
        self.structure = None
        self._df = None
        self._seq = None

    @property
    def df(self):
        if self._df is None:
            self.generate_dataframe()
        return self._df

    @df.setter
    def df(self, df):
        self._df = df

    @property
    def sequences(self):
        """Tries to load the sequences from the PDB file if present.
        It is also possible to use a PPBuilder for this, but sometimes
        it breaks up a sequence into smaller parts and they stop
        matching with the FASTA files.
        Returns a dictionary of the form chain_id => sequence"""
        if self._seq is not None:
            return self._seq

        if self.pdb_file is None:
            return None

        self._seq = {}
        with open(self.pdb_file, "rU") as handle:
            for record in Bio.SeqIO.parse(handle, "pdb-seqres"):
                # record.name often is <unknown value>, but ID is not
                # so we define the name as pdb id + chain id
                chain_id = record.id #sometimes it's "A" (3RRF) sometimes it's "109T:A"
                if ":" in chain_id:
                    chain_id = chain_id.split(":")[1]
                name = "{}_{}".format(self.pdb.id, chain_id)
                self._seq.update({name:record.seq})

        return self._seq

    @property
    def pdb(self):
        return self.__pdb

    @pdb.setter
    def pdb(self, pdb):
        self.pdb_file = None
        if isinstance(pdb, Bio.PDB.Structure.Structure):
            self.__pdb = pdb
        elif isinstance(pdb, str):
            self.pdb_file = pdb
            with warnings.catch_warnings():
                if self.suppress_bio_warnings:
                    if not Protein._warned_about_suppressing_bio:
                        warnings.warn("suppress_bio_warnings=True, ignoring Bio's warnings.")
                        Protein._warned_about_suppressing_bio = True
                    warnings.simplefilter("ignore")
                parser = PDBParser()
                # Infer pdb_id from filename
                pdb_id = pdb.split("/")[-1][:-4]
                pdb = parser.get_structure(pdb_id, pdb)
                self.__pdb = pdb
        else:
            raise Exception("""A Bio.PDB.Structure.Structure or a
                            valid path must be used""")

    @pdb.deleter
    def pdb(self):
        del self.__pdb

    def generate_structure(self, filter_rows):
        """
        Generate a structure model from selected rows of the dataframe.

        Parameters
        ----------
        filter_rows: function
            A filter function that is applied to each row of the dataframe

        Returns
        -------
        structure: StructureModel
        """
        rows = self.df.loc[
            self.df.apply(filter_rows, axis=1),
            ["full_id", "coord"]].reset_index(drop=True)
        ids, coords = rows["full_id"], rows["coord"]
        self.structure = StructureModel(ids, coords)
        return self.structure

    def generate_graph(self, model, model_params):
        self.graph = model.generate_graph(self, model_params)
        return self.graph

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

    def get_residues_dict(self, filter_res=lambda x: True,
                     filter_attr=lambda x: x):
        """Get residues data from PDB as an array of dictionaries.
        This method uses Bio.PDB.Structure.Structure.get_residues() method to
        iterate through residues. In addition, filters can
        be applied to choose which residues and which attributes to retrieve.

        ----------
        Parameters
        filter_res : function
            Function applied to each residue, used to filter attributes.
        filter_attr : function
            Function applied to each residue, it must return True or False.
            If True the residue is returned, if False, the residue is filtered.

        Returns
        -------
        numpy.Array
            Array of residue dicts.

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
        atoms = [filter_attr(res.__dict__) for res in
                 self.pdb.get_residues() if filter_res(res.__dict__)]
        return atoms

    def get_residues(self, filter_res=lambda x: True,
                     filter_attr=lambda x: x):
        """Get residues data from PDB as an array of Bio.PDB.Residue.Residue.
        This method uses Bio.PDB.Structure.Structure.get_residues() method to
        iterate through residues.  In addition, filters can
        be applied to choose which residues and which attributes to retrieve.

        ----------
        Parameters
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
        hoh = prot.get_residues(filter_res=lambda x: x.resname == "HOH")
        hoh[0].__dict__

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
        resnames = prot.get_residues(filter_attr=lambda x: x.resname)
        """
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

    def get_conservation_features(self, path, chain_list = None):
        """
        Calculates conservation features from aligning the sequence in the consurf file
        located in `path`. Features are those specified in alignment.join_conservation_data.
        Returns a dict that maps this protein's residue full ids to conservation features.

        Parameters
        ----------
        path: string
            Full path to consurf grades file.
        chain_list: list or None
            List of chains to be considered for alignment, or None if you want to align all of
            them. Each chain is aligned separately.
        Returns
        -------
        all_features: dict
            Maps res_full_ids to features. Format is taken from alignment.join_conservation_data
        """
        all_features = {}
        if chain_list is None:
            chain_list = set([r.parent.id for r in self.get_residues()])
        for chain in chain_list:
            valid_residues = [r for r in self.get_residues()
                if valid_amino_3(r.resname) and r.parent.id == chain]
            self_sequence = ''.join([amino_1code(r.resname) for r in valid_residues])
            features = {r.full_id:dict() for r in valid_residues}
            # Join feature through alignment of sequences.
            features = alignment.join_conservation_data(self_sequence, features, path)
            all_features.update(features)
        return all_features

    def add_residue_features(self, features):
        """
        Adds residue-level features to the dataframe.

        Parameters
        ----------
        features: dict
            Dictionary mapping residue full ids to new features.

        Returns
        -------
        None
        """
        self.df = pd.concat(
            [self.df, self.df.apply(
                lambda row: features[row["res_full_id"]] if row["res_full_id"] in features else None,
                axis = 1, result_type="expand")],
            axis=1)

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
        return self.df.head(n)

    @staticmethod
    def __get_coordinates(coord, raise_error=True):
        if hasattr(coord, "__getitem__"):
            return coord[0], coord[1], coord[2]
        else:
            if raise_error:
                raise TypeError("""Non-suscriptable type when parsing file.
                                Please, check if there are null values in PDB
                                or use raise_error=False""")
            else:
                return None, None, None

    def generate_dataframe(self,
                           columns=["bfactor", "chain", "coord",
                                    "disordered_flag", "element",
                                    "full_id", "res_full_id",
                                    "mass", "resname", "occupancy"],
                           split_coordinates=True, raise_error=False):
        """ Generate a Pandas DataFrame from the PDB
        Parameters
        ----------
        columns : list
            list of column names to subset DataFrame
        split_coordinates: bool
            whether to return three extra columns x, y, z for coordinates or not
        raise_error: bool
            when trying to split coordinates you can choose to raise error or
            to generate null values. Default is False.

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
                        atom_i["res_full_id"] = residue.full_id
                        atom_i["chain"] = chain.id
                        atom_i["model"] = str(model.id)
                        full_atom.append(atom_i)
                        del(atom_i)
        df = pd.DataFrame(full_atom).loc[:, columns]
        if split_coordinates:
            df["x"], df["y"], df["z"] = zip(*df.coord.apply(
                lambda x: self.__get_coordinates(x, raise_error)))
        self._df = df

    def select_chains(self, chain_list):
        """Discards rows from the dataframe that are not in the chainlist."""
        self.df = self.df[self.df.chain.isin(chain_list)]

    def _filter_het_rows(self, allowed_ligands, discard_water, keep_normal_atoms):
        """Provides a pandas filter to handle ligands, water and normal atoms"""
        return lambda res: ((keep_normal_atoms and (res[3][0][0:2] not in ["H_", "W"]))
                            or (res[3][0][0:2] == "H_" and res[3][0][2:] in allowed_ligands)
                            or (res[3][0] == "W" and not discard_water))

    def select_ligands(self, allowed_ligands, discard_water=True):
        """Keep only heterogen atoms that are ligands in provided list.
        `allowed_ligands` must be hetIDs or ligand IDs as found in the PDB,
        e.g. ATP, BFS, GOL, etc.
        `discard_water` controls whether to filter out water atoms.
        """
        if isinstance(allowed_ligands, str):
            allowed_ligands = [allowed_ligands]

        allowed_rows = self.df.res_full_id.apply(
            self._filter_het_rows(allowed_ligands, discard_water, True))

        self.df = self.df.loc[allowed_rows]

    def extract_ligands(self, allowed_ligands, remove_from_df=True):
        """Extract and return ligand rows from dataframe, optionally
        removing from the dataframe as well (`remove_from_df`)"""
        if isinstance(allowed_ligands, str):
            allowed_ligands = [allowed_ligands]

        is_ligand = self.df.res_full_id.apply(
            self._filter_het_rows(allowed_ligands, True, False))

        ligand_rows = self.df.loc[is_ligand]
        if remove_from_df:
            self.df = self.df.loc[~is_ligand]

        return ligand_rows

    def discard_ligands(self):
        """Discards rows from the dataframe that correspond to ligands or
        heterogen atoms.
        For more information read about the HET section in PDB files:
        https://www.wwpdb.org/documentation/file-format-content/format33/sect4.html"""
        het_rows = self.df.res_full_id.apply(
            lambda res: res[3][0] == "W" or res[3][0][0:2] == "H_")
        self.df = self.df.loc[~het_rows]

    def discard_empty_coordinates(self):
        """Discards rows from the dataframe with missing coordinates.
        This is necessary when building structure models and the like.
        More info on missing coords:
        pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/missing-coordinates-and-biological-assemblies"""
        self.df = self.df.loc[~self.df.coord.isnull()]

    def add_distance_to_target_feature(self, target_rows):
        """Adds a `distance` feature to the dataframe which holds
        the minimum distance of each atom to the atoms of target_rows.
        This can be useful as a target for supervised learning."""
        target_coords = target_rows.coord.to_list()

        self.df["distance"] = self.df.coord.apply(
            lambda atom: min(map(lambda target_atom: np.linalg.norm(atom-target_atom), target_coords))
        )