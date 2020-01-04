from biocoso.io import Reader
from flatten_dict import flatten
import pandas as pd
import json


class SiftsProcessor:
    def __init__(self, path):
        """This class takes the raw Sifts files and implements a method that returns a Pandas DataFrame"""
        self.json_file = Reader(path).from_json()

    def process_dict(self, db):
        """
        This method iterates through the .json file, get all the relevant information and structurate it
        :param db: list of databases wanted to be returned
        :return: pandas DataFrame
        """
        output_dict = {}
        result = []
        for key in self.json_file.keys():
            pdb_id = key
            output_dict["pdb_id"] = pdb_id
            for i in self.json_file[pdb_id].keys():
                if i in db:
                    for j in self.json_file [pdb_id][i].keys():
                        df = self.process_subdict_(self.json_file [pdb_id][i][j])
                        df["db"] = i
                        df["pdb"] = pdb_id
                        df["protein_id"] = j
                        df = pd.melt(df, id_vars=["pdb", "db", "protein_id", "description", "identifier", "chain_id"])
                        result.append(df)
        return pd.concat(result).reset_index(drop=True)

    def process_subdict_(self, entry_dict):
        """
        Helper method used by process_dict
        :param entry_dict:
        :return:
        """
        try:
            description = entry_dict["description"]
        except Exception:
            description = ""
        try:
            identifier = entry_dict["identifier"]
        except Exception:
            identifier = ""
        mappings = entry_dict["mappings"]
        new_mapping = []
        for chain in mappings:
            flattened_dict = flatten(chain, reducer=self.underscore_reducer)
            new_mapping.append(flattened_dict)

        df = pd.DataFrame(new_mapping)
        df["description"] = description
        df["identifier"] = identifier

        return df

    @staticmethod
    def underscore_reducer(k1, k2):
        """Static function used by flatten_dict.flatten"""
        if k1 is None:
            result = k2
        else:
            result = k1 + "_" + k2
        return result


class ConsurfDB:

    @staticmethod
    def parse_pdbaa_list(input_path, write_json= False, output_path='consurf_equivalence_dict.json'):
        """
        This method reads the pdbaa_list.nr file and returns a processed json ready to be used to request pdbs from
        ConsurfDB. It adds each key as a value of itself, which makes it easier to make requests to ConsurfDB.
        :param input_path: str, path where .nr file is stored
        :param write_json: bool, if True then write .json file
        :param output_path: str, path where .json file will be stored
        :return: dict, .nr file processed
        """
        data = {}

        with open(input_path) as ifile:
            for line in ifile.readlines():
                h, l = line.rstrip('\n').rstrip('.').split(':')
                l = l.replace(' ', '').split(',')
                l = [x for x in l if len(x) > 0]
                l.append(h)
                data[h] = l
        if write_json:
            with open(output_path, 'w') as ofile:
                json.dump(data, ofile)
        return data