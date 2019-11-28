import json
import pandas as pd
import os
""" Module for reading different kind of files"""


class Writer:
    """Class used for writing different things easily"""
    def __init__(self, path):
        """

        :param path: path can be a file path or a directory path
        """
        self.path = path

    def write_str(self, string_obj):

        with open(self.path, "w") as f:
            f.write(string_obj)

    def create_directory(self, use_dirname=False):
        path = self.path
        if use_dirname:
            path = os.path.dirname(self.path)
        if not os.path.exists(path):
            os.makedirs(path)


class Reader:
    def __init__(self, path):
        self.path = path

    def from_json(self, to_df = False):
        with open(self.path) as f:
            obj_from_file = json.load(f)
            if to_df:
                obj_from_file = pd.DataFrame(obj_from_file)
        return obj_from_file
