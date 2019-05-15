import requests
from biocoso.io import Writer
import os
import re
from biocoso.url import *
import requests
import pandas as pd
from bs4 import BeautifulSoup


class StandardDownloader:

    def __init__(self, ids_list=None, url_object=None, base_path="download"):
        """
        This class implements method for generic downloading from url objects
        :param ids_list: list of ids to be downloaded
        :param url_object: object from biocoso.url
        :param base_path: path where files will be downloaded
        """
        self.ids_list = ids_list
        self.url_object = url_object
        self.base_path = base_path
        self.not_downloaded = []
        files_in_bp = (f for f in os.listdir(self.base_path) if os.path.isfile(os.path.join(self.base_path, f)))
        self.downloaded_in_path = [re.findall("(.*)?[.]", i)[0].strip() for i in files_in_bp]
        self.downloaded_in_path.sort()

    def request_and_write(self, create_dir_if_not_exists=False, use_dirname=False):
        for id in self.ids_list:
            try:
                if id not in self.downloaded_in_path:
                        sess = requests.Session()
                        request = requests.Request(method=self.url_object.method.upper(),
                                                   url=self.url_object.url.replace(self.url_object.mock_id, id))
                        response = sess.send(request.prepare())
                        if response.status_code == 200:
                            if create_dir_if_not_exists:
                                Writer(os.path.join(self.base_path, id + self.url_object.file_ext)).\
                                    create_directory(use_dirname=use_dirname)

                            Writer(os.path.join(self.base_path, id + self.url_object.file_ext)).write_str(response.text)
                            print("ID:", id, "succesfully written")
                        else:
                            self.not_downloaded.append(id)
                            print("ID:", id, "couldn't be written")
                            print("Status Code:", response.status_code)
                            print("Response:", response.text)
                else:
                    print("Previously downloaded id:", id)
            except Exception as e:
                print(e)
        Writer(os.path.join(self.base_path, "not_downloaded.txt")).write_str(str(self.not_downloaded))


class PdbDownloader:

    def __init__(self, ids_list=None, io_reader=None, pdb_url=PdbUrl("pdb_1"),
                 pdb_redo_url=PdbRedoUrl("pdbredo_3"), base_path="download"):
        """
        This class implements methods for downloading .pdb files
        Example:
        from biocoso.downloader import PdbDownloader
        ids_list = ['19HC', '1A05','1A0J','1A0M','1A12','1A1X','1A27','1A2Z','1A3A']
        dw = PdbDownloader(ids_list)
        dw.request_and_write()

        :param ids_list: list, PDB files ids to be downloaded
        :param io_reader: IOReader object that has a "to_ids_list" method which returns a list of ids
        :param pdb_url: url.PdbUrl object
        :param pdb_redo_url: url.PdbRedoUrl object
        :param base_path: path where files will be downloaded
        """
        #TODO: make IOReader object
        #TODO: migrate some functionality to abstract Downloader class
        self.pdb = pdb_url
        self.pdb_redo = pdb_redo_url
        self.base_path = base_path
        self.not_downloaded = []
        self.downloaded_in_path = [re.findall("(.*)?[.]",i)[0].strip() for i in os.listdir(self.base_path)]
        self.downloaded_in_path.sort()
        if io_reader is None:
            self.ids_list = ids_list
            self.type = None
        else:
            self.ids_list = io_reader.to_ids_list()
            self.type = io_reader.type
        self.url_to_complete = None

    def request_and_write(self, try_redo=True):
        """
        This method downloads each pdb from self.ids_list and save them to self.base_path. You can set whether to get
        files from PDB or from PDB-REDO, or first try to get it from PDB-REDO and, if it doesn't exists, then get it
        from PDB. This is its default behaviour.
        :param try_redo: bool, if true first try to download data from PDB-REDO
        :return: None
        """

        for id in self.ids_list:
            if id.strip() not in self.downloaded_in_path:
                try:
                    if try_redo:
                        response, file_ext = self.try_redo_then_pdb_(id)
                        if response.status_code == 200:
                            Writer(os.path.join(self.base_path, id + file_ext)).write_str(response.text)
                            print("ID:", id, "succesfully written")
                    else:
                        print("Entre")
                        new_url = self.pdb.url.replace(self.pdb.mock_id, id.strip().lower())
                        print(new_url)
                        response = requests.get(new_url)
                        if response.status_code == 200:
                            Writer(os.path.join(self.base_path, id + self.pdb.file_ext)).write_str(response.text)
                            print("ID:", id, "succesfully written")
                        else:
                            self.not_downloaded.append(id)
                            print("ID:", id, "couldn't be written")
                except Exception as e:
                    print(e)
            else:
                print("Previously downloaded id:", id)
        Writer(os.path.join(self.base_path, "not_downloaded.txt")).write_str(str(self.not_downloaded))

    def try_redo_then_pdb_(self, id):
        """
        Auxiliar method that first try to download the ID data from REDO and then from RCSB
        :param id: str, id of pdb file
        :return: request response, file extension
        """
        new_url = self.pdb_redo.url.replace(self.pdb_redo.mock_id,id.strip().lower())
        response = requests.get(new_url)
        file_ext = self.pdb_redo.file_ext
        if response.status_code != 200:
            new_url = self.pdb.url.replace(self.pdb.mock_id, id.strip().lower())
            response = requests.get(new_url)
            file_ext = self.pdb.file_ext
        return response, file_ext


class CapriDownloader:

    def __init__(self, capri_pdb_path = os.path.join("capri", "pdb"), capri_brk_path = os.path.join("capri", "brk"),
                 url_object = CapriUrl()):
        """
        Class used for downloading data from Capri dataset
        :param capri_pdb_path: path where capri pdbs will be downloaded
        :param capri_brk_path: path where capri brk files will be downloaded
        """
        self.pdb_path = capri_pdb_path
        self.brk_path = capri_brk_path
        self.capri_url = url_object

    def get_historical_data(self, drop_na=True):
        html = requests.get(self.capri_url.url)
        soup = BeautifulSoup(html.text, "html.parser")
        table = soup.find_all("tr")
        table_results = []
        for i in range(3, len(table) - 1):
            content_list = table[i].find_all("td")
            table_results.append(self.get_contents_(content_list))
        results = pd.DataFrame(table_results)
        if drop_na:
            results = results[~((results.pdb_url.isna()) | (results.starting_file_url.isna()))]
        return results

    @staticmethod
    def get_contents_(td_list):
        t_id = td_list[0].contents[0]
        try:
            protein_name = td_list[2].contents[0]
        except:
            protein_name = None
        try:
            final_pdb = td_list[4].contents[0]
            pdb_regex = re.findall('href="(.*)">?(.*)</a', str(final_pdb))
            pdb_url = pdb_regex[0][0]
            pdb_name = pdb_regex[0][1]
        except:
            pdb_url = None
            pdb_name = None
        try:
            starting_file = td_list[6].contents[0]
            sfile_regex = re.findall('href="(.*)">?(.*)</a', str(starting_file))
            sfile_url = sfile_regex[0][0].replace("./", "https://www.ebi.ac.uk/msd-srv/capri/")
            sfile_name = sfile_regex[0][1]
        except:
            sfile_url = None
            sfile_name = None
        return {"id": t_id, "protein_name": protein_name, "pdb_url": pdb_url, "pdb_name": pdb_name,
                "starting_file_url": sfile_url, "starting_file_name": sfile_name}


class ConsurfDBDownloader(StandardDownloader):

    def __init__(self, consurf_equivalence_dict, ids_list=None, consurf_url_object=ConsurfDBUrl(), base_path="ConsurfDB"):
        """
        Class used for downloading data from ConsurfDB
        :param ids_list: ids of pdb files and chains in the format PDB_CHAIN. For example
        :param consurf_url_object:
        :param pdbaa_list_path:
        """
        self.consurf_equivalence_dict = consurf_equivalence_dict
        StandardDownloader.__init__(self, ids_list=None, url_object=consurf_url_object, base_path=base_path)
        self.ids_list, self.not_matched, self.more_than_one_match = self.get_processed_ids_list_(ids_list)


    def get_processed_ids_list_(self, ids_list):
        new_list = []
        not_matched = []
        more_than_one_match = []
        for pdb_chain in ids_list:
            pdb_ref = [k for k in self.consurf_equivalence_dict.keys() if pdb_chain in self.consurf_equivalence_dict[k]]
            if not pdb_ref:
                print("PDB id wasn't matched")
                not_matched.append(pdb_chain)
            elif len(pdb_ref) > 1:
                print("There is more than one match, using the first one:", pdb_ref)
                more_than_one_match.append(pdb_chain)
                pdb_ref = pdb_ref[0].replace("_","/")
                new_list.append(pdb_ref)
                new_list = list(set(new_list))
            else:
                pdb_ref = pdb_ref[0].replace("_", "/")
                new_list.append(pdb_ref)
                new_list = list(set(new_list))
        return new_list, not_matched, more_than_one_match
