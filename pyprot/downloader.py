import requests
from pyprot.io import Writer
import os
import re
from pyprot.url import *
import requests
import pandas as pd

class StandardDownloader:

    def __init__(self, ids_list=None, url_object=None, base_path="download", force_download=False):
        """
        This class implements method for generic downloading from url objects
        Parameters
        ----------
        ids_list: list
            list of ids to be downloaded
        url_object: pyprot.url
        base_path: str
            path where files will be downloaded
        force_download: bool
            force download of files, even if they match a file
        """
        self.ids_list = ids_list
        self.url_object = url_object
        self.base_path = base_path
        self.force_download = force_download

        # Make a map of id => filename for ids that have already been downloaded.
        files_in_bp = (f for f in os.listdir(self.base_path)
            if os.path.isfile(os.path.join(self.base_path, f)))

        self.downloaded_files = {fn.split(".")[0].strip().lower() : os.path.join(self.base_path, fn)
            for fn in files_in_bp}

    def do_make_request(self, id):
        """
        Makes a specific id request.
        :param id: the id of the object being requested
        :return: (file_extension, file_content) or (None, None) on error
        """

        sess = requests.Session()
        request = requests.Request(method=self.url_object.method.upper(),
                                    url=self.url_object.url.replace(self.url_object.mock_id, id))
        response = sess.send(request.prepare())
        if response.status_code == 200:
            return self.url_object.file_ext, response.text

        return None, None

    def request_and_write(self, create_dir_if_not_exists = False, use_dirname = False, **kwargs):
        """Request IDs and write them to disk.
        This method downloads each id from self.ids_list and saves them to self.base_path.
        The request is made using self.do_make_request, and only the file extension and content
        that it returns are used.

        Parameters
        ----------
        create_dir_if_not_exists: bool
            (default: False)
        use_dirname : bool
            (default: False)
        **kwargs:
            additional keywords arguments to be passed on to do_make_request

        Returns
        -------
        List
            Each element has the filename of the corresponding element in the ids_list,
            or None if that id couldn't be downloaded.
        """
        files_list = []
        for id in self.ids_list:
            if id.lower() in self.downloaded_files.keys() and not self.force_download:
                print("Previously downloaded id:", id)
                files_list.append(self.downloaded_files[id.lower()])
                continue

            file_ext, file_text = self.do_make_request(id, **kwargs)
            if file_text is not None:
                filename = os.path.join(self.base_path, id + file_ext)
                if create_dir_if_not_exists:
                    Writer(os.path.join(filename)).\
                        create_directory(use_dirname=use_dirname)

                Writer(filename).write_str(file_text)
                files_list.append(filename)
            else:
                files_list.append(None)
                print("ID:", id, "couldn't be written")
                print("Status Code:", response.status_code)
                print("Response:", response.text)

        return files_list



class PdbDownloader(StandardDownloader):

    def __init__(self, ids_list=None, io_reader=None, pdb_url=PdbUrl("pdb_1"),
                 pdb_redo_url=PdbRedoUrl("pdbredo_3"), base_path="download"):
        """
        This class implements methods for downloading .pdb files
        Example:
        from pyprot.downloader import PdbDownloader
        ids_list = ['19HC', '1A05','1A0J','1A0M','1A12','1A1X','1A27','1A2Z','1A3A']
        dw = PdbDownloader(ids_list)
        dw.request_and_write()

        :param ids_list: list, PDB files ids to be downloaded
        :param io_reader: IOReader object that has a "to_ids_list" method which returns a list of ids
        :param pdb_url: url.PdbUrl object
        :param pdb_redo_url: url.PdbRedoUrl object
        :param base_path: path where files will be downloaded
        """
        super().__init__(ids_list=ids_list, base_path = base_path)
        self.type = None
        if io_reader is not None:
            self.ids_list = io_reader.to_ids_list()
            self.type = io_reader.type

        #TODO: make IOReader object

        self.pdb = pdb_url
        self.pdb_redo = pdb_redo_url
        self.url_to_complete = None

    def do_make_request(self, id, try_redo = True):
        """
        Download a pdb that matches the `id`.
        You can set whether to get it from PDB or from PDB-REDO, or first try to get it from
        PDB-REDO and, if it doesn't exists, then get it from PDB. This is its default behaviour.
        :param id: the id of the PDB being requested.
        :param try_redo: bool, if true first try to download data from PDB-REDO
        :return: (file_extension, file_content) or (None, None) on error.
        """

        if try_redo:
            response, file_ext = self.try_redo_then_pdb_(id)
            if response.status_code == 200:
                return file_ext, response.text
        else:
            new_url = self.pdb.url.replace(self.pdb.mock_id, id.strip().lower())
            response = requests.get(new_url)
            if response.status_code == 200:
                return self.pdb.file_ext, response.text

        return None, None

    def try_redo_then_pdb_(self, id):
        """
        Auxiliary method that first try to download the ID data from REDO and then from RCSB
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
