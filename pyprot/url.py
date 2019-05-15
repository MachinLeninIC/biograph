
class Url:
    def __init__(self, url_id=None):
        """
        Base class used for making different URLs
        :param url_id: id of endpoint to be used from URL
        """
        self.base_url = None
        self.url_id = url_id
        self.base_description = None
        self.urls = None

    @property
    def base_url(self):
        return self.__base_url

    @base_url.setter
    def base_url(self, url):
        self.__base_url = url

    @property
    def url_id(self):
        return self.__url_id

    @url_id.setter
    def url_id(self, url_id):
        self.__url_id = url_id

    @property
    def base_description(self):
        return self.__base_description

    @base_description.setter
    def base_description(self, descr):
        self.__base_description = descr

    def config_url(self, url_id):
        """
        Set values to attributes
        :param url_id: id of endpoint to be used from URL
        :return: url_id, description, url, mock_id, file_ext
        """
        if url_id is None:
            url_id = None
            description = None
            url = None
            file_ext = None
            method = None
            mock_id = self.urls["mock_id"]
        else:
            assert url_id in self.urls["endpoints"].keys(), "Url id not implemented, please check urls.endpoints"
            url_id = url_id
            description = self.urls["endpoints"][url_id]["description"]
            url = self.urls["endpoints"][url_id]["url"]
            mock_id = self.urls["mock_id"]
            file_ext = self.urls["endpoints"][url_id]["file_ext"]
            method = self.urls["endpoints"][url_id]["method"]
        return url_id, description, url, mock_id, file_ext, method


class ConsurfDBUrl(Url):
    def __init__(self, url_id="consurf_1"):
        """
        Class used to call the URL: http://bental.tau.ac.il/new_ConSurfDB/DB/
        """
        Url.__init__(self)
        self.urls = {"base_url":"http://bental.tau.ac.il/new_ConSurfDB/DB/",
                     "base_description":"ConSurfDB, proteins conservations",
                     "mock_id":"prot_id/chain_id",
                     "endpoints":{
                         "consurf_1":{
                             "url": "http://bental.tau.ac.il/new_ConSurfDB/DB/prot_id/chain_id/consurf.grades",
                             "description":"-",
                             "file_ext":".cs",
                             "method":"get"
                         }
                     }
                }

        self.base_url = self.urls["base_url"]
        self.base_description = self.urls["base_description"]
        self.url_id, self.description, self.url, self.mock_id, self.file_ext, self.method = self.config_url(url_id)


class PdbRedoUrl(Url):
    def __init__(self, url_id = None):
        """
        Class used to call the URL: https://pdb-redo.eu/db
        :param url_id: id of endpoint to be used from URL
        """
        Url.__init__(self, url_id=url_id)
        self.urls = {"base_url":"https://pdb-redo.eu/db/",
                     "base_description":"PDB-REDO database",
                     "mock_id":"9xyz",
                     "endpoints":{
                         "pdbredo_1":{
                             "url":"https://pdb-redo.eu/db/9xyz/9xyz_0cyc.pdb.gz",
                             "description":"Atomic coordinates:Initial model",
                             "file_ext":".pdb.gz",
                             "method": "get"},
                         "pdbredo_2":{
                             "url":"https://pdb-redo.eu/db/9xyz/9xyz_besttls.pdb.gz",
                             "description":"Atomic coordinates:Re-refined (only) model",
                             "file_ext": ".pdb.gz",
                             "method": "get"},
                         "pdbredo_3":{
                             "url":"https://pdb-redo.eu/db/9xyz/9xyz_final.pdb",
                             "description":"Atomic coordinates:Re-refined & rebuilt structure model",
                             "file_ext": ".pdb",
                             "method": "get"},
                         "pdbredo_4":{
                             "url":"https://pdb-redo.eu/db/9xyz/9xyz_final_tot.pdb",
                             "description":"Atomic coordinates:Re-refined & rebuilt structure model with total B-factors (PDB)",
                             "file_ext": ".pdb",
                             "method": "get"},
                         "pdbredo_5":{
                             "url":"https://pdb-redo.eu/db/9xyz/9xyz_final.cif",
                             "description":"Atomic coordinates:Re-refined & rebuilt structure model with total B-factors (mmCIF)",
                             "file_ext": ".cif",
                             "method": "get"}
                     }
                }

        self.base_url = self.urls["base_url"]
        self.base_description = self.urls["base_description"]
        self.url_id, self.description, self.url, self.mock_id, self.file_ext, self.method = self.config_url(url_id)


class PdbUrl(Url):

    def __init__(self, url_id="pdb_1"):
        """
        Class used to call the URL: https://files.rcsb.org/download
        :param url_id: id of endpoint to be used from URL
        """
        Url.__init__(self, url_id=url_id)
        self.urls = {"base_url":"https://files.rcsb.org/download/",
                     "base_description":"rscb.org - PDB",
                     "mock_id":"123abc",
                     "endpoints":{
                         "pdb_1": {"url":"https://files.rcsb.org/download/123abc.pdb",
                                   "description":"Biological Assembly File in PDB	",
                                   "file_ext":".pdb",
                                   "method": "get"}
                     }
                }
        self.base_url = self.urls["base_url"]
        self.base_description = self.urls["base_description"]
        self.url_id, self.description, self.url, self.mock_id, self.file_ext, self.method = self.config_url(url_id)


class SiftsUrl(Url):
    def __init__(self, url_id = "sifts_1"):
        """
        Class used to call the URL: http://www.ebi.ac.uk/pdbe/docs/sifts/
        :param url_id: id of endpoint to be used from URL
        """
        Url.__init__(self, url_id=url_id)
        self.urls = {"base_url":"http://www.ebi.ac.uk/pdbe/api/mappings/",
                     "base_description":"""
                                            Structure Integration with Function, Taxonomy and Sequence (SIFTS) is a 
                                            project in the PDBe-KB resource for residue-level mapping between UniProt 
                                            and PDB entries. SIFTS also provides annotation from the IntEnz, GO, 
                                            InterPro, Pfam, CATH, SCOP, PubMed, Ensembl and Homologene resources. 
                                            The information is updated and released every week concurrently with the 
                                            release of new PDB entries and is widely used by resources such as RCSB PDB,
                                             PDBj, PDBsum, Pfam, SCOP and InterPro. 
                                             """,
                     "mock_id":"123abc",
                     "endpoints": {
                         "sifts_1": {"url": "http://www.ebi.ac.uk/pdbe/api/mappings/123abc",
                                     "description": """
                                                    Mappings (as assigned by the SIFTS process) from PDB structures to 
                                                    UniProt, Pfam, InterPro, CATH, SCOP, IntEnz, GO, Ensembl and HMMER 
                                                    accessions (and vice versa). """,
                                     "file_ext": ".json",
                                     "method":"get"},
                         "sifts_2": {"url": "http://www.ebi.ac.uk/pdbe/api/mappings/uniprot/123abc",
                                     "description": """
                                                        Mappings (as assigned by the SIFTS process) from PDB structures 
                                                        to UniProt.  """,
                                     "file_ext": ".json",
                                     "method":"get"},
                         "sifts_3": {"url": "http://www.ebi.ac.uk/pdbe/api/mappings/pfam/123abc",
                                     "description": """
                                                        Mappings (as assigned by the SIFTS process) from PDB structures 
                                                        to Pfam.   """,
                                     "file_ext": ".json",
                                     "method": "get"}
                     }
            }
        self.base_url = self.urls["base_url"]
        self.base_description = self.urls["base_description"]
        self.url_id, self.description, self.url, self.mock_id, self.file_ext, self.method = self.config_url(url_id)


class CapriUrl(Url):
    def __init__(self, url_id="capri_1"):
        """
        Class used to call the URL: https://www.ebi.ac.uk/msd-srv/capri/pdb_ids.html
        :param url_id: id of endpoint to be used from URL
        """
        Url.__init__(self, url_id=url_id)
        self.urls = {"base_url": "http://www.ebi.ac.uk/pdbe/api/mappings/",
                     "base_description": """CAPRI communitywide experiment on the comparative evaluation of 
                     protein-protein docking for structure prediction""",
                     "mock_id": "123abc",
                     "endpoints": {
                         "capri_1": {"url": "https://www.ebi.ac.uk/msd-srv/capri/pdb_ids.html",
                                     "description": """ CAPRI target publications""",
                                     "file_ext": ".html",
                                     "method": "get"},
                         "capri_2": {"url": "https://www.ebi.ac.uk/msd-srv/capri/orig/123abc.brk",
                                     "description": """ CAPRI starting file""",
                                     "file_ext": ".brk",
                                     "method": "get"}
                     }
                     }
        self.base_url = self.urls["base_url"]
        self.base_description = self.urls["base_description"]
        self.url_id, self.description, self.url, self.mock_id, self.file_ext, self.method = self.config_url(url_id)
