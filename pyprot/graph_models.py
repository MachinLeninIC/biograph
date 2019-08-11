import networkx as nx
from itertools import combinations
from getcontacts.contact_calc.compute_contacts import compute_contacts
from getcontacts.contact_calc import transformations
from getcontacts.Applications.contact_network_analysis import create_graph
import tempfile
import sys

from pyprot.constants import valid_amino_3

class StructureGraphGenerator:
    def __init__(self):
        self.G = nx.Graph()

    def generate_graph(self, obj, params):
        if obj.structure is None:
            raise Exception("""To execute StructureGraphGenerator first is
                                needed that Protein has a structure model""")

        step = params["step"]
        core_edges = set()
        for t in obj.structure._get_core(step):
            for e in combinations(t, 2):
                core_edges.add(e)
        in_surf = {e: False for e in core_edges}
        for t in obj.structure.get_simplices_by_step(step):
            for e in combinations(t, 2):
                in_surf[e] = True
        self.G.add_edges_from(
            [(e[0], e[1],
                {'weight': obj.structure.get_edge_length(e),
                'in_surf': in_surf[e]}) for e in core_edges])

        # We'll need to keep track of atom ids.
        for node_id in self.G.nodes:
            self.G.nodes[node_id]["full_id"] = obj.structure.point_ids[node_id]

        return self.G

    def add_features(self, dataframe, columns = ["bfactor", "score", "color",
            "color_confidence_interval_high", "color_confidence_interval_low",
            "score_confidence_interval_high", "score_confidence_interval_low"]):
        """
        Add features to the structure graph based on a dataframe.
        Dataframe must have full_id with standard Bio.PDB format. Features
        are added based on the atom ids that were selected when the
        Protein.structure was created.

        Parameters
        ----------
        dataframe : pd.DataFrame
            Dataframe to get features from.
        columns : list
            Columns to be included.

        Returns
        -------
        nx.Graph
        """
        columns = [col for col in columns if col in dataframe.columns]

        # In this graph, we identify each node by the ID of the point that
        # was passed to the Delaunay algorithm.
        for node_id in self.G.nodes:
            full_id = self.G.nodes[node_id]["full_id"]

            # Get a single atom row as a dict.
            features = dataframe[dataframe.full_id == full_id].reset_index(drop=True)
            features = {key:val[0] for key,val in features.to_dict().items()}

            for col in columns:
                self.G.nodes[node_id][col] = features[col]

        return self.G

    '''def export_graphs(self, as_networkx, surface_graph):
        if as_networkx:
            results["full_networkx_graph"] = G
            if as_adjancency:
                results["full_adjacency"] = [i for i in G.adjacency()]

        if surface_graph:
            G_surf = nx.Graph([e for e in G.edges(data=True) if e[-1]['in_surf']])
            results["surface_networkx_graph"] = G_surf
            if as_adjancency:
                results["surface_adjacency"] = [i for i in G_surf.adjacency()]'''

class StaticContactGraphGenerator:
    def __init__(self):
        self.G = nx.Graph()
    def _parse_geometric_criteria(self, geom):
        """
        Set sensible defaults for the geometric parameters of the model.
        Values are taken from the argparsers in the getcontacts library.
        """
        # Distances are in angstroms, angles are in degrees.
        geom.setdefault("SALT_BRIDGE_CUTOFF_DISTANCE", 4.0)
        geom.setdefault("PI_CATION_CUTOFF_DISTANCE", 6.0)
        geom.setdefault("PI_CATION_CUTOFF_ANGLE", 60)
        geom.setdefault("PI_STACK_CUTOFF_DISTANCE", 7.0)
        geom.setdefault("PI_STACK_CUTOFF_ANGLE", 30)
        geom.setdefault("PI_STACK_PSI_ANGLE", 45)
        geom.setdefault("T_STACK_CUTOFF_DISTANCE", 5.0)
        geom.setdefault("T_STACK_CUTOFF_ANGLE", 30)
        geom.setdefault("T_STACK_PSI_ANGLE", 45)
        geom.setdefault("HBOND_CUTOFF_DISTANCE", 3.5)
        geom.setdefault("HBOND_CUTOFF_ANGLE", 180)
        geom.setdefault("HBOND_RES_DIFF", 1)
        geom.setdefault("VDW_EPSILON", 0.5)
        geom.setdefault("VDW_RES_DIFF", 2)
        return geom

    def _parse_params(self, params):
        """
        Set sensible defaults for the parameters of the model. Values are taken
        from the get static contacts scripts in the getcontacts library.
        """
        params.setdefault("itypes", ["sb", "pc", "ps", "ts", "vdw", "hb", "hp"])
        params.setdefault("geom_criteria", dict())
        params["geom_criteria"] = self._parse_geometric_criteria(params["geom_criteria"])
        params.setdefault("distout", False)
        params.setdefault("ligand", "")
        params.setdefault("solv", None)
        params.setdefault("lipid", "")
        params.setdefault("sele1", "protein")
        params.setdefault("sele2", params["sele1"])
        params.setdefault("save_contact_filename", None)
        return params

    def generate_graph(self, protein, params):
        """
        Generate a contact static graph from the PDB file.
        Parameters
        ----------
        protein: the Protein object. Either it was created with a file name as
        argument or you prove a 'pdb_file' in params.
        params: parameters for the model, including geometry criteria. It has
        sensible defaults. You can view the full list of parameters in the
        getcontact package.

        Returns
        -------
        graph: a networkx.Graph instance containing the full static contact graph
        for the given protein.
        """
        params = self._parse_params(params)
        pdb_file = protein.pdb_file
        pdb_file = pdb_file if "pdb_file" not in params else params["pdb_file"]
        if pdb_file is None:
            raise Exception("You need to provide a pdb_file as a parameter.")

        # Unfortunately, the library requires a filename to operate.
        # And it actually calls open().
        if params["save_contact_filename"] is not None:
            contacts_filename = params["save_contact_filename"]
        else:
            tempdir = tempfile.TemporaryDirectory()
            contacts_filename = "{}/{}_contacts.tsv".format(tempdir.name,
                pdb_file.replace(".pdb", ""))

        # The library also manipulates the stdout file descriptors, which
        # triggers unsupported behavior in iPython-like environments (e.g.
        # Jupyter). Therefore we need to mask it for a bit.
        _stdout = sys.stdout
        sys.stdout = tempfile.TemporaryFile(mode="w+")
        # If instead of the trajectory you pass it the topology again, then
        # it calculates static contacts.
        # Since we are not using trajectories, cores is always 1, beg and end
        # are 0, and stride is 1.
        compute_contacts(pdb_file, pdb_file, contacts_filename,
            params["itypes"], params["geom_criteria"],
            1, 0, 0, 1,
            params["distout"], params["ligand"], params["solv"], params["lipid"],
            params["sele1"], params["sele2"])

        # Give back the stdout.
        sys.stdout= _stdout
        lines = []
        with open(contacts_filename, "r") as handle:
            line = handle.readline()
            while line:
                lines.append(line)
                line = handle.readline()
        # We want all interaction types, so we pass None.
        atom_interactions, num_frames = transformations.parse_contacts(lines, None)
        residue_interactions = transformations.res_contacts(atom_interactions)
        interaction_counts = transformations.gen_counts(residue_interactions)

        for (res1, res2), num_interactions in interaction_counts.items():
            self.G.add_node(res1)
            self.G.add_node(res2)
            self.G.add_edge(res1, res2, weight=num_interactions)

        return self.G

    def add_features(self, dataframe, columns = ["bfactor", "score", "color",
            "color_confidence_interval_high", "color_confidence_interval_low",
            "score_confidence_interval_high", "score_confidence_interval_low"]):
        """
        Add features to the static contact graph based on a dataframe.
        Dataframe must have res_full_id with standard Bio.PDB format
        and resname with three-letter aminoacid code.

        Parameters
        ----------
        dataframe : pd.DataFrame
            Dataframe to get features from.
        columns : list
            Columns to be included. Must be numeric. If many values
            are repeated for the same residue (e.g. if the dataframe
            is at the atom level), they are averaged.

        Returns
        -------
        nx.Graph
        """
        features = dataframe.groupby(["res_full_id", "resname"]).mean().reset_index()
        columns = [col for col in columns if col in features.columns]

        count_missing = 0

        for index, row in features.iterrows():
            # Sometimes there are water molecules or things that are not aminoacids
            # in the dataframe.
            if not valid_amino_3(row["resname"]):
                continue

            #res_full_id has format (pdb_name, 0, chain, (atom, residue_inumber, t))
            #residue identifier in get_contact is chain:amino_3code:residue_inumber
            res_id = "{}:{}:{}".format(row["res_full_id"][2], row["resname"], row["res_full_id"][3][1])

            if res_id not in self.G.nodes:
                print("Missing node here: {} - {} aka {}".format(
                    row["res_full_id"], row["resname"], res_id))
                count_missing+=1
                continue

            for col in columns:
                self.G.nodes[res_id][col] = row[col]

        print("Missing nodes: {}".format(count_missing))
        return self.G
