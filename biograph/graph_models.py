import networkx as nx
import pandas as pd
from itertools import combinations
from getcontacts.contact_calc.compute_contacts import compute_contacts
from getcontacts.contact_calc import transformations
from getcontacts.Applications.contact_network_analysis import create_graph
import tempfile
import sys

from biograph.constants import valid_amino_3

class GraphModel:
    def __init__(self):
        self.G = nx.Graph()
    def generate_graph(self, protein, params):
        raise NotImplementedError
    def add_features(self, dataframe, columns):
        raise NotImplementedError
    @staticmethod
    def graph_to_dataframe(G):
        return pd.DataFrame.from_dict({node_idx: {feature: val for feature,val in G.nodes[node_idx].items()}
              for node_idx in G.nodes}, orient="index")
    @staticmethod
    def get_diffused_graph(G, aggregator = None, keys = None, steps = 1):
        """Diffuses some `keys` features across graph neighbors that are
        `steps` apart and afterwards process the groups using `aggregator`.
        Diffused features end with "_n" with n being the steps to the node.

        Parameters
        ----------
        aggregator: function or None
            A function that receives a dictionary of list of diffused
            features (e.g. {mass_1 => [1,12], mass_2 => [8,1,13]}) and
            returns a new dictionary.
            If None, non-numeric (i.e. not complex, int, float or bool)
            features are discarded and numeric features are averaged.
        keys: list
            Features to be diffused. All nodes in graph must have them.
        steps: int
            Number of steps from each node to be used in diffusion.

        Returns
        -------
        networkx.Graph with same structure as before but with additional
        diffused features.
        """
        if keys is None:
            some_node = list(G.nodes)[0]
            keys = G.nodes[some_node].keys()

        diffused_graph = G.copy()
        # Setup feature holders.
        for node_idx in diffused_graph.nodes:
            for key in keys:
                for dist in range(1, steps+1):
                    diffused_graph.nodes[node_idx]["{}_{}".format(key, dist)] = list()

        for node_idx in G.nodes:
            # nx doesn't offer a bfs traversal that also yields distances,
            # so we need to keep track of them.
            distances = {node_idx: 0}
            new_keys = set()
            for origin, neighbor in nx.bfs_edges(G, node_idx):
                distances[neighbor] = distances[origin] + 1
                if distances[neighbor] > steps:
                    break
                for key in keys:
                    neighbor_value = diffused_graph.nodes[neighbor][key]
                    diffused_key = "{}_{}".format(key, distances[neighbor])
                    new_keys.add(diffused_key)
                    try:
                        diffused_graph.nodes[node_idx][diffused_key].append(neighbor_value)
                    except KeyError:
                        print(neighbor_value)
                        print(keys)
                        print(diffused_key)
                        print(diffused_graph.nodes[node_idx])
                        return

        # Process feature lists with an aggregator.
        for node_idx in diffused_graph.nodes:
            node_features = diffused_graph.nodes[node_idx].copy()

            diffused_features = {key:val
                for key, val in node_features.items()
                if key in new_keys}

            if aggregator is None:
                # Basic aggregator. If the type is numeric-like then do a mean, otherwise
                # discard it.
                diffused_features = {k:v for k,v in diffused_features.items()
                    if isinstance(v[0], (int, float, complex))}
                for diffused_key in diffused_features.keys():
                    feature_list = diffused_features[diffused_key]
                    diffused_features[diffused_key] = sum(feature_list)/len(feature_list)

            else:
                diffused_features = aggregator(diffused_features)

            # We can't replace the dict object, so here:
            diffused_graph.nodes[node_idx].clear()
            diffused_graph.nodes[node_idx].update(diffused_features)

            # Add back the original node-specific values ("distance zero" and ignored
            # columns).
            for key in node_features.keys():
                if key not in new_keys:
                    diffused_graph.nodes[node_idx][key] = node_features[key]

        return diffused_graph


class StructureGraphGenerator(GraphModel):
    def generate_graph(self, obj, params):
        """
        Generates a graph from the structure model of the protein.
        Protein.generate_graph double dispatches to this method.
        Parameters
        ----------
        obj: Protein
        params: dict
            - step: step of the structure generator used to create the graph
            - surface_only (optional, default False): whether to discard interior
              points and edges or not.
        """
        if obj.structure is None:
            raise Exception("To execute StructureGraphGenerator first it is"
                            "needed for the Protein to have a structure model")

        step = params["step"]
        core_edges = set()
        for t in obj.structure._get_core(step):
            for e in combinations(t, 2):
                core_edges.add(e)
        in_surf = {e: False for e in core_edges}
        for t in obj.structure.get_simplices_by_step(step):
            for e in combinations(t, 2):
                in_surf[e] = True
        if params.get("surface_only", False):
            self.G.add_edges_from([(
                e[0],
                e[1],
                {'weight': obj.structure.get_edge_length(e),
                'in_surf': in_surf[e]})
                for e in core_edges if in_surf[e]])
        else:
            self.G.add_edges_from([(
                e[0],
                e[1],
                {'weight': obj.structure.get_edge_length(e),
                'in_surf': in_surf[e]})
                for e in core_edges])

        # We'll need to keep track of atom ids.
        for node_id in self.G.nodes:
            self.G.nodes[node_id]["full_id"] = obj.structure.point_ids[node_id]

        if params.get("surface_only", False):
            # Remove interior points
            nodes_to_remove = []
            for idx, adj_dict in self.G.adjacency():
                if len(adj_dict) == 0:
                    nodes_to_remove.append(idx)
            for idx in nodes_to_remove:
                self.G.remove_node(idx)

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

class StaticContactGraphGenerator(GraphModel):
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

    def add_features(self, dataframe, agg = dict(), columns = [
            "bfactor", "score", "color",
            "color_confidence_interval_high", "color_confidence_interval_low",
            "score_confidence_interval_high", "score_confidence_interval_low"]):
        """
        Add features to the static contact graph based on a dataframe.
        Dataframe must have res_full_id with standard Bio.PDB format
        and resname with three-letter aminoacid code.
        Since each node in the graph represents a residue, rows that match
        the residue full id are grouped and mean() aggregated to insert the
        values in the graph node. For object types, selects the first element.
        You can change this behavior by using the `agg` parameter.

        Parameters
        ----------
        dataframe : pd.DataFrame
            Dataframe to get features from.
        agg : dict
            Dictionary mapping column names to functions or function names
            (see pandas.DataFrame.aggregate for more info.)
            By default 'mean' for numeric types and first-element for
            object types.
        columns : list
            Columns to be included. Must be numeric. If many values
            are repeated for the same residue (e.g. if the dataframe
            is at the atom level), they are averaged.

        Returns
        -------
        nx.Graph
        """
        # By default behavior is 'mean' for numeric datatypes and 'first element'
        # for all others (objects etc.)
        aggfn = {
            col: "mean" if pd.api.types.is_numeric_dtype(dtype) else lambda L: L[0]
            for col, dtype in dataframe.dtypes.to_dict().items()
            if col in columns
        }
        aggfn.update(agg)
        features = dataframe.groupby(["res_full_id", "resname"]).agg(aggfn).reset_index()
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
