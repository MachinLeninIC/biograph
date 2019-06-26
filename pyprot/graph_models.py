import networkx as nx
from itertools import combinations
from getcontacts.contact_calc.compute_contacts import compute_contacts
from getcontacts.contact_calc import transformations
from getcontacts.Applications.contact_network_analysis import create_graph
import tempfile
import sys

class StructureGraphGenerator:
    def __init__(self):
        self.G = nx.Graph()

    def generate_graph(self, obj, params):
        if obj.structure is None:
            raise Warning("""To execute StructureGraphGenerator first is
                                needed that Protein has a structure model""")
        else:
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
        params.setdefault("itypes", "all") # it is not actually optional
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
        sys.stdout = tempfile.TemporaryFile()
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
        with open(contacts_filename, "r") as handle:
            # We want all interaction types, so we pass None.
            atom_interactions = transformations.parse_contacts(handle, None)
        residue_interactions = transformations.res_contacts(atom_interactions)
        interaction_counts = transformations.gen_counts(residue_interactions)

        for (res1, res2), num_interactions in interaction_counts.items():
            self.G.add_node(res1)
            self.G.add_node(res2)
            self.G.add_edge(res1, res2, weight=num_interactions)

        return self.G