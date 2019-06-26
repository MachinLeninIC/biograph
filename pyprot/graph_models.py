import networkx as nx
from itertools import combinations
from getcontacts.contact_calc.compute_contacts import compute_contacts
from getcontacts.contact_calc import transformations
from getcontacts.Applications.contact_network_analysis import create_graph
import tempfile

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

class StaticConctactGraphGenerator:
    def __init__(self):
        self.G = nx.Graph()

    def _parse_params(self, params):
        params.setdefault("itypes", None)
        params.setdefault("geom_criteria", )
        params.setdefault("cores", 2)
        params.setdefault("beg", )
        params.setdefault("end", )
        params.setdefault("stride", )
        params.setdefault("distout", )
        params.setdefault("ligand", )
        params.setdefault("solv", )
        params.setdefault("lipid", )
        params.setdefault("sele1", )
        params.setdefault("sele2",  )
        params.setdefault("save_contact_filename", None)
        return params

    def generate_graph(self, protein, params):
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

        # If instead of the trajectory you pass it the topology again, then
        # it calculates static contacts.
        compute_contacts(pdb_file, pdb_file, contacts_filename,
            params["itypes"], params["geom_criteria"], params["cores"],
            params["beg"], params["end"], params["stride"], params["distout"],
            params["ligand"], params["solv"], params["lipid"], params["sele1"],
            params["sele2"])
        with open(contacts_filename, "r") as handle:
            atom_interactions = transformations.parse_contacts(handle, params["itypes"])
        residue_interactions = transformations.res_contacts(atom_interactions)
        interaction_counts = transformations.gen_counts(residue_interactions)

        for (res1, res2), num_interactions in interaction_counts.items():
            self.G.add_node(res1)
            self.G.add_node(res2)
            self.G.add_edge(res1, res2, weight=num_interactions)

        return self.G