import networkx as nx
from itertools import combinations


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
