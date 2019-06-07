import networkx as nx
from itertools import combinations


class StructureGraphGenerator:
    def __init__(self):
        self.G = nx.Graph()

    def generate_graph(self, structure, core_step, surface_step):
        core_edges = set()
        for t in structure._get_core(core_step):
            for e in combinations(t, 2):
                core_edges.add(e)
        in_surf = {e: False for e in core_edges}
        for t in structure.get_simplices_by_step(surface_step):
            for e in combinations(t, 2):
                in_surf[e] = True
        self.G.add_edges_from(
            [(e[0], e[1],
              {'weight': structure.get_edge_length(e),
               'in_surf': in_surf[e]}) for e in core_edges])

    def export_graphs(self, as_networkx, surface_graph):
        if as_networkx:
            results["full_networkx_graph"] = G
            if as_adjancency:
                results["full_adjacency"] = [i for i in G.adjacency()]

        if surface_graph:
            G_surf = nx.Graph([e for e in G.edges(data=True) if e[-1]['in_surf']])
            results["surface_networkx_graph"] = G_surf
            if as_adjancency:
                results["surface_adjacency"] = [i for i in G_surf.adjacency()]


## contactos (se)
## 1- Tenes dos residuos y hay que ver si los atomos estan a una cierta distancia
## 2- Ver si hay un enlace quimico entre los puntos (probablemente sea lo mismo que 1)
## 3-



    def export_graph(self, as_networkx=True, as_adjancency=True,
                     surface_graph=True, step="default"):
        # TODO: comentar y mover a GraphModel
        results = {}
        if step == "default":
            step = self.persistent_hom_params["fat_step"]
        else:
            assert type(step) is int, "step must be 'default' or int"

        core_edges = set()
        for t in self._get_core(step):
            for e in combinations(t, 2):
                core_edges.add(e)
        in_surf = {e: False for e in core_edges}
        for t in self.get_simplices_by_step(step):
            for e in combinations(t, 2):
                in_surf[e] = True
        G = nx.Graph()
        G.add_edges_from([(e[0], e[1], {'weight': self.get_edge_length(e), 'in_surf':in_surf[e]}) for e in core_edges])

        if as_networkx:
            results["full_networkx_graph"] = G
            if as_adjancency:
                results["full_adjacency"] = [i for i in G.adjacency()]

        if surface_graph:
            G_surf = nx.Graph([e for e in G.edges(data=True) if e[-1]['in_surf']])
            results["surface_networkx_graph"] = G_surf
            if as_adjancency:
                results["surface_adjacency"] = [i for i in G_surf.adjacency()]

        return results
