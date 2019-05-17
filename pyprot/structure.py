import numpy as np
from scipy.spatial import Delaunay
from itertools import combinations, chain
import os
from pyprot.constants import PERSEUSPATH
import subprocess
import networkx as nx
from pyprot.io import Writer
import pyprot


class StructureModel:

    def __init__(self, points, graph_id="structure"):
        self.__points = points
        self.graph_id = graph_id
        self.simplices_order = None
        self.persistent_hom_params = {}
        self.delaunay = None
        self._initiate_delaunay()
        self._set_simplices_order()

    @property
    def points(self):
        return self.__points

    @points.setter
    def points(self, points):
        assert isinstance(points, np.ndarray), """A Numpy array matrix must be
                                                        used for points"""
        self.__points = points

    @points.deleter
    def points(self):
        del self.__points

    def _initiate_delaunay(self):
        """Method used to instantiate scipy.spatial.Delaunay using self.points.
        In addition, simplices are sorted. For internal use only.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        self.delaunay = Delaunay(points=self.points)
        self.delaunay.simplices.sort()

    def _set_simplices_order(self, order="auto", mode="size"):
        """Generate numpy array with simplices ordered by a criteria.
        Default mode is order simplices by size

        Parameters
        ----------
        order : str or (N,) numpy ndarray of ints
            How to order simplices. If str, it must be 'auto' and the Function
            will use mode to decide how to order. If numpy ndarray, the orden
            of the simplices
        mode : str
            mode: criteria used to order simplices. Default is to order
            simplices by size

        Returns
        -------
        None
        """
        simplices_order = None
        if order == "auto":
            if mode == "size":
                combs = chain.from_iterable(combinations(list(range(4)), 2))
                pair_indices = np.fromiter(combs, int, count=12)
                pair_indices = pair_indices.reshape(6, 2).T
                points_diff = self.points[self.delaunay.simplices][:, pair_indices[0]] - \
                                                 self.points[self.delaunay.simplices][:, pair_indices[1]]
                simplices_sizes = np.sort(np.linalg.norm(points_diff, axis=-1),
                                          axis=1)
                # Lexicographic comparison
                simplices_order = np.lexsort(simplices_sizes.T)
        else:
            simplices_order = order

        self.simplices_order = simplices_order

    def diffuse(self, diffuser_matrix, steps=1):
        """Diffuse information using a matrix diffuser.

        Parameters
        ----------
        diffuser_matrix : numpy Array
            matrix used to diffuse information
        steps : int
            number of times the information is diffused using diffuser_matrix

        Returns
        -------
        numpy Array
            Matrix diffused
        """
        assert steps >= 1, "steps must be equal or bigger than 1"

        for _ in range(steps):
            diff_matrix = np.array(np.matmul(diffuser_matrix, self.points))
        return diff_matrix

    def count_faces(self, simplices):
        """Given a set of simplices, return the unique faces, and the number of
         ocurrences of each one

        Parameters
        ----------
        simplices: numpy array
             array of simplices

        Returns
        -------
        numpy array
            unique faces
        numpy array
            number of time each face appears

        """
        faces = self.get_faces(simplices)
        unique_faces, counted_faces = np.unique(np.concatenate(faces, axis=0),
                                                axis=0, return_counts=True)
        return unique_faces, counted_faces

    def get_simplices_by_step(self, step):
        """Given a step returns the faces that make the surface at that step

        Parameters
        ----------
        step : int
            Step of Perseus's simtop method

        Returns
        -------
        numpy array
            Faces that make the surface
        """
        simplices = self._get_core(step)
        unique_faces, counted_faces = self.count_faces(simplices)
        return unique_faces[counted_faces == 1]

    def get_surface(self, step):
        """ DEPRECATED
        Given a fat step returns the faces that make the surface
        :param fat_step: int
        :return: numpy array, faces that make the surface
        """
        simplices = self.delaunay.simplices[self.simplices_order[:step]]
        unique_faces, counted_faces = self.count_faces(simplices)
        return unique_faces[counted_faces == 1]

    def get_edge_length(self, e):
        """ Given an edge, return its length

        Parameters
        ----------
        e : tuple or list
            tuple or list with two elements that represents two points in the
            structure. These points define an edge.

        Returns
        -------
        float
            length of edge
        """
        # TODO : check return DocBlock
        return np.linalg.norm(self.delaunay.points[e[0]]-self.delaunay.points[e[1]])

    def _get_min_steps(self, order):
        """Method used to get minimum number of steps so as simplices cover all
         the points in the structure

        Parameters
        ----------
        order : numpy array
            indices used to order delaunay simplices

        Returns
        -------
        int
            Minimum number of steps to cover all points.
        """
        min_n_steps = None
        for step in range(self.delaunay.nsimplex):
            simplices_ = self.delaunay.simplices[order[:step]]
            if np.unique(simplices_).shape[0] == self.delaunay.npoints:
                min_n_steps = step
                break
        return min_n_steps

    def export_graph(self, as_networkx=True, as_adjancency=True,
                     surface_graph=True, step="default"):
        # TODO: comentar y mover a GraphModel
        """
        :param as_networkx:
        :param as_adjancency:
        :param surface_graph:
        :param step:
        :return: dict
        """
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

    def calculate_depth(self, surface_graph):
        """Calculates depth on the graph surface.

        Parameters
        ----------
        surface_graph : nx.Graph
            Graph of surface.

        Returns
        -------
         tuple of two dictionaries keyed by node
            The first dictionary stores distance from one of the source nodes.
            The second stores the path from one of the sources to that node.
        """
        level = list(self.get_vertices(self.get_simplices_by_step(step=-1)))
        return nx.multi_source_dijkstra(surface_graph, level)

    @staticmethod
    def _faces(s):
        return [tuple(np.delete(s, i)) for i in range(s.shape[0])]

    def _get_core(self, step):
        """Get simplices up to a given step

        Parameters
        ----------
        step : int
            step in the Perseus' simtop algorithm

        Returns
        -------
        np.Array
            Array of simplices.

        """
        # TODO: check return DockBlock
        return self.delaunay.simplices[self.simplices_order[:step]]

    def _step_to_length(self, i):
        # TODO: check
        return max(self.edge_size(self.delaunay.simplices[self.simplices_order[i]]))

    @staticmethod
    def read_betti(betti_filename, step=None):
        # TODO: move to other module
        """
        Reads a Betti file
        :param betti_filename: str, betti file path
        :param fat_step:
        :param step:
        :return:
        """
        ifile = open(betti_filename)
        betti = [list(map(int, x.split())) for x in ifile.readlines() if
                 x.split()]
        ifile.close()
        betti = [x for x in betti if x[0] <= step]
        return betti

    @staticmethod
    def get_faces(simplices):
        """From a K dim array of n-simplices gets a cube of facesself. That is,
        given an array of shape (K,n+1) returns an array of shape (K, n+1, n)
        of the n+1 faces of each of the K simplices.

        Parameters
        ----------
        simplices: numpy array
             array of simplices

        Returns
        -------
        numpy array
            Cube of faces of the n-simplices
        """
        dim = simplices.shape[1]
        return np.stack([np.delete(simplices, i, axis=1) for i in range(dim)],
                        axis=1)

    @staticmethod
    def get_vertices(simplices):
        """Returns vertices involved in given array of n-simplices.

        Parameters
        ----------
        simplices: numpy array
             array of simplices

        Returns
        -------
        numpy array
            Unique vertices in given simplices
        """
        return np.unique(simplices)


class Perseus:
    def __init__(self):
        MODULEDIR = os.path.dirname(os.path.abspath(__file__))
        self.PERSEUSPATH = os.path.join(MODULEDIR, "perseus", "perseus")

    def calculate_fatcore(self, structure, topology=(1, 0, 0),
                          perseus_input_filename="ordered_simplices.pers",
                          output_dir="perseus", output_filename="pdb"):
        """
        Function used to calculate fatcore. This function uses Perseus's
        simtop method.
        Perseus's documentation can be found in:
        http://people.maths.ox.ac.uk/nanda/perseus/index.html

        Parameters
        ----------
        topology : tuple or list
            topology is used to define different steps of interest.
        output_dir : str
            directory where output will be written
        output_filename : str
            name for output file

        Returns
        -------
        StructureModel
            StructureModel with persistent_hom_params computed
        """
        # TODO: ver con Lean
        structure_ = structure
        if isinstance(structure_, pyprot.protein.Protein):
            structure = structure.structure

        # TODO: validar que structure tenga datos
        # TODO: validar que sea Protein o StructeModel

        dir_path = os.path.join(output_dir, output_filename)
        Writer(dir_path).create_directory()

        with open(os.path.join(dir_path, perseus_input_filename),
                  'w') as ofile:
            # These are the dimensions of the array passed to perseus
            ofile.write('3\n1\n')
            input = np.append(
                structure.delaunay.simplices[structure.simplices_order],
                np.arange(1, 1 + structure.delaunay.nsimplex).reshape(-1, 1),
                axis=1)
            np.savetxt(ofile, input, delimiter=' ', newline='\n', fmt='%i')
        min_step = structure._get_min_steps(order=structure.simplices_order)
        to_run = (self.PERSEUSPATH, 'simtop',
                  os.path.join(dir_path, perseus_input_filename),
                  os.path.join(dir_path, "perseus_output"))
        popen = subprocess.Popen(to_run, stdout=subprocess.PIPE)
        popen.wait()
        with open(os.path.join(output_dir, output_filename,
                               'perseus_output_betti.txt'), 'r') as ifile:
            small_step = min_step
            fat_step = min_step  # -1
            big_step = min_step  # -1
            for line in ifile.readlines()[1:]:
                # first line is blank
                i, h0, h1, h2, _ = [int(x) for x in line.split()]
                if i >= min_step and h0 == topology[0]:

                    # only modifies the fat_step once
                    if small_step == min_step:
                        small_step = i

                    # only modifies the fat_step once
                    if fat_step == min_step and h1 == topology[1]:
                        fat_step = i
                    if h1 == topology[1] and h2 == topology[2]:
                        big_step = i
                        break
                # TODO: check why is this assignment here
                lh0 = h0

        if isinstance(structure_, pyprot.protein.Protein):
            structure_.structure.persistent_hom_params["min_step"] = min_step
            structure_.structure.persistent_hom_params["small_step"] = small_step
            structure_.structure.persistent_hom_params["fat_step"] = fat_step
            structure_.structure.persistent_hom_params["big_step"] = big_step
        else:
            structure_.persistent_hom_params["min_step"] = min_step
            structure_.persistent_hom_params["small_step"] = small_step
            structure_.persistent_hom_params["fat_step"] = fat_step
            structure_.persistent_hom_params["big_step"] = big_step
        return structure_
