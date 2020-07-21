import numpy as np
from scipy.spatial import Delaunay
from itertools import combinations, chain
import os
import subprocess
import networkx as nx
from biograph.io import Writer
import biograph
import matplotlib.pyplot as plt
import pandas as pd
import gudhi
from mpl_toolkits import mplot3d

class CoarseGrainedFilters:
    """A class to specify helper functions to filter rows of
    the dataframe when creating a structure model"""
    @staticmethod
    def alpha_carbon_filter():
        """Simplest filter that only keeps the alpha carbon for each aminoacid"""
        return lambda row: row["full_id"][4][0] == "CA"

class StructureModel:
    """StructureModel is a 3D representation of a protein. It has different
    tools to model and visualize the surface.

    Parameters
    ----------
    point_ids : list or list-like
        Identifier for each item in points. Useful for understanding the
        results of the structure model.
    points : list, np.ndarray or pd.Series
        List or list-like of lists in which each element has three coordinates.
        In case of providing a pd.Series, the elements must be np.ndarrays.

    Attributes
    ----------
    point_ids: list-like
        A list-like structure for identifying each coordinate row (see `points`).
    points : np.ndarray
        Point coordinate matrix.
    simplices_order : np.ndarray
        Order of simplices as computed by Delaunay
    persistent_hom_params : dict
        Dict to store the result of persistent homology parameters
    delaunay : scipy.spatial.Delaunay
        Delaunay object used for triangulation of points.
    """

    def __init__(self, point_ids, points):
        self.point_ids = point_ids
        self.points = points
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
        if isinstance(points, pd.Series):
            if not isinstance(points[0], np.ndarray):
                raise Exception("If providing a pd.Series as points the elements must be numpy arrays.")
            self.__points = np.stack(points.values)
        elif isinstance(points, np.ndarray):
            self.__points = points
        elif isinstance(points, list):
            self.__points = np.stack(np.array(point) for point in points)
        else:
            raise Exception("points must be a list, a numpy matrix or a pd.Series of numpy arrays")

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
                points_diff = self.points[self.delaunay.simplices][:, pair_indices[0]] - self.points[self.delaunay.simplices][:, pair_indices[1]]
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
        return np.linalg.norm(self.delaunay.points[e[0]] - self.delaunay.points[e[1]])

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
        numpy.Array
            Array of simplices.

        """
        # TODO: check return DockBlock
        return self.delaunay.simplices[self.simplices_order[:step]]

    def _step_to_length(self, i):
        # TODO: check
        return max(self.edge_size(self.delaunay.simplices[self.simplices_order[i]]))

    def get_step_for_topology(self, topology=[1, 0, 0], recalculate=False):
        """Get step of the filtration to match the desired topology while
        having all the points. Result is cached for performance, though
        it can be recalculated with recalculate=True.
        """
        if str(topology) in self.persistent_hom_params and not recalculate:
            return self.persistent_hom_params[str(topology)]

        persistent_homology = PersistentHomology(self.delaunay.simplices[self.simplices_order])

        step = persistent_homology.get_step_for_topology(
            min_step=self._get_min_steps(self.simplices_order),
            betti_numbers=topology
        )

        self.persistent_hom_params[str(topology)] = step
        return step


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

    def plot(self, step, view_position, view_init_elev=20, view_init_azim=35,
             cmap="viridis"):
        """Plot surface of structure.

        Parameters
        ----------
        step : int
            step in the Delaunay triangulation to define surface. To plot the
            surface the simplices until that step will be used.
        view_position : int
            Used predefined position to view the surface. The following
            positions are available:
        view_init_elev : int
            If you want a personalized view position you can set the elevation
            using this parameter, it is then passed to matplotlib3d api,
            view_init function. First, you need to set view_position = 0.
        view_init_azim : int
            If you want a personalized view position you can set the azimuth
            using this parameter, it is then passed to matplotlib3d api,
            view_init function. First, you need to set view_position = 0.
        cmap : str
            `cmap` used to color the plot.

        Returns
        -------
        None
        """
        view_init_elev, view_init_azim = 35, view_position * 45

        triangles = self.get_simplices_by_step(step)
        ax = plt.axes(projection='3d')
        ax.view_init(view_init_elev, view_init_azim)
        ax.plot_trisurf(self.points[:, 0], self.points[:, 1],
                        self.points[:, 2], triangles=triangles, cmap=cmap)


class PersistentHomology:
    """Class that for a given filtration finds the cut to match a target
    topology."""
    def __init__(self, simplices):
        self.simplices = simplices

    def get_step_for_topology(self, min_step=0, betti_numbers=[1, 0, 0]):
        """
        Using gudhi, obtain the minimum step of the filtration that matches
        the topology that is higher than `min_step`.
        Parameters
        ----------
        min_step: int
            Default 0. The minimum step to use. It is recommended to use
            the minimum step in which all the desired nodes are in the
            simplicial complex.

        betti_numbers: list of int
            List like (b0, b1, b2) where:
            b0 is the first betti number, the number of connected components.
            b1 is the second betti number, is the number of holes.
            b2 is the third betti number, is the number of cavities.

            At each step of the simplicial tessellation betti numbers are
            computed using Gudhi. This parameter regulates when to save
            steps of interest.
            For example, if b0 is 1, then when there is a single connected
            component.

        Returns
        -------
        """
        simplex_tree = gudhi.SimplexTree()
        for i, tetrahedron in enumerate(self.simplices):
            simplex_tree.insert(tetrahedron, i)

        simplex_tree.compute_persistence()
        for birth in range(min_step, self.simplices.shape[0]):
            if simplex_tree.persistent_betti_numbers(birth, birth + 1) == betti_numbers:
                return birth

        raise ValueError("Desired topology could not be reached.")
