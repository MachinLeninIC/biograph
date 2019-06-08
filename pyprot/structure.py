import numpy as np
from scipy.spatial import Delaunay
from itertools import combinations, chain
import os
import subprocess
import networkx as nx
from pyprot.io import Writer
import pyprot
import matplotlib.pyplot as plt
import pandas as pd


class StructureModel:
    """StructureModel is a 3D representation of a protein. It has different
    tools to model and visualize the surface.

    Parameters
    ----------
    points : list or list-like
        List or list-like of lists in which each element has three coordinates.

    Attributes
    ----------
    points : list or list-like
        List or list-like of lists in which each element has three coordinates.
    simplices_order : np.ndarray
        Order of simplices as computed by Delaunay
    persistent_hom_params : dict
        Dict to store the result of persistent homology parameters
    delaunay : scipy.spatial.Delaunay
        Delaunay object used for triangulation of points.
    """

    def __init__(self, points):
        self.__points = points
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
        if view_position == 0:
            None
        elif view_position == 1:
            view_init_elev, view_init_azim = 35, 0
        elif view_position == 2:
            view_init_elev, view_init_azim = 35, 45
        elif view_position == 3:
            view_init_elev, view_init_azim = 35, 90
        elif view_position == 4:
            view_init_elev, view_init_azim = 35, 135
        elif view_position == 5:
            view_init_elev, view_init_azim = 35, 180
        elif view_position == 6:
            view_init_elev, view_init_azim = 35, 225
        elif view_position == 7:
            view_init_elev, view_init_azim = 35, 270
        elif view_position == 8:
            view_init_elev, view_init_azim = 35, 315
        else:
            raise Exception("View not implemented")
        triangles = self.get_simplices_by_step(step)
        ax = plt.axes(projection='3d')
        ax.view_init(view_init_elev, view_init_azim)
        ax.plot_trisurf(self.points[:, 0], self.points[:, 1],
                        self.points[:, 2], triangles=triangles, cmap=cmap)


class Perseus:
    """ Class that wraps the usage of Perseus command line interface. With
    Perseus a persistent homology-based approach is used to model the surface
    and core of proteins. In addition, this can be used to model de depth of
    a given residue in the protein.
    """

    def __init__(self):
        MODULEDIR = os.path.dirname(os.path.abspath(__file__))
        self.__perseuspath = os.path.join(MODULEDIR, "perseus", "perseus")

    def execute_persistent_hom(self, structure, topology=(1, 0, 0),
                               perseus_input_filename="ordered_simplices.pers",
                               output_dir="perseus", output_filename="pdb"):
        """
        Function used to execute persistent homology. This function uses Perseus's
        simtop method.
        Perseus's documentation can be found in:
        http://people.maths.ox.ac.uk/nanda/perseus/index.html

        Parameters
        ----------
        topology : tuple or list
            topology is used to define different steps of interest.
            First element, is the first betti number (b0), is the number of
            connected components.
            Second element, is the second betti number (b1), is the number of
            holes.
            Third element, is the third betti number (b2), is the number of
            cavities.
            At each step of the simplicial tessellation bettin numbers are
            computed, using Perseus. This parameter regulates when to save
            steps of interest.
            For example, if b0 is 1, then when there is a connected component

        output_dir : str
            directory where output will be written
        output_filename : str
            name for output file

        Returns
        -------
        StructureModel
            StructureModel with persistent_hom_params computed
            StructureModel.persistent_hom_params is a dict with the following
            steps computed:
            b0_step: step at which all points in the structure are covered by
            the simplices.
            b1_step: when all points in the structure are covered by the
            simplices and topology[0] condition is met.
            Typically, it might be of interest the step at which there is only
            one connected component.
            b2_step: when all points in the structure are covered by the
            simplices and topology[0] and topology[1] conditions are met.
            For example, you might be interested in knowing when there is only
            one connected component and no holes.
            b3_step: when all points in the structure are covered by the
            simplices and topology[0] and topology[1] and topology[2] conditions
            are met.
            It may be the case that you are interested in knowing when there
            is only one connected component, no holes and no cavities."""
        structure_ = structure
        if isinstance(structure_, pyprot.protein.Protein):
            structure = structure.structure

        # TODO: validar que structure tenga datos
        # TODO: validar que sea Protein

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
        b0_step = structure._get_min_steps(order=structure.simplices_order)
        to_run = (self.__perseuspath, 'simtop',
                  os.path.join(dir_path, perseus_input_filename),
                  os.path.join(dir_path, "perseus_output"))
        popen = subprocess.Popen(to_run, stdout=subprocess.PIPE)
        popen.wait()
        with open(os.path.join(output_dir, output_filename,
                               'perseus_output_betti.txt'), 'r') as ifile:
            b1_step = b0_step
            b2_step = b0_step
            b3_step = b0_step
            for line in ifile.readlines()[1:]:
                # first line is blank
                i, h0, h1, h2, _ = [int(x) for x in line.split()]
                if i >= b0_step and h0 == topology[0]:

                    # only modifies b0_step once
                    if b1_step == b0_step:
                        b1_step = i

                    # only modifies the b2_step online once
                    if b2_step == b1_step and h1 == topology[1]:
                        b2_step = i
                    # save b3_step and break
                    if h1 == topology[1] and h2 == topology[2]:
                        b3_step = i
                        break

        if isinstance(structure_, pyprot.protein.Protein):
            structure_.structure.persistent_hom_params["b0_step"] = b0_step
            structure_.structure.persistent_hom_params["b1_step"] = b1_step
            structure_.structure.persistent_hom_params["b2_step"] = b2_step
            structure_.structure.persistent_hom_params["b3_step"] = b3_step
        else:
            structure_.persistent_hom_params["b0_step"] = b0_step
            structure_.persistent_hom_params["b1_step"] = b1_step
            structure_.persistent_hom_params["b2_step"] = b2_step
            structure_.persistent_hom_params["b3_step"] = b3_step
        return structure_

    @staticmethod
    def read_betti(betti_filename, step):
        """Read a betti file, the output of Perseus. If you had run Perseus
        before you can just read that result, instead of running the process
        again.

        Parameters
        ----------
        betti_filename : str
            full path to betti file
        step : int
            step until which read the file

        Returns
        -------
        pandas DataFrame
            Output of Perseus' simtop. For example, a line in the file can be
            12 14 4 7 0
            which indicates that when all the cells with birth time less than
            or equal to 12 are included, then there are 14 connected components
            4 tunnels, 7 cavities and no higher dimensional generators of
            homology. The numbers 14, 4 and 7 in this context are called the
            zeroth, first and second Betti numbers of the 12-th subcomplex in
            the persistence filtration.
        """
        ifile = open(betti_filename)
        betti = [list(map(int, x.split())) for x in ifile.readlines() if
                 x.split()]
        ifile.close()
        betti = [x for x in betti if int(x[0]) <= step]
        betti = pd.DataFrame(betti, columns=["birth_time",
                                             "0th_dim(connected_components)",
                                             "1th_dim(tunnels)",
                                             "2th_dim(cavities)",
                                             "3th_dim"])
        return betti
