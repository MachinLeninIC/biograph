### Acá dejo pedazos de código que, por el momento, no los incluímos en el módulo principal. Si eventualmente abrimos
# este paquete deberíamos volarlo...

def calculate_core(self, verbose=True):
    if self.delaunay is None:
        self.initiate_delaunay_()
        self.set_simplices_order()

    boundary, triangles = self.get_boundary_matrix()
    removed_simplices = np.zeros(self.delaunay.nsimplex, dtype='int')
    vertex_occupation = {v: [s for s, simplex in enumerate(self.delaunay.simplices) if v in simplex] for v in
                         range(self.delaunay.npoints)}
    boundary_non_zero = np.array(boundary.nonzero())
    triangles_count = np.bincount(boundary_non_zero[0])
    external = boundary_non_zero[1][triangles_count[boundary_non_zero[0]] == 1]
    external_available = external[removed_simplices[external] == 0]
    surface_vertices = list(np.unique(triangles[triangles_count == 1].flatten()))
    surface_edges = {}  # e: 0 for e in combinations(range(self.npoints), 2)}
    for i, t in enumerate(triangles):  # [triangles_count == 1]:
        for e in self.faces_(t):
            surface_edges[e] = 1 if triangles_count[i] == 1 else 0

    second_time = False
    while external_available.size != 0:
        k = self.simplices_order[external_available].argmax()  # chooses bigest tetrahedron (I hope)
        s = external_available[k]
        simplex = self.delaunay.simplices[s]
        min_vertex_occupation = min([len(vertex_occupation[v]) for v in simplex])
        boundary_of_simplex = [t for t, ss in boundary_non_zero.T if ss == s]
        in_surface_triangles = triangles[[t for t in boundary_of_simplex if triangles_count[t] == 1]]
        in_surface_vertices = list(np.unique(in_surface_triangles.flatten()))
        insurface_edges = [e for e in combinations(simplex, 2) if surface_edges[e] == 1]
        is_removable = \
            (len(in_surface_triangles) == 1 and len(insurface_edges) == 3 and len(in_surface_vertices) == 3) \
            or (len(in_surface_triangles) == 2 and len(insurface_edges) == 5) \
            or (len(
                in_surface_triangles) == 3 and min_vertex_occupation > 1)
        if is_removable:
            removed_simplices[s] = -1
            for t in boundary_of_simplex:
                triangles_count[t] = triangles_count[t] - 1
            for v in simplex:
                vertex_occupation[v].remove(s)
            external = boundary_non_zero[1][triangles_count[boundary_non_zero[0]] == 1]
            for e in combinations(simplex, 2):
                surface_edges[e] = 1
            second_time = False  # if a simplex is removed, it is not a dummy step
            removed_simplices = np.minimum(removed_simplices, np.zeros(removed_simplices.size,
                                                                       dtype=removed_simplices.dtype))  # forgets fixed simplices
        else:
            removed_simplices[s] = 1  # simplex skiped
        external_available = external[removed_simplices[external] == 0]

    if verbose:
        print(('done, ended with ', len(removed_simplices[removed_simplices == -1]), ' cells removed and ',
               len(removed_simplices[removed_simplices == 1]), ' cells visited but fixed'))
    return removed_simplices  # , surfaceVertices, external_available



    def get_boundary_matrix(self, simplices=None):
        if simplices is None and self.delaunay is None:
            self.initiate_delaunay_()
            self.set_simplices_order()
        simplices = simplices if type(simplices) != type(None) else self.delaunay.simplices
        m = simplices.shape[0]
        d = simplices[0].shape[0]-1
        if d == 1:
            faces = np.zeros(self.delaunay.npoints, dtype='int32')
            n=self.delaunay.npoints
            boundary = dok_matrix((n,m), dtype='int32')
            for j,s in enumerate(simplices):
                boundary[s[1], j] =  1
                boundary[s[0], j] = -1
        else:
            b = {}
            for j,s in enumerate(simplices):
                sign = 1
                for t in self.faces_(s):
                    if t in b:
                        b[t].append((j,sign))
                    else:
                        b[t]  = [(j,sign)]
                    sign = -sign
            faces = np.array(list(b.keys()), dtype='int32')
            n = faces.shape[0]
            boundary=dok_matrix((n,m), dtype='int32')
            for i,t in enumerate(b.keys()):
                for j,s in b[t]:
                    boundary[i,j] = s
        return boundary, faces

