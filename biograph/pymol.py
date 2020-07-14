import pymol
from pymol import cmd
from pymol.wizard import Wizard
from chempy import cpv
from pymol.cgo import *
import sys

def makePrimitive(cgo, name):
    az = cmd.get('auto_zoom', quiet=1)
    cmd.set('auto_zoom', 0, quiet=1)
    cmd.load_cgo(cgo, name)
    cmd.set('auto_zoom', az, quiet=1)

def point(p):
    x, y, z = p
    return [COLOR, 1, 1, 1, SPHERE, float(x), float(y), float(z), 0.5]

def line(p1, p2):
    x1, y1, z1 = p1
    x2, y2, z2 = p2
    return [CYLINDER, float(x1), float(y1), float(z1), float(x2), float(y2), float(z2), 0.25, 1, 1, 1, 1, 1, 1]

def triangle(corner1, corner2, corner3, normal):
    planeObj = []
    planeObj.extend(point(corner1))
    planeObj.extend(point(corner2))
    planeObj.extend(point(corner3))
    planeObj.extend(line(corner1, corner2))
    planeObj.extend(line(corner2, corner3))
    planeObj.extend(line(corner3, corner1))

    planeObj.extend([COLOR, 0.8, 0.8, 0.8])
    planeObj.extend([BEGIN, TRIANGLE_STRIP])
    #planeObj.append(NORMAL)
    #planeObj.extend(normal)
    for corner in [corner1, corner2, corner3, corner1]:
        planeObj.append(VERTEX)
        planeObj.extend(corner)
    planeObj.append(END)
    return planeObj

def planeFromPoints(point1, point2, point3, facetSize):
    corner1 = point1
    corner2 = point2
    corner3 = point3
    normal = cpv.cross_product(corner1, corner2)
    return triangle(corner1, corner2, corner3, normal)


class PlaneWizard(Wizard):

    def __init__(self):
        Wizard.__init__(self)

        # some attributes to do with picking
        self.pick_count = 0
        self.object_count = 0
        self.object_prefix = "pw"

        # the plane facet size (the 'radius' of the section of plane we show)
        self.facetSize = 5

        self.selection_mode = cmd.get_setting_legacy("mouse_selection_mode")
        cmd.set("mouse_selection_mode",0) # set selection mode to atomic
        cmd.deselect()

    def draw_triangle(self, point1, point2, point3):
        plane = planeFromPoints(point1, point2, point3, self.facetSize)

        planeName = "plane-%02d" % self.object_count
        self.object_count += 1
        makePrimitive(plane, planeName)
        cmd.show("cgo", "plane*")

    def gudhi_topology(self, pdb_name):
        from biograph.protein import Protein
        #from biograph.structure import Perseus
        from itertools import combinations

        print(pdb_name)

        protein = Protein.fetch(pdb_name, base_path="/tmp")
        protein.df = protein.df[~protein.df.coord.isnull()]
        structure = protein.generate_structure(lambda row: row["full_id"][4][0] == "CA")
        b3_step = structure.get_step_for_topology(topology=[1,0,0])

        core = protein.structure.get_simplices_by_step(b3_step)

        for i, tetrahedron in enumerate(core):
            #time.sleep(1)
            for face in combinations(tetrahedron, 3):
                point1 = structure.points[face[0]]
                point2 = structure.points[face[1]]
                point3 = structure.points[face[2]]
                self.draw_triangle(point1, point2, point3)

    def reset(self):
        cmd.delete(self.object_prefix + "*")
        cmd.delete("sele*")
        cmd.delete("_indicate*")
        cmd.unpick()
        cmd.refresh_wizard()

    def delete_all(self):
        cmd.delete("plane*")

    def cleanup(self):
        cmd.set("mouse_selection_mode",self.selection_mode) # restore selection mode
        self.reset()
        self.delete_all()

    def get_prompt(self):
        self.prompt = None
        return self.prompt


    def get_panel(self):
        return [
            [ 1, 'Plane Wizard',''],
            [ 2, 'Delete All Planes' , 'cmd.get_wizard().delete_all()'],
        ]

# Running on PyMol
wiz = PlaneWizard()
# add arg
cmd.extend("delaunay", wiz.gudhi_topology)
# make this the active wizard
cmd.set_wizard(wiz)


