from Bio.PDB import PDBParser
from pyprot.protein import Protein
from pyprot.structure import StructureModel
from pyprot.structure import Perseus
import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import gc
gc.collect()
parser = PDBParser()
pdb = parser.get_structure("1a3z","/home/usuario/Documentos/github/pyprot/tests/1A3Z.pdb")
type(pdb)

prot = Protein(pdb)


prot.get_residues(filter_attr=lambda x: x["resname"], res_as_dict=True, filter_res=lambda x: x["resname"] == "HOH")

class OurClass:

    def __init__(self, a):
        self.OurAtt = a

    @property
    def OurAtt(self):
        return self.__OurAtt

    @OurAtt.setter
    def OurAtt(self, val):
        if val < 0:
            self.__OurAtt = 0
        elif val > 1000:
            self.__OurAtt = 1000
        else:
            self.__OurAtt = val

x = OurClass(10000000)
x.OurAtt

atoms = prot.get_atoms(filter_attr=lambda x: x.coord,filter_atoms=lambda x: x.name == "CA", atom_as_dict=False)
#atoms = [list(i.coord) for i in prot.pdb.get_atoms()]

struc = StructureModel(atoms)




per = Perseus()
sutr = per.calculate_fatcore(prot)

# %%
x, y, z = ca_atoms[:,0], ca_atoms[:,1], ca_atoms[:,2]
ax = plt.axes(projection='3d')
ax.scatter3D(x,y,z, 'viridis')

# %%
ax = plt.axes(projection='3d')
ax.plot_trisurf(x,y,z, cmap="viridis")
# %%
triangles = struc.get_simplices_by_step(357)
plt.axes(projection='3d').plot_trisurf(x,y,z, triangles=triangles ,cmap="viridis")

# %%

triangles = struc.get_simplices_by_step(389)
plt.axes(projection='3d').plot_trisurf(x,y,z, triangles=triangles ,cmap="viridis")
# %%
triangles = struc.get_simplices_by_step(389)
ax = plt.axes(projection='3d')
ax.view_init(40, 35)
ax.plot_trisurf(x,y,z, triangles=triangles ,cmap="viridis")
# %%
triangles = struc.get_simplices_by_step(389)
ax = plt.axes(projection='3d')
ax.view_init(20, 35)
ax.plot_trisurf(x,y,z, triangles=triangles ,cmap="viridis")
# %%
triangles = struc.get_simplices_by_step(550)
ax = plt.axes(projection='3d')
ax.view_init(20, 35)
ax.plot_trisurf(x,y,z, triangles=triangles ,cmap="viridis")
