from Bio.PDB import PDBParser
from biograph.protein import Protein
from biograph.structure import StructureModel
from biograph.structure import Perseus
import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
% matplotlib inline
import gc
gc.collect()
parser = PDBParser()
pdb = parser.get_structure("1a3z","/home/usuario/Documentos/github/biograph/tests/1A3Z.pdb")
import pandas as pd
prot = Protein(pdb, generate_dataframe=False)
prot.head()
prot.df.head()
df = prot.generate_dataframe()
df


pd.DataFrame(full_atom)
"resname"
ch = [i for i in prot.pdb.get_chains()][0]
ch.id
ch = [i for i in prot.pdb.get_models()][0]
ch.id

model = [i for i in prot.pdb.get_models()][0]
[i for i in model.get_chains()]
atoms = prot.get_atoms(filter_atoms=lambda x: x["name"] == "CA")


atoms[0]
pd.DataFrame(list(atoms))
atoms = prot.get_atoms(filter_attr=lambda x: x.coord,filter_atoms=lambda x: x.name == "CA", atom_as_dict=False)
#atoms = [list(i.coord) for i in prot.pdb.get_atoms()]

# struc = StructureModel(atoms)
# struc

prot.generate_structure(atoms)

p = Perseus()
import os
cwd = os.getcwd()
path = cwd + "/perseus/pdb/perseus_output_betti.txt"
path
p.read_betti(path, 500)

prot.structure.plot(300,10)
prot.structure.plot(300,1, view_init_elev=35, view_init_azim=45)
prot.structure.plot(300,2, view_init_elev=35, view_init_azim=90)
prot.structure.plot(300,3, view_init_elev=35, view_init_azim=135)
prot.structure.plot(300,4, view_init_elev=35, view_init_azim=180)
prot.structure.plot(300,1, view_init_elev=35, view_init_azim=225)
prot.structure.plot(300,1, view_init_elev=35, view_init_azim=270)
prot.structure.plot(300,1, view_init_elev=35, view_init_azim=315)
prot.structure.plot(300,1, view_init_elev=35, view_init_azim=360)
prot.structure.plot(300,1, view_init_elev=35, view_init_azim=0)
per = Perseus()
sutr = per.calculate_fatcore(prot)

ca_atoms = atoms
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
