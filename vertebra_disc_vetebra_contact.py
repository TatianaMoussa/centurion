import fedoo as fd
import numpy as np


fd.ModelingSpace("2D")


mesh_disc = fd.mesh.rectangle_mesh(
   nx=21, ny=6,
   x_min=-20, x_max=20,
   y_min=0, y_max=10,
   elm_type="quad4",
   name="disc"
)


mesh_vert = fd.mesh.rectangle_mesh(
   nx=21, ny=11,
   x_min=-20, x_max=20,
   y_min=10, y_max=30,
   elm_type="quad4",
   name="vertebre"
)


mesh = fd.Mesh.stack(mesh_disc, mesh_vert)


pairs = mesh.find_coincident_nodes(tol=1e-8)
if len(pairs) > 0:
   mesh.merge_nodes(pairs)


mesh.element_sets["disc"] = np.arange(0, mesh_disc.n_elements)
mesh.element_sets["vert"] = np.arange(mesh_disc.n_elements, mesh_disc.n_elements + mesh_vert.n_elements)


mat_disc = fd.constitutivelaw.ElasticIsotrop(1.5, 0.48, name="disc_mat")
mat_vert = fd.constitutivelaw.ElasticIsotrop(12000.0, 0.30, name="vert_mat")


material = fd.constitutivelaw.Heterogeneous(
   (mat_disc, mat_vert),
   ("disc", "vert")
)


wf = fd.weakform.StressEquilibrium(material)
assembly = fd.Assembly.create(wf, mesh)


pb = fd.problem.Linear(assembly)
pb.set_solver("cg")


bottom = mesh.find_nodes("Y", 0.0, tol=1e-8)
top = mesh.find_nodes("Y", 30.0, tol=1e-8)


pb.bc.add("Dirichlet", bottom, "Disp", 0.0)
pb.bc.add("Dirichlet", top, "DispY", -1.5)


pb.apply_boundary_conditions()
pb.solve()


results = pb.get_results(assembly, output_list=["Disp", "Stress", "Strain"])


results.plot("Stress", "vm", "Node", show=True, show_edges=False, show_nodes=False)
