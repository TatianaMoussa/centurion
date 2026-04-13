import fedoo as fd
import numpy as np

profil = fd.mesh.disk_mesh( radius=15.0,
    nr=12,
    nt=36,
    elm_type="quad4",
    ndim=None,
    name="disk2d"
)
mesh = fd.mesh.extrude(profil, 10.0, 8)
print(f"element type: '{mesh.elm_type}'")

fd.ModelingSpace("3D")

material = fd.constitutivelaw.ElasticIsotrop(2.0, 0.45)
wf = fd.weakform.StressEquilibrium(material)
assembly = fd.Assembly.create(wf, mesh)


pb = fd.problem.Linear(assembly)
pb.set_solver("cg")

zmin = mesh.bounding_box.zmin
zmax = mesh.bounding_box.zmax

bottom = mesh.find_nodes("Z", zmin, tol=1e-8)
top = mesh.find_nodes("Z", zmax, tol=1e-8)

pb.bc.add("Dirichlet", bottom, "Disp", 0.0)
pb.bc.add("Dirichlet", top, "DispZ", -0.8)

pb.apply_boundary_conditions()
pb.solve()

res = pb.get_results(assembly, ["Stress", "Disp"])
res.plot("Disp", "Z", "Node", show_edges=False, scale=2)