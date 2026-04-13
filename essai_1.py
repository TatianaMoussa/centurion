import fedoo as fd
import numpy as np

fd.ModelingSpace("3D")
print('mon code marche')
Rext = 15.0
Rnuc = 6.0
H = 10.0

profil = fd.mesh.disk_mesh(
    radius=Rext,
    nr=12,
    nt=36,
    elm_type="quad4",
    ndim=None,
    name="disk2d"
)

mesh = fd.mesh.extrude(profil, H, 8)

conn = mesh.elements
nodes = mesh.nodes

centers = nodes[conn].mean(axis=1)
r = np.sqrt(centers[:, 0]**2 + centers[:, 1]**2)
print("hiiii")
elem_nucleus = np.where(r <= Rnuc)[0]
elem_annulus = np.where(r > Rnuc)[0]

mesh_nucleus = mesh.extract_elements(elem_nucleus, name="nucleus")
mesh_annulus = mesh.extract_elements(elem_annulus, name="annulus")

material_nucleus = fd.constitutivelaw.ElasticIsotrop(1.0, 0.49)
material_annulus = fd.constitutivelaw.ElasticIsotrop(4.0, 0.45)

wf_nucleus = fd.weakform.StressEquilibrium(material_nucleus)
wf_annulus = fd.weakform.StressEquilibrium(material_annulus)

assembly_nucleus = fd.Assembly.create(wf_nucleus, mesh_nucleus)
assembly_annulus = fd.Assembly.create(wf_annulus, mesh_annulus)

assembly = assembly_nucleus + assembly_annulus

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

results_nucleus = pb.get_results(assembly_nucleus,output_list=["Stress", "Disp", "Strain"])

results_annulus = pb.get_results(assembly_annulus,output_list=["Stress", "Disp", "Strain"])

plotter = results_nucleus.plot("stress", "vm",show=False)
results_annulus.plot("stress", "vm",plotter=plotter)

plotter.show()