import fedoo as fd
import numpy as np

fd.ModelingSpace("2D")

L_disc = 20.0
H_disc = 8.0
L_vert = 30.0
H_vert = 12.0

mesh_vert_bot = fd.mesh.rectangle_mesh(
    nx=31, ny=13,
    x_min=-L_vert/2, x_max=L_vert/2,
    y_min=0.0, y_max=H_vert,
    elm_type="quad4",
    name="vertebre_bas"
)

mesh_disc = fd.mesh.rectangle_mesh(
    nx=21, ny=9,
    x_min=-L_disc/2, x_max=L_disc/2,
    y_min=H_vert, y_max=H_vert + H_disc,
    elm_type="quad4",
    name="disc"
)

mesh_vert_top = fd.mesh.rectangle_mesh(
    nx=31, ny=13,
    x_min=-L_vert/2, x_max=L_vert/2,
    y_min=H_vert + H_disc, y_max=2*H_vert + H_disc,
    elm_type="quad4",
    name="vertebre_haut"
)

mesh_1 = fd.Mesh.stack(mesh_vert_bot, mesh_disc)
mesh = fd.Mesh.stack(mesh_1 , mesh_vert_top)

pairs = mesh.find_coincident_nodes(tol=1e-8)
if len(pairs) > 0:
    mesh.merge_nodes(pairs)

n_bot = mesh_vert_bot.n_elements
n_disc = mesh_disc.n_elements
n_top = mesh_vert_top.n_elements
mesh.element_sets["vert_bot"] = np.arange(0, n_bot)
mesh.element_sets["disc"] = np.arange(n_bot, n_bot + n_disc)
mesh.element_sets["vert_top"] = np.arange(n_bot + n_disc, n_bot + n_disc + n_top)

mat_vert = fd.constitutivelaw.ElasticIsotrop(12000.0, 0.30, name="vert_mat")
mat_disc = fd.constitutivelaw.ElasticIsotrop(1.5, 0.48, name="disc_mat")

material = fd.constitutivelaw.Heterogeneous(
    (mat_vert, mat_disc, mat_vert),
    ("vert_bot", "disc", "vert_top")
)

wf = fd.weakform.StressEquilibrium(material)
assembly = fd.Assembly.create(wf, mesh)

pb = fd.problem.Linear(assembly)
pb.set_solver("cg")

bottom = mesh.find_nodes("Y", 0.0, tol=1e-8)
top = mesh.find_nodes("Y", 2*H_vert + H_disc, tol=1e-8)

pb.bc.add("Dirichlet", bottom, "Disp", 0.0)
pb.bc.add("Dirichlet", top, "DispY", -1.2)

pb.apply_boundary_conditions()
pb.solve()

results = pb.get_results(assembly, output_list=["Disp", "Stress", "Strain"])
results.plot(
    "Stress",
    "vm",
    "Node",
    show=True,
    show_nodes=False,
    show_edges=False
)
