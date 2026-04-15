import fedoo as fd
import numpy as np

fd.ModelingSpace("2D")

filename = 'disk_rectangle_contact' #file to save results

#---- Create geometries --------------
mesh_rect = fd.mesh.rectangle_mesh(nx=11, ny=21,
                                   x_min=0, x_max=1, y_min=0, y_max=1,
                                   elm_type = 'quad4', name = 'Domain'
                                  )
mesh_rect.element_sets['rect'] = np.arange(0, mesh_rect.n_elements)
mesh_disk = fd.mesh.disk_mesh(radius=0.5, nr=6, nt=6, elm_type = 'quad4')
mesh_disk.nodes+=np.array([1.5,0.48]) # translate disk on the right
mesh_disk.element_sets['disk'] = np.arange(0,mesh_disk.n_elements)

mesh = fd.Mesh.stack(mesh_rect,mesh_disk)

#node sets for boundary conditions
nodes_left = mesh.find_nodes('X',0)
nodes_right = mesh.find_nodes('X',1)

nodes_bc = mesh.find_nodes('X>1.5')
nodes_bc = list(set(nodes_bc).intersection(mesh.node_sets['boundary']))

#---- Define contact --------------
#slave surface = right face of rectangle mesh
nodes_contact = nodes_right
surf = fd.mesh.extract_surface(mesh.extract_elements('disk'))
contact = fd.constraint.Contact(nodes_contact, surf)
contact.contact_search_once = True
contact.eps_n = 5e5

#---- Material properties --------------
props = np.array([200e3, 0.3, 1e-5, 300, 1000, 0.3])
# E, nu, alpha (non used), Re, k, m
material_rect = fd.constitutivelaw.Simcoon("EPICP", props)
material_disk = fd.constitutivelaw.ElasticIsotrop(50e3, 0.3) #E, nu
material = fd.constitutivelaw.Heterogeneous(
    (material_rect, material_disk),
    ('rect', 'disk')
    )

#---- Build problem --------------
wf = fd.weakform.StressEquilibrium(material, nlgeom = True)
solid_assembly = fd.Assembly.create(wf, mesh)

assembly = fd.Assembly.sum(solid_assembly, contact)

pb = fd.problem.NonLinear(assembly)
results = pb.add_output(filename,
                        solid_assembly,
                        ['Disp', 'Stress', 'Strain', 'Fext']
                        )

pb.bc.add('Dirichlet',nodes_left, 'Disp',0)
pb.bc.add('Dirichlet',nodes_bc, 'Disp', [-0.4,0.2])

pb.set_nr_criterion("Displacement", tol = 5e-3, max_subiter = 5)

#---- Solve problem in two steps: load, unload --------------
pb.nlsolve(dt = 0.005, tmax = 1, update_dt = True, interval_output = 0.01)

pb.bc.remove(-1) #remove last boundary contidion
pb.bc.add('Dirichlet',nodes_bc, 'Disp', [0,0])

pb.nlsolve(dt = 0.005, tmax = 1, update_dt = True, interval_output = 0.01)


# =============================================================
# Example of plots with pyvista - uncomment the desired plot
# =============================================================

# ------------------------------------
# Simple plot with default options
# ------------------------------------
results.plot('Stress', component='vm', data_type='Node')
results.plot('Disp', component = 0, data_type='Node')

# ------------------------------------
# Write movie with default options
# ------------------------------------
results.write_movie(filename,
                    'Stress',
                    component = 'XX',
                    data_type = 'Node',
                    framerate = 24,
                    quality = 5,
                    clim = [-3e3, 3e3]
                   )