

import fedoo as fd
import numpy as np

fd.ModelingSpace("2D")

mesh_rect = fd.mesh.rectangle_mesh(nx=11, ny=21,
                                   x_min=0, x_max=1, y_min=0, y_max=1,
                                   elm_type='quad4', name='Domain')
                                   
mesh_disk = fd.mesh.disk_mesh(radius=0.5, nr=6, nt=6, elm_type='quad4')
mesh_disk.nodes += np.array([1.5, 0.48])
mesh = fd.Mesh.stack(mesh_rect, mesh_disk)

surf = fd.mesh.extract_surface(mesh)
ipc_contact = fd.constraint.IPCContact(
    mesh, surface_mesh=surf,
    friction_coefficient=0.3,     
    use_ccd=True,                
)

material = fd.constitutivelaw.ElasticIsotrop(200e3, 0.3)
wf = fd.weakform.StressEquilibrium(material, nlgeom=True)
solid_assembly = fd.Assembly.create(wf, mesh)
assembly = fd.Assembly.sum(solid_assembly, ipc_contact)

pb = fd.problem.NonLinear(assembly)
res = pb.add_output('results', solid_assembly, ['Disp', 'Stress'])
# ... add BCs ...
pb.nlsolve(dt=0.005, tmax=1)

