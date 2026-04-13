SetFactory("OpenCASCADE");

R_nucleus = 8.0;   
R_outer   = 18.0;  
H         = 10.0;  
lc = 2.0;

Cylinder(1) = {0, 0, 0,   0, 0, H,   R_nucleus};

Cylinder(2) = {0, 0, 0,   0, 0, H,   R_outer};

BooleanDifference(3) = { Volume{2}; Delete; }{ Volume{1}; };

Cylinder(4) = {0, 0, 0, 0, 0, H,   R_nucleus};

BooleanFragments{ Volume{3}; Delete; }{ Volume{4}; Delete; }

Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;

Physical Volume("Nucleus") = {4};
Physical Volume("Annulus") = {3};

Mesh 3;
Save "disque_2zones.msh";











