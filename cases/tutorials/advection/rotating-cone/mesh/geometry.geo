SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 1, 1, 0};
Mesh.Algorithm = 6;
Mesh.MeshSizeMin = 0.05;
Mesh.MeshSizeMax = 0.05;
Physical Surface("fluid") = {1};
Physical Curve("domain") = {1, 2, 3, 4};
