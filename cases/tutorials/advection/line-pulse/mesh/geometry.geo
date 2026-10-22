Point(1) = {0.0, 0.0, 0.0, 0.025};
Point(2) = {0.17, 0.0, 0.0, 0.035};
Point(3) = {0.43, 0.0, 0.0, 0.020};
Point(4) = {0.71, 0.0, 0.0, 0.040};
Point(5) = {1.0, 0.0, 0.0, 0.025};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};

Physical Curve("fluid") = {1, 2, 3, 4};
Physical Point("left") = {1};
Physical Point("right") = {5};
