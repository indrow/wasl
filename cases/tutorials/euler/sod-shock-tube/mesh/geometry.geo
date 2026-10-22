Point(1) = {0.0, 0.0, 0.0, 0.01};
Point(2) = {1.0, 0.0, 0.0, 0.01};

Line(1) = {1, 2};

Physical Curve("fluid") = {1};
Physical Point("left") = {1};
Physical Point("right") = {2};
