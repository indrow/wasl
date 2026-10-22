Point(1) = {0.0, 0.0, 0.0, 0.025};
Point(2) = {1.0, 0.0, 0.0, 0.025};

Line(1) = {1, 2};

Physical Curve("fluid") = {1};
Physical Point("left") = {1};
Physical Point("right") = {2};
