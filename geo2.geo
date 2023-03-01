/*
For generaling quad elemets with gmsh go to
Tools -> options -> mesh -> general ->  Suvdivision algorithm -> all quads
*/
a =  3;
b =  5;
c =  7;
d =  10;
e = 15;
f = 25;
delta = 0.1;
//+
Point(1)  = {15.8, 1.2, 0.0, delta};
Point(2)  = {32.8, 1.2, 0.0, delta};
Point(3)  = {32.8, 7.3, 0.0, delta};
Point(4)  = {31.3, 7.3, 0.0, delta};
Point(5)  = {31, 8.5, 0.0, delta};
Point(6)  = {30, 8.5, 0.0, delta};
Point(7)  = {30.1, 7.3, 0.0, delta};
Point(8)  = {30.3, 6.1, 0.0, delta};
Point(9)  = {24.7, 5.8, 0.0, delta};
Point(10) = {24.7, 5.6, 0.0, delta};
Point(11) = {24.2, 5.6, 0.0, delta};
Point(12) = {24.2, 5.3, 0.0, delta};
Point(13) = {23.8, 5.3, 0.0, delta};
Point(14) = {23.8, 5.1, 0.0, delta};
Point(15) = {23.4, 5.1, 0.0, delta};
Point(16) = {23.4, 4.9, 0.0, delta};
Point(17) = {23, 4.9, 0.0, delta};
Point(18) = {23, 4.6, 0.0, delta};
Point(19) = {22.5, 4.6, 0.0, delta};
Point(20) = {22.5, 4.3, 0.0, delta};
Point(21) = {20.3, 4.3, 0.0, delta};
Point(22) = {20.3, 4, 0.0, delta};
Point(23) = {19.8, 4, 0.0, delta};
Point(24) = {19.8, 3.8, 0.0, delta};
Point(25) = {19.4, 3.8, 0.0, delta};
Point(26) = {19.4, 3.5, 0.0, delta};
Point(27) = {19, 3.5, 0.0, delta};
Point(28) = {19, 3.3, 0.0, delta};
Point(29) = {18.6, 3.3, 0.0, delta};
Point(30) = {18.6, 3.1, 0.0, delta};
Point(31) = {18.2, 3, 0.0, delta};
Point(32) = {18.2, 2.8, 0.0, delta};
Point(33) = {17.7, 2.8, 0.0, delta};
Point(34) = {17.7, 2.6, 0.0, delta};
Point(35) = {17.4, 2.5, 0.0, delta};
Point(36) = {17.4, 2.3, 0.0, delta};
Point(37) = {16.9, 2.3, 0.0, delta};
Point(38) = {16.9, 2, 0.0, delta};
Point(39) = {16.5, 2, 0.0, delta};
Point(40) = {16.5, 1.8, 0.0, delta};
Point(41) = {16.1, 1.8, 0.0, delta};
Point(42) = {16.1, 1.5, 0.0, delta};
Point(43) = {15.8, 1.5, 0.0, delta};
//+
Line(1)  = {1,  2};
Line(2)  = {2,  3};
Line(3)  = {3,  4};
Line(4)  = {4,  5};
Line(5)  = {5,  6};
Line(6)  = {6,  7};
Line(7)  = {7,  8};
Line(8)  = {8,  9};
Line(9)  = {9, 10};
Line(10) = {10,11};
Line(11) = {11,12};
Line(12) = {12,13};
Line(13) = {13,14};
Line(14) = {14,15};
Line(15) = {15,16};
Line(16) = {16,17};
Line(17) = {17,18};
Line(18) = {18,19};
Line(19) = {19,20};
Line(20) = {20,21};
Line(21) = {21,22};
Line(22) = {22,23};
Line(23) = {23,24};
Line(24) = {24,25};
Line(25) = {25,26};
Line(26) = {26,27};
Line(27) = {27,28};
Line(28) = {28,29};
Line(29) = {29,30};
Line(30) = {30,31};
Line(31) = {31,32};
Line(32) = {32,33};
Line(33) = {33,34};
Line(34) = {34,35};
Line(35) = {35,36};
Line(36) = {36,37};
Line(37) = {37,38};
Line(38) = {38,39};
Line(39) = {39,40};
Line(40) = {40,41};
Line(41) = {41,42};
Line(42) = {42,43};
Line(43) = {43,1};
//+
Curve Loop(10) = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43};
//+
Plane Surface(100) = {10};

// To generate quadrangles instead of triangles, we can simply add
Recombine Surface{100};
Transfinite Curve {1 } = f  Using Progression 1;
Transfinite Curve {2 } = f  Using Progression 1;
Transfinite Curve {3 } = c  Using Progression 1;
Transfinite Curve {4 } = c  Using Progression 1;
Transfinite Curve {5 } = c  Using Progression 1;
Transfinite Curve {6 } = c  Using Progression 1;
Transfinite Curve {7 } = c  Using Progression 1;
Transfinite Curve {8 } = e  Using Progression 1;
Transfinite Curve {9 } = a  Using Progression 1;
Transfinite Curve {10} = b  Using Progression 1;
Transfinite Curve {11} = a  Using Progression 1;
Transfinite Curve {12} = b  Using Progression 1;
Transfinite Curve {13} = a  Using Progression 1;
Transfinite Curve {14} = b  Using Progression 1;
Transfinite Curve {15} = a  Using Progression 1;
Transfinite Curve {16} = b  Using Progression 1;
Transfinite Curve {17} = a  Using Progression 1;
Transfinite Curve {18} = b  Using Progression 1;
Transfinite Curve {19} = a  Using Progression 1;
Transfinite Curve {20} = d  Using Progression 1;
Transfinite Curve {21} = a  Using Progression 1;
Transfinite Curve {22} = b  Using Progression 1;
Transfinite Curve {23} = a  Using Progression 1;
Transfinite Curve {24} = b  Using Progression 1;
Transfinite Curve {25} = a  Using Progression 1;
Transfinite Curve {26} = b  Using Progression 1;
Transfinite Curve {27} = a  Using Progression 1;
Transfinite Curve {28} = b  Using Progression 1;
Transfinite Curve {29} = a  Using Progression 1;
Transfinite Curve {30} = b  Using Progression 1;
Transfinite Curve {31} = a  Using Progression 1;
Transfinite Curve {32} = b  Using Progression 1;
Transfinite Curve {33} = a  Using Progression 1;
Transfinite Curve {34} = b  Using Progression 1;
Transfinite Curve {35} = a  Using Progression 1;
Transfinite Curve {36} = b  Using Progression 1;
Transfinite Curve {37} = a  Using Progression 1;
Transfinite Curve {38} = b  Using Progression 1;
Transfinite Curve {39} = a  Using Progression 1;
Transfinite Curve {40} = b  Using Progression 1;
Transfinite Curve {41} = a  Using Progression 1;
Transfinite Curve {42} = b  Using Progression 1;
Transfinite Curve {43} = a  Using Progression 1;
//Uncomment the following line to try
// the Frontal-Delaunay algorithms for quads:
//
//Mesh.Algorithm = 8;


//+
Transfinite Curve {1} = 20 Using Progression 1;
Transfinite Curve {2} = 20 Using Progression 1;
Transfinite Curve {3} = 5 Using Progression 1;
Transfinite Curve {4} = 5 Using Progression 1;
Transfinite Curve {5} = 5 Using Progression 1;
Transfinite Curve {6} = 5 Using Progression 1;
Transfinite Curve {7} = 5 Using Progression 1;
Transfinite Curve {8} = 10 Using Progression 1;
Transfinite Curve {9} = 2 Using Progression 1;
Transfinite Curve {10} = 4 Using Progression 1;


// The default recombination algorithm might leave some triangles in the mesh,
// if recombining all the triangles leads to badly shaped quads. In such cases,
// to generate full-quad meshes, you can either subdivide the resulting hybrid
// mesh (with Mesh.SubdivisionAlgorithm = 1), or use the full-quad recombination
// algorithm, which will automatically perform a coarser mesh followed by
// recombination, smoothing and subdivision. Uncomment the following line to try
// the full-quad algorithm:
//
//Mesh.RecombinationAlgorithm = 3; // or 3

Physical Surface(1000) = {100};
