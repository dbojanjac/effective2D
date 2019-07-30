/*********************************************************************
 *
 * CROSS PATCH unit cell in 2D by Dario Bojanjac
 *
 * Subdomains in Unit Cell mesh are annotated by:
 *    1: Cross Patch inside unit cell (colored red)
 *    2: Outside material 1^2 (colored blue)
 *
 *  Input parameters:
 *    h: cross leg width
 *    l: cross leg length
 *
 * GMSH 3.0.6 version
 *********************************************************************/

h = 0.2;
l = 0.25;

// Meshing parameters
lc_patch = 15 * 1e-3;
lc_outer = 20 * 1e-3;

// Cross Patch
Point(1) = {0.5 - h / 2, 0.5 - h / 2, 0, lc_patch};
Point(2) = {0.5 - h / 2, 0.5 - l, 0, lc_patch};
Point(3) = {0.5 + h / 2, 0.5 - l, 0, lc_patch};
Point(4) = {0.5 + h / 2, 0.5 - h / 2, 0, lc_patch};

Point(5) = {0.5 + l, 0.5 - h / 2, 0, lc_patch};
Point(6) = {0.5 + l, 0.5 + h / 2, 0, lc_patch};
Point(7) = {0.5 + h / 2, 0.5 + h / 2, 0, lc_patch};

Point(8) = {0.5 + h / 2, 0.5 + l, 0, lc_patch};
Point(9) = {0.5 - h / 2, 0.5 + l, 0, lc_patch};
Point(10) = {0.5 - h / 2, 0.5 + h / 2, 0, lc_patch};

Point(11) = {0.5 - l, 0.5 + h / 2, 0, lc_patch};
Point(12) = {0.5 - l, 0.5 - h / 2, 0, lc_patch};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 1};

Line Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
Plane Surface(1) = {1};
Physical Surface(1) = {1};

Color Red{ Surface{ 1 }; }
//------------------------------------------------------------------------------


// Unite cell
Point(101) = {0, 0, 0, lc_outer};
Point(102) = {1, 0, 0, lc_outer};
Point(103) = {1, 1, 0, lc_outer};
Point(104) = {0, 1, 0, lc_outer};

Line(101) = {101, 102};
Line(102) = {102, 103};
Line(103) = {103, 104};
Line(104) = {104, 101};

Line Loop(100) = {101, 102, 103, 104};
Plane Surface(2) = {100, -1};
Physical Surface(2) = {2};

Color Blue{ Surface{ 2 }; }

Mesh.MshFileVersion = 2;

