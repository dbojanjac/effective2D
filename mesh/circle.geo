/*********************************************************************
 *
 * CIRCLE PATCH unit cell in 2D by Dario Bojanjac
 *
 * Subdomains in Unit Cell mesh are annotated by:
 *    1: Circle Patch with radius r inside unit cell (colored red)
 *    2: Outside material 1^2 (colored blue)
 *
 *  Input parameters:
 *    r: circle patch radius
 *
 * GMSH 3.0.6 version
 *********************************************************************/

r = 0.3;

// Meshing parameters
lc_patch = 10 * 1e-3;
lc_centre = 30 * 1e-3;
lc_outer = 30 * 1e-3;

// Circle Patch
Point(1) = {0.5 + r, 0.5, 0, lc_patch};
Point(2) = {0.5, 0.5 + r, 0, lc_patch};
Point(3) = {0.5 - r, 0.5, 0, lc_patch};
Point(4) = {0.5, 0.5 - r, 0, lc_patch};
Point(5) = {0.5, 0.5, 0, lc_centre};

Circle(1) = {1, 5, 2};
Circle(2) = {2, 5, 3};
Circle(3) = {3, 5, 4};
Circle(4) = {4, 5, 1};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Point{5} In Surface{1};
Physical Surface(1) = {1};

Color Red{ Surface{ 1 }; }
//------------------------------------------------------------------------------


// Unite cell
Point(11) = {0, 0, 0, lc_outer};
Point(12) = {1, 0, 0, lc_outer};
Point(13) = {1, 1, 0, lc_outer};
Point(14) = {0, 1, 0, lc_outer};

Line(11) = {11, 12};
Line(12) = {12, 13};
Line(13) = {13, 14};
Line(14) = {14, 11};

Line Loop(2) = {11, 12, 13, 14};
Plane Surface(2) = {2, -1};
Physical Surface(2) = {2};

Color Blue{ Surface{ 2 }; }

Mesh.MshFileVersion = 2;

