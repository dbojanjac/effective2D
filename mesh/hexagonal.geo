/*********************************************************************
 *
 *  HEXAGONAL square unit cell in 2D by Dario Bojanjac
 *
 *********************************************************************/

// Hexagon Parameters

 a = 0.2;
 del = 0.00;
 v = a * Sqrt( 3 ) / 2;
 h = 0.5 * (1 - 2 * a * Sqrt( 3 ));
 Printf("Razmak h iznosi %f", h);

// Meshing parameters
 lc1 = 3 * 1e-3;
 lc2 = lc1;

// Hexagon 1

 Point(11) = {0, h/2 + del, 0, lc1};
 Point(12) = {a/2, h/2 + del, 0, lc1};
 Point(13) = {a, h/2 + v, 0, lc1};
 Point(14) = {a/2, h/2 + 2*v - del, 0, lc1};
 Point(15) = {0, h/2 + 2*v - del, 0, lc1};

 Line(11) = {11, 12};
 Line(12) = {12, 13};
 Line(13) = {13, 14};
 Line(14) = {14, 15};
 Line(15) = {15, 11};

 Line Loop(1) = {11, 12, 13, 14, 15};

// Hexagon 2

 Point(21) = {0.5 - a, 0, 0, lc1};
 Point(22) = {0.5 + a, 0, 0, lc1};
 Point(23) = {0.5 + a/2, v - del, 0, lc1};
 Point(24) = {0.5 - a/2, v - del, 0, lc1};

 Line(21) = {21, 22};
 Line(22) = {22, 23};
 Line(23) = {23, 24};
 Line(24) = {24, 21};

 Line Loop(2) = {21, 22, 23, 24};

// Hexagon 3

 Point(31) = {1 - a/2, h/2 + del, 0, lc1};
 Point(32) = {1, h/2 + del, 0, lc1};
 Point(33) = {1, h/2 + 2*v -del, 0, lc1};
 Point(34) = {1 - a/2, h/2 + 2*v -del, 0, lc1};
 Point(35) = {1 - a, h/2 + v, 0, lc1};

 Line(31) = {31, 32};
 Line(32) = {32, 33};
 Line(33) = {33, 34};
 Line(34) = {34, 35};
 Line(35) = {35, 31};

 Line Loop(3) = {31, 32, 33, 34, 35};

// Hexagon 4

 Point(41) = {1 - a/2, 1 - h/2 -2*v +del, 0, lc1};
 Point(42) = {1, 1 - h/2 - 2*v + del, 0, lc1};
 Point(43) = {1, 1 - h/2 -del, 0, lc1};
 Point(44) = {1 - a/2, 1 - h/2 - del, 0, lc1};
 Point(45) = {1 - a,1 - h/2 - v, 0, lc1};

 Line(41) = {41, 42};
 Line(42) = {42, 43};
 Line(43) = {43, 44};
 Line(44) = {44, 45};
 Line(45) = {45, 41};

 Line Loop(4) = {41, 42, 43, 44, 45};

// Hexagon 5

 Point(51) = {0.5 - a/2, 1 - v + del, 0, lc1};
 Point(52) = {0.5 + a/2, 1 - v + del, 0, lc1};
 Point(53) = {0.5 + a, 1, 0, lc1};
 Point(54) = {0.5 - a, 1, 0, lc1};

 Line(51) = {51, 52};
 Line(52) = {52, 53};
 Line(53) = {53, 54};
 Line(54) = {54, 51};

 Line Loop(5) = {51, 52, 53, 54};

// Hexagon 6

 Point(61) = {0, 1 - h/2 - 2*v + del, 0, lc1};
 Point(62) = {a/2, 1 - h/2 - 2*v + del, 0, lc1};
 Point(63) = {a, 1 - h/2 - v, 0, lc1};
 Point(64) = {a/2, 1 - h/2 - del, 0, lc1};
 Point(65) = {0, 1 - h/2 - del, 0, lc1};

 Line(61) = {61, 62};
 Line(62) = {62, 63};
 Line(63) = {63, 64};
 Line(64) = {64, 65};
 Line(65) = {65, 61};

 Line Loop(6) = {61, 62, 63, 64, 65};

// Hexagon 7

 Point(71) = {0.5 - a/2, 0.5 - v + del, 0, lc1};
 Point(72) = {0.5 + a/2, 0.5 - v + del, 0, lc1};
 Point(73) = {0.5 + a, 0.5, 0, lc1};
 Point(74) = {0.5 + a/2, 0.5 + v - del, 0, lc1};
 Point(75) = {0.5 - a/2, 0.5 + v - del, 0, lc1};
 Point(76) = {0.5 - a, 0.5, 0, lc1};

 Line(71) = {71, 72};
 Line(72) = {72, 73};
 Line(73) = {73, 74};
 Line(74) = {74, 75};
 Line(75) = {75, 76};
 Line(76) = {76, 71};

 Line Loop(7) = {71, 72, 73, 74, 75, 76};
 Line Loop(11) = {76, 75, 74, 73, 72, 71};

 Plane Surface(1) = {1};
 Plane Surface(2) = {2};
 Plane Surface(3) = {3};
 Plane Surface(4) = {4};
 Plane Surface(5) = {5};
 Plane Surface(6) = {6};
 Plane Surface(7) = {7};


// Outer box

 Point(1) = {0, 0, 0, lc1};
 Point(2) = {1, 0, 0, lc1};
 Point(3) = {1, 1, 0, lc1};
 Point(4) = {0, 1, 0, lc1};

 Line(1) = {1, 21};
 Line(2) = {22, 2};
 Line(3) = {2, 32};
 Line(4) = {33, 42};
 Line(5) = {43, 3};
 Line(6) = {3, 53};
 Line(7) = {54, 4};
 Line(8) = {4, 65};
 Line(9) = {61, 15};
 Line(10) = {11, 1};


 Line Loop(10) = {1, -24, -23, -22, 2, 3, -31, -35, -34, -33, 4, -41, -45, -44, -43, 5, 6, -52, -51, -54, 7, 8, -64, -63, -62, -61, 9, -14, -13, -12, -11, 10};

 Plane Surface(10) = {10, 11};

 Physical Surface(1) = {1, 2, 3, 4, 5, 6, 7};
 Physical Surface(2) = {10};

 Color Blue{ Surface{ 1 }; }
 Color Red{ Surface{ 2 }; }
