w  = 0.306 * 3;
h  = 1.47 * 0.5 * 3;
r  = w/4;
lc = 0.1;


// on définit les points
Point(1) = {-w/2, -h/4, 0, lc}; // g, b
Point(2) = {0, -h/4, 0, lc}; // b
Point(3) = {w/2, -h/4, 0, lc}; // d, b
Point(4) = {w/2, h/4, 0, lc}; // d, h
Point(5) = {0, h/4, 0, lc}; // h
Point(6) = {-w/2, h/4, 0, lc}; // g, h

Point(7) = {-w/2+r, -h/4, 0, lc}; // g, b
Point(8) = {w/2-r, -h/4, 0, lc}; // d, b
Point(9) = {w/2-r, h/4, 0, lc}; // d, h
Point(10) = {-w/2+r, h/4, 0, lc}; // g, h

// on définit les lignes
Line(1) = {6,1}; // g
Line(2) = {1,3}; // d
Line(3) = {3,4}; // h 
Line(4) = {4,6}; // b

Line(7) = {10,7}; // g
Line(8) = {8,9}; // d


// cercle centré en {5} de diamètre {4,3}
Circle(5) = {1,2,3}; // b
Circle(6) = {4,5,6}; // h

Circle(9) = {7,2,8}; // b
Circle(10) = {9,5,10}; // h


// on définit les courbes
Curve Loop(1) = {1, 5, 3, 6};
Curve Inner(2) = {7, 9, 8, 10};
Plane Surface(1) = {1, -2};



