w  = 0.306 * 3;
h  = 1.47 * 0.5 * 3;
r  = w/4;
lc = 0.1;


// on définit les points
Point(1) = {0, h/2, 0, lc}; 
Point(2) = {0, h/4, 0, lc}; 
Point(3) = {0, -h/4, 0, lc}; 
Point(4) = {0, -h/2, 0, lc}; 
Point(5) = {w/2, -h/4, 0, lc}; 
Point(6) = {w/2, h/4, 0, lc}; 

// on définit les lignes
Line(1) = {1,4}; 
Line(2) = {5,6}; 

//Circle(3) = {4,3,5}; 
//Circle(4) = {6,2,1}; 
Line(3) = {4,5};
Line(4) = {6,1};

// on définit les courbes
Curve Loop(1) = {1, 3, 2, 4};
Plane Surface(1) = {1};



