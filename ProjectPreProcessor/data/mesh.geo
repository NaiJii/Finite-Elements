w  = 1;
h  = 1;
r  = w/4;
lc = 0.1;

// coin inférieur gauche
Point(1) = {-w/2, -h/2, 0, lc}; 
// coin inférieur droit
Point(2) = { w/2, -h/2, 0, lc};
// coin supérieur droit
Point(3) = { w/2,  h/2, 0, lc};
// coin supérieur gauche
Point(4) = {-w/2,  h/2, 0, lc};
// coin supérieur gauche du trou
Point(5) = {-w/2,  r  , 0, lc};
// centre supérieur du trou
Point(6) = {   0,  r  , 0, lc};
// centre inférieur du trou
Point(7) = {   0, -r  , 0, lc};
// coin inférieur gauche du trou
Point(8) = {-w/2, -r  , 0, lc};
// centre 
Point(9) = {   0,  0  , 0, lc};

// ligne du bas
Line(1) = {1,2};
// ligne de droite
Line(2) = {2,3};
// ligne du haut
Line(3) = {3,4};
// ligne supérieure de gauche
Line(4) = {4,5};
// ligne supérieur du trou
Line(5) = {5,6};
// cercle centré en {9} de diamètre {7,6}
Circle(6) = {7,9,6};
// ligne inférieure du trou
Line(7) = {7,8};
// ligne inférieure de gauche 
Line(8) = {8,1};

Curve Loop(1) = {1:5, -6, 7:8};
// on prend chaque ligne 1 à 5, on inverse la ligne 6 car elle est inversée, et on prend les lignes 7 et 8.
Plane Surface(1) = {1};
// on fait une surface plane avec la courbe 1.



