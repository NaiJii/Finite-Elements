w  = 1;
h  = 1;
r  = w/4;
lc = 0.1;

// coin inf�rieur gauche
Point(1) = {-w/2, -h/2, 0, lc}; 
// coin inf�rieur droit
Point(2) = { w/2, -h/2, 0, lc};
// coin sup�rieur droit
Point(3) = { w/2,  h/2, 0, lc};
// coin sup�rieur gauche
Point(4) = {-w/2,  h/2, 0, lc};
// coin sup�rieur gauche du trou
Point(5) = {-w/2,  r  , 0, lc};
// centre sup�rieur du trou
Point(6) = {   0,  r  , 0, lc};
// centre inf�rieur du trou
Point(7) = {   0, -r  , 0, lc};
// coin inf�rieur gauche du trou
Point(8) = {-w/2, -r  , 0, lc};
// centre 
Point(9) = {   0,  0  , 0, lc};

// ligne du bas
Line(1) = {1,2};
// ligne de droite
Line(2) = {2,3};
// ligne du haut
Line(3) = {3,4};
// ligne sup�rieure de gauche
Line(4) = {4,5};
// ligne sup�rieur du trou
Line(5) = {5,6};
// cercle centr� en {9} de diam�tre {7,6}
Circle(6) = {7,9,6};
// ligne inf�rieure du trou
Line(7) = {7,8};
// ligne inf�rieure de gauche 
Line(8) = {8,1};

Curve Loop(1) = {1:5, -6, 7:8};
// on prend chaque ligne 1 � 5, on inverse la ligne 6 car elle est invers�e, et on prend les lignes 7 et 8.
Plane Surface(1) = {1};
// on fait une surface plane avec la courbe 1.



