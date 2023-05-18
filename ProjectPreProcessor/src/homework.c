#include "fem.h"

//
// Ici, vous pouvez définir votre géométrie :-)
//  (1) Raffiner intelligemment.... (yes )
//  (2) Construire la geometrie avec OpenCascade 
//  (3) Construire la geometrie avec les outils de GMSH
//  (4) Obtenir la geometrie en lisant un fichier .geo de GMSH


double geoSize(double x, double y) {

    femGeo* theGeometry = geoGetGeometry();
    return theGeometry->h * (1.0 - 0.5 * x);
}


void geoMeshGenerate() {

    femGeo* theGeometry = geoGetGeometry();
    geoSetSizeCallback(geoSize);   

    /*
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
    */

    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;
    
    int ierr;
    double r = w/4;
    double t = w/8;

    int rect[] = {2,gmshModelOccAddRectangle(-w / 2,-h / 4,0.0,w,h / 2,-1,0.0,&ierr) };
    int diskTop[] = {2,gmshModelOccAddDisk(0.0,h / 4.0,0.0,r,r,-1,NULL,0,NULL,0,&ierr) };
    int diskBottom[] = {2,gmshModelOccAddDisk(0.0,-h / 4.0,0.0,r,r,-1,NULL,0,NULL,0,&ierr) };
    int rectIn[] = {2,gmshModelOccAddRectangle(-w / 2 - t,-h / 4,0.0,w - t*2,h / 2,-1,0.0,&ierr) };
    int diskTopIn[] = { 2,gmshModelOccAddDisk(0.0,h / 4.0,0.0,r-t,r-t,-1,NULL,0,NULL,0,&ierr) };
    int diskBottomIn[] = { 2,gmshModelOccAddDisk(0.0,h / 4.0,0.0,r-t,r-t,-1,NULL,0,NULL,0,&ierr) };

    gmshModelOccCut(rect,2,disk,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr); 
    gmshModelOccCut(rect,2,slit,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr); 
    gmshModelOccSynchronize(&ierr); 

    if (theGeometry->elementType == FEM_QUAD) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshOptionSetNumber("Mesh.RecombineAll",1,&ierr);
        gmshOptionSetNumber("Mesh.Algorithm",8,&ierr);  
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm",1.0,&ierr); 
        gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  
        gmshModelMeshGenerate(2,&ierr);  }
  
    if (theGeometry->elementType == FEM_TRIANGLE) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshModelMeshGenerate(2,&ierr);  }
 
    return;
}


void geoMeshGenerateGeo() {

    femGeo* theGeometry = geoGetGeometry();
    geoSetSizeCallback(geoSize);   


    /*
    4 ------------------ 3
    |                    |
    |                    |
    5 ------- 6          |
               \         |
                )        |
               /         |
    8 ------- 7          |
    |                    |
    |                    |
    1 ------------------ 2
    */

    int ierr;
    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;
    double r = w/4;
    double lc = theGeometry->h;

    int p1 = gmshModelGeoAddPoint(-w/2, -h/2, 0., lc, 1, &ierr);
    int p2 = gmshModelGeoAddPoint( w/2, -h/2, 0., lc, 2, &ierr);
    int p3 = gmshModelGeoAddPoint( w/2,  h/2, 0., lc, 3, &ierr);
    int p4 = gmshModelGeoAddPoint(-w/2,  h/2, 0., lc, 4, &ierr);
    int p5 = gmshModelGeoAddPoint(-w/2,    r, 0., lc, 5, &ierr);
    int p6 = gmshModelGeoAddPoint(0.,      r, 0., lc, 6, &ierr);
    int p7 = gmshModelGeoAddPoint(0.,     -r, 0., lc, 7, &ierr);
    int p8 = gmshModelGeoAddPoint(-w/2,   -r, 0., lc, 8, &ierr);
    int p9 = gmshModelGeoAddPoint(0.,     0., 0., lc, 9, &ierr); // center of circle


    int l1 = gmshModelGeoAddLine(p1, p2, 1, &ierr);
    int l2 = gmshModelGeoAddLine(p2, p3, 2, &ierr);
    int l3 = gmshModelGeoAddLine(p3, p4, 3, &ierr);
    int l4 = gmshModelGeoAddLine(p4, p5, 4, &ierr);
    int l5 = gmshModelGeoAddLine(p5, p6, 5, &ierr);
    int l6 = gmshModelGeoAddCircleArc(p7, p9, p6, 6, 0., 0., 0., &ierr); // NB : the direction of the curve is reversed
    int l7 = gmshModelGeoAddLine(p7, p8, 7, &ierr);
    int l8 = gmshModelGeoAddLine(p8, p1, 8, &ierr);

    int lTags[] = {l1, l2, l3, l4, l5, -l6, l7, l8}; // NB : "-l6" because the curve is reversed 
    int c1[] = {1};
    c1[0] = gmshModelGeoAddCurveLoop(lTags, 8, 1, 0, &ierr);  
    int s1 = gmshModelGeoAddPlaneSurface(c1, 1, 1, &ierr);
    gmshModelGeoSynchronize(&ierr);


    if (theGeometry->elementType == FEM_QUAD) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshOptionSetNumber("Mesh.RecombineAll",1,&ierr);
        gmshOptionSetNumber("Mesh.Algorithm",8,&ierr);  
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm",1.0,&ierr); 
        gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  
        gmshModelMeshGenerate(2,&ierr);  }
  
    if (theGeometry->elementType == FEM_TRIANGLE) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshModelMeshGenerate(2,&ierr);  }

 //   gmshFltkRun(&ierr);
}


void geoMeshGenerateGeoFile(const char *filename){
    femGeo* theGeometry = geoGetGeometry();
    int ierr;
    gmshOpen(filename, &ierr); ErrorGmsh(ierr);
    if (theGeometry->elementType == FEM_QUAD) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshOptionSetNumber("Mesh.RecombineAll",1,&ierr);
        gmshOptionSetNumber("Mesh.Algorithm",8,&ierr); 
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm",1.0,&ierr); 
        gmshModelGeoMeshSetRecombine(2,1,45,&ierr); 
        gmshModelMeshGenerate(2,&ierr);  }
 
    if (theGeometry->elementType == FEM_TRIANGLE) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshModelMeshGenerate(2,&ierr);  }
    return;
}
