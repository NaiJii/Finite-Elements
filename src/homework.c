#include "fem.h"
#include <math.h>


// Il faut un fifrelin generaliser ce code.....
//  (1) Ajouter l'axisymétrique !    (mandatory)
//  (2) Ajouter les conditions de Neumann !   (mandatory)  
//  (3) Ajouter les conditions en normal et tangentiel !   (strongly advised)
//  (4) Et remplacer le solveur plein par un truc un fifrelin plus subtil  (mandatory)

// TODO:
// 1. ?
// 2. Voir les dernières slides
// 3. Voir les dernières slides
// 4. Ils parlent d'un solveur creux j'imagine, comme vu au cours.

double *femElasticitySolve(femProblem *theProblem) {
    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theMesh = theGeometry->theElements;
     
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,iEdge,i,j,d,map[4],mapX[4],mapY[4];

    // enlève les erreurs sur Windows
    memset(x, 0, sizeof(x));
    memset(y, 0, sizeof(y));
    memset(phi, 0, sizeof(phi));
    memset(dphidxsi, 0, sizeof(dphidxsi));
    memset(dphideta, 0, sizeof(dphideta));
    memset(dphidx, 0, sizeof(dphidx));
    memset(dphidy, 0, sizeof(dphidy));
    memset(map, 0, sizeof(map));
    memset(mapX, 0, sizeof(mapX));
    memset(mapY, 0, sizeof(mapY));
    
    int nLocal = theMesh->nLocalNode;

    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;      
    double rho = theProblem->rho;
    double g   = theProblem->g;
    double **A = theSystem->A;
    double *B  = theSystem->B;  
    
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; j++) {
            map[j]  = theMesh->elem[iElem*nLocal+j];
            mapX[j] = 2*map[j];
            mapY[j] = 2*map[j] + 1;
            x[j]    = theNodes->X[map[j]];
            y[j]    = theNodes->Y[map[j]];} 
        
        for (iInteg=0; iInteg < theRule->n; iInteg++) {    
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];  
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            
            double dxdxsi = 0.0;
            double dxdeta = 0.0;
            double dydxsi = 0.0; 
            double dydeta = 0.0;
            for (i = 0; i < theSpace->n; i++) {  
                dxdxsi += x[i]*dphidxsi[i];       
                dxdeta += x[i]*dphideta[i];   
                dydxsi += y[i]*dphidxsi[i];   
                dydeta += y[i]*dphideta[i]; }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            
            for (i = 0; i < theSpace->n; i++) {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; }            
            for (i = 0; i < theSpace->n; i++) { 
                for(j = 0; j < theSpace->n; j++) {
                    A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + 
                                            dphidy[i] * c * dphidy[j]) * jac * weight;                                                                                            
                    A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + 
                                            dphidy[i] * c * dphidx[j]) * jac * weight;                                                                                           
                    A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + 
                                            dphidx[i] * c * dphidy[j]) * jac * weight;                                                                                            
                    A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + 
                                            dphidx[i] * c * dphidx[j]) * jac * weight; }}
             for (i = 0; i < theSpace->n; i++) {
                B[mapY[i]] -= phi[i] * g * rho * jac * weight; 
                // Calculer les dérivées normales et tangentes
                double nx = -dphidx[i];
                double ny = -dphidy[i];
                double tx = dphidy[i];
                double ty = -dphidx[i];

                //Conditions en normal et tangentiel
                int nodeIndex = map[i];
                femBoundaryCondition** conditions = theProblem->conditions;
                int nConditions = theProblem->nBoundaryConditions;

                for (j = 0; j < nConditions; j++) {
                    femBoundaryCondition* condition = conditions[j];
                    if (condition->type == NEUMANN_N) {
                        // Condition en normal
                        double normalValue = condition->value;
                        B[mapX[i]] += nx * normalValue * jac * weight;
                        B[mapY[i]] += ny * normalValue * jac * weight;
                    }
                    else if (condition->type == NEUMANN_T) {
                        // Condition en tangentiel
                        double tangentValue = condition->value;
                        B[mapX[i]] += tx * tangentValue * jac * weight;
                        B[mapY[i]] += ty * tangentValue * jac * weight;
                    }
                }
             }

             //Conditions de Neumann
             int nEdges = theGeometry->theEdges->nElem;
             int nLocalEdge = theGeometry->theEdges->nLocalNode;
             int nBoundaryConditions = theProblem->nBoundaryConditions;
             femBoundaryCondition** conditions = theProblem->conditions;
             femNodes* edgesNodes = theGeometry->theEdges->nodes;

             for (iEdge = 0; iEdge < nEdges; iEdge++) {
                 for (i = 0; i < nLocalEdge; i++) {
                     map[i] = theGeometry->theEdges->elem[iEdge * nLocalEdge + i];
                     mapX[i] = 2 * map[i];
                     mapY[i] = 2 * map[i] + 1;
                     x[i] = edgesNodes->X[map[i]];
                     y[i] = edgesNodes->Y[map[i]];
                 }

                 double edgeLength = sqrt(pow(x[1] - x[0], 2) + pow(y[1] - y[0], 2));
                 double* edgeNormal = malloc(2 * sizeof(double));
                 edgeNormal[0] = (y[1] - y[0]) / edgeLength;
                 edgeNormal[1] = -(x[1] - x[0]) / edgeLength;

                 for (iInteg = 0; iInteg < theRule->n; iInteg++) {
                     double xsi = theRule->xsi[iInteg];
                     double eta = theRule->eta[iInteg];
                     double weight = theRule->weight[iInteg];
                     femDiscretePhi2(theSpace, xsi, eta, phi);
                     femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);

                     double jac = edgeLength / 2.0;
                     double normalJac = jac * weight;

                     for (i = 0; i < theSpace->n; i++) {
                         double neumannValue = 0.0;
                         for (j = 0; j < nBoundaryConditions; j++) {
                             femBoundaryCondition* condition = conditions[j];
                             if (condition->type == NEUMANN_N) {
                                 neumannValue += condition->value * phi[i];
                             }
                             else if (condition->type == NEUMANN_T) {
                                 double tangentComponent = condition->value * edgeNormal[0];
                                 neumannValue += tangentComponent * dphidxsi[i] +
                                     tangentComponent * dphideta[i];
                             }
                         }
                         B[mapX[i]] += neumannValue * normalJac * edgeNormal[0];
                         B[mapY[i]] += neumannValue * normalJac * edgeNormal[1];
                     }
                 }

                 free(edgeNormal);
             }

        }
    } 
  
    int *theConstrainedNodes = theProblem->constrainedNodes;     
    for (int i=0; i < theSystem->size; i++) {
        if (theConstrainedNodes[i] != -1) {
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femFullSystemConstrain(theSystem,i,value); }}
                            
    return femFullSystemEliminate(theSystem);
}



//BROUILLION:
// Calcul du vecteur des contraintes Neumann
//double neumannX = 0.0;
///double neumannY = 0.0;
//for (i = 0; i < theSpace->n; i++) {
    //neumannX += phi[i];
    //neumannY += phi[i];// *theProblem->conditions[map[i]]->value* dphidy[i]}

// Ajout des contributions des conditions de Neumann dans le système
//for (i = 0; i < theSpace->n; i++) {
    //B[mapX[i]] += neumannX * jac * weight;
    //B[mapY[i]] += neumannY * jac * weight;}