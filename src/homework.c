#include "fem.h"
#include <math.h>


// Il faut un fifrelin generaliser ce code.....
//  (1) Ajouter l'axisymétrique !    (mandatory)
//  (2) Ajouter les conditions de Neumann !   (mandatory)  
//  (3) Ajouter les conditions en normal et tangentiel !   (strongly advised)
//  (4) Et remplacer le solveur plein par un truc un fifrelin plus subtil  (mandatory)

int cmp(void* v, const void* a, const void* b) {
    double* vec = (double*)v; 
	const double x = vec[*(const int*)a];
    const double y = vec[*(const int*)b];
    return y - x; // TODO: check if this is correct
}

// renumeration des noeuds
void femMeshRenum(femMesh* theMesh, femRenumType renumType) {  
    if (renumType == RENUM_NONE)
        return;
    
    const int num = theMesh->nodes->nNodes;
    int* renumNodes = malloc(num * sizeof(int));

    if (!renumNodes) {
        Error("Memory allocation failed in femMeshRenum");
    }

    for (int i = 0; i < num; ++i)
		renumNodes[i] = i;

    double* v = NULL;

    if (renumType == RENUM_X) {
        v = theMesh->nodes->X;
    }
    else if (renumType == RENUM_Y) {
        v = theMesh->nodes->Y;
	}
    else {
        Error("No such renumType in femMeshRenum");
	}

    qsort_s(renumNodes, num, sizeof(int), cmp, v);

    for (int i = 0; i < num; ++i) {
		const int idx = renumNodes[i];
        theMesh->_enum[idx] = i;
    }

	free(renumNodes);
}

double* femElasticitySolve(femProblem* theProblem) {
    femFullSystem* theSystem = theProblem->system;
    femIntegration* theRule = theProblem->rule;
    femDiscrete* theSpace = theProblem->space;
    femGeo* theGeometry = theProblem->geometry;
    femNodes* theNodes = theGeometry->theNodes;
    femMesh* theMesh = theGeometry->theElements;

    int nLocal = theMesh->nLocalNode;
    int problemCase = theProblem->planarStrainStress;

    double* x = malloc(sizeof(double) * nLocal);
    double* y = malloc(sizeof(double) * nLocal);
    double* phi = malloc(sizeof(double) * nLocal);
    double* dphidxsi = malloc(sizeof(double) * nLocal);
    double* dphideta = malloc(sizeof(double) * nLocal);
    double* dphidx = malloc(sizeof(double) * nLocal);
    double* dphidy = malloc(sizeof(double) * nLocal);
    if (!x || !y || !phi || !dphidxsi || !dphideta || !dphidx || !dphidy)
        Error("Allocation error");

    int iElem, iInteg, iEdge, i, j, d;
    int* map = malloc(sizeof(int) * nLocal);
    int* mapX = malloc(sizeof(int) * nLocal);
    int* mapY = malloc(sizeof(int) * nLocal);
    if (!map || !mapX || !mapY)
        Error("Allocation error");

    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;      
    double rho = theProblem->rho;
    double g   = theProblem->g;
    double **A = theSystem->A;
    double *B  = theSystem->B;  
    femRenumType renumType = RENUM_NONE;

    femMeshRenum(theMesh, renumType);

    // band solver ... //    
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j = 0; j < nLocal; j++) {
            map[j] = theMesh->elem[iElem * nLocal + j];
            mapX[j] = 2 * map[j];
            mapY[j] = 2 * map[j] + 1;
            x[j] = theNodes->X[map[j]];
            y[j] = theNodes->Y[map[j]];
        }

        for (iInteg = 0; iInteg < theRule->n; iInteg++) {
            double xsi = theRule->xsi[iInteg];
            double eta = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];
            femDiscretePhi2(theSpace, xsi, eta, phi);
            femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);

            double dxdxsi = 0.0;
            double dxdeta = 0.0;
            double dydxsi = 0.0;
            double dydeta = 0.0;
            // axisymmetry
            double xphi = 0.0;
            double yphi = 0.0; // unused because of gravity

            for (i = 0; i < theSpace->n; i++) {
                dxdxsi += x[i] * dphidxsi[i];
                dxdeta += x[i] * dphideta[i];
                dydxsi += y[i] * dphidxsi[i];
                dydeta += y[i] * dphideta[i];
                xphi += x[i] * phi[i];
                yphi += y[i] * phi[i]; // unused because of gravity
            }

            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);

            for (i = 0; i < theSpace->n; i++) {
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
            }

            for (i = 0; i < theSpace->n; i++) {
                for (j = 0; j < theSpace->n; j++) {

                    if (problemCase == AXISYM) {
                        A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] * xphi +
                            dphidy[i] * c * dphidy[j] * xphi +
                            phi[i] * (b * dphidx[j] + a * phi[j] / xphi) +
                            dphidx[i] * b * phi[j]) * jac * weight;
                        A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] * xphi +
                            dphidy[i] * c * dphidx[j] * xphi +
                            phi[i] * b * dphidy[j]) * jac * weight;
                        A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] * xphi +
                            dphidx[i] * c * dphidy[j] * xphi +
                            dphidy[i] * b * phi[j]) * jac * weight;
                        A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] * xphi +
                            dphidx[i] * c * dphidx[j] * xphi) * jac * weight;
                    }
                    else {
                        A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] +
                            dphidy[i] * c * dphidy[j]) * jac * weight;
                        A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] +
                            dphidy[i] * c * dphidx[j]) * jac * weight;
                        A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] +
                            dphidx[i] * c * dphidy[j]) * jac * weight;
                        A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] +
                            dphidx[i] * c * dphidx[j]) * jac * weight;
                    }
                }
            }


            for (i = 0; i < theSpace->n; i++) {
                if (problemCase == AXISYM)
                    B[mapY[i]] -= phi[i] * g * rho * jac * weight * xphi;
                else
                    B[mapY[i]] -= phi[i] * g * rho * jac * weight;
            }
        }
    }

    int nConditions = theProblem->nBoundaryConditions;
    for (int c = 0; c < nConditions; c++) {
        femBoundaryCondition* cond = theProblem->conditions[i];
        if (!cond)
            Error("No condition");

        femBoundaryType type = cond->type;
        femMesh* mesh = cond->domain->mesh;
        if (!mesh)
            Error("No mesh");

        int* elem = cond->domain->elem;
        int nElem = cond->domain->nElem;
        double* X = mesh->nodes->X;
        double* Y = mesh->nodes->X;


    }

    int *theConstrainedNodes = theProblem->constrainedNodes;     
    for (int i=0; i < theSystem->size; i++) {
        if (theConstrainedNodes[i] != -1) {
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femFullSystemConstrain(theSystem,i,value); }}
                            
    return femFullSystemEliminate(theSystem);
}
