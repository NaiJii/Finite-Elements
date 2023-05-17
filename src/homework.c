#include "fem.h"
#include <math.h>

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

// Il faut un fifrelin generaliser ce code.....
//  (1) Ajouter l'axisymétrique !    (mandatory)
//  (2) Ajouter les conditions de Neumann !   (mandatory)  
//  (3) Ajouter les conditions en normal et tangentiel !   (strongly advised)
//  (4) Et remplacer le solveur plein par un truc un fifrelin plus subtil  (mandatory)

double* vec = NULL;
int cmp(const void* a, const void* b) {
    //double* vec = (double*)v; 
	const double x = vec[*(const int*)a];
    const double y = vec[*(const int*)b];
    return y - x; // TODO: check if this is correct
}

// renumeration des noeuds
void femMeshRenumber(femMesh* theMesh, femRenumType renumType) {     
    const int num = theMesh->nodes->nNodes;
    int* renumNodes = malloc(num * sizeof(int));

    if (!renumNodes) {
        Error("Memory allocation failed in femMeshRenum");
    }

    for (int i = 0; i < num; ++i)
		renumNodes[i] = i;

    if (renumType == RENUM_X) {
        vec = theMesh->nodes->X;
    }
    else if (renumType == RENUM_Y) {
        vec = theMesh->nodes->Y;
	}
    else if (renumType != RENUM_NONE) {
        Error("No such renumType in femMeshRenum");
    }

    qsort(renumNodes, num, sizeof(int), cmp);

    for (int i = 0; i < num; ++i) {
        theMesh->number[renumNodes[i]] = i;
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

    femMeshRenumber(theMesh, FEM_XNUM);
    int bandWidth = femComputeBand(theMesh);

    printf("band = %d\n, size = %d\n", bandWidth, theSystem->size);

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
                            dphidy[i] * c * dphidy[j] * xphi + phi[i] * (b * dphidx[j] + a * phi[j] / xphi) +
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
        femBoundaryCondition* cond = theProblem->conditions[c];
        if (!cond)
            Error("No condition");

        printf("Domain: %s", cond->domain->name);

        femBoundaryType type = cond->type;
        femMesh* mesh = cond->domain->mesh;
        if (!mesh)
            Error("No mesh");

        int* elem = cond->domain->elem;
        int nElem = cond->domain->nElem;
        int nNodes = nElem + 1;
        double* X = mesh->nodes->X;
        double* Y = mesh->nodes->X;
        double* N = cond->domain->normals;
        double* T = cond->domain->tangents;

        if (cond->domain->nComponents == 0) {
            cond->domain->nComponents = mesh->nElem + 1;

            for (int e = 0; e < nElem; e++) {
                // both nodes of segment n
                int a = mesh->elem[elem[e] * 2];
                int b = mesh->elem[elem[e] * 2 + 1];
                double aPos[2] = { X[a], Y[a] };
                double bPos[2] = { X[b], Y[b] };
                // node n
                N[e * 2] += aPos[1] - bPos[1];
                N[e * 2 + 1] += bPos[0] - aPos[0];
                T[e * 2] += bPos[0] - aPos[0];
                T[e * 2 + 1] += bPos[1] - aPos[1];
                // node n+1
                N[e * 2 + 2] += aPos[1] - bPos[1];
                N[e * 2 + 3] += bPos[0] - aPos[0];
                T[e * 2 + 2] += bPos[0] - aPos[0];
                T[e * 2 + 3] += bPos[1] - aPos[1];
            }

            // toutes les tangentes et normales sont definies, on les normalise
            for (int n = 0; n < nNodes; n++) {

                double normN = sqrt(pow(N[n * 2], 2.0) + pow(N[n * 2 + 1], 2.0));
                if (normN == 0.0)
                    Error("0,0 position for T");

                N[n * 2] /= normN;
                N[n * 2 + 1] /= normN;

                double normT = sqrt(pow(T[n * 2], 2.0) + pow(T[n * 2 + 1], 2.0));
                if (normT == 0.0)
                    Error("0,0 position for T");

                T[n * 2] /= normT;
                T[n * 2 + 1] /= normT;

                // print the normal and tangent components
                //printf("T: %f %f \n N: %f %f \n", T[n * 2], T[n * 2 + 1], N[n * 2], N[n * 2 + 1]);
            }

            for (int k = 0; k < nElem + 1; k++) {
                printf("normale du noeud %d = (%f ; %f)\n", k, N[2 * k], N[2 * k + 1]);
                printf("tangente du noeud %d = (%f ; %f)\n", k, T[2 * k], T[2 * k + 1]);
            }
        }


        //if (type == DIRICHLET_N || type == DIRICHLET_T ||
        //    type == NEUMANN_N || type == NEUMANN_T) {
        //}

#if 0
        if (cnd->type == DIRICHLET_N || cnd->type == DIRICHLET_T) {
            // Itération sur tous les noeuds du domaine
            for (int j = 0; j < nElem + 1; j++) {
                // technique pour récupérer tous les noeuds sur base des éléments
                if (j == nElem) {
                    node0 = bndMesh->elem[2 * bndElem[j - 1] + 1];
                }
                else { node0 = bndMesh->elem[2 * bndElem[j]]; }

                // Combinaison linéaire des lignes et colonnes de la matrice pour changer (U, V) en (N, T) (voir rapport)
                double A_U, A_V, B_U, B_V, nx, ny, tx, ty;
                nx = normales[2 * j];
                ny = normales[2 * j + 1];
                tx = tangentes[2 * j];
                ty = tangentes[2 * j + 1];

                B_U = B[2 * node0];
                B_V = B[2 * node0 + 1];
                B[2 * node0] = nx * B_U + ny * B_V;
                B[2 * node0 + 1] = tx * B_U + ty * B_V;
                // Modification des lignes
                for (int k = 0; k < theSystem->size; k++) {
                    A_U = A[2 * node0][k];
                    A_V = A[2 * node0 + 1][k];
                    A[2 * node0][k] = nx * A_U + ny * A_V;
                    A[2 * node0 + 1][k] = tx * A_U + ty * A_V;
                }
                // Modification des colonnes
                for (int k = 0; k < theSystem->size; k++) {
                    A_U = A[k][2 * node0];
                    A_V = A[k][2 * node0 + 1];
                    A[k][2 * node0] = nx * A_U + ny * A_V;
                    A[k][2 * node0 + 1] = tx * A_U + ty * A_V;
                }
            }
        }
#endif

        if (type == DIRICHLET_N || type == DIRICHLET_T) {
            for (int n = 0; n < nNodes; n++) {
                int node = mesh->elem[elem[MIN(n, nElem - 1)]];
                int shift = type == DIRICHLET_N ? 0 : 1;
                femFullSystemConstrain(theSystem, node * 2 + shift, cond->value);
            }
        }

        if (type >= NEUMANN_X && type <= NEUMANN_T) {
            for (int e = 0; e < nElem; e++) {
                int a = mesh->elem[2 * elem[e]];
                int b = mesh->elem[2 * elem[e] + 1];
                double jac = sqrt((X[a] - X[b]) * (X[a] - X[b]) + (Y[a] - Y[b]) * (Y[a] - Y[b])) * 0.5;

                if (type == NEUMANN_X || type == NEUMANN_Y) {
                    int shift = type == NEUMANN_X ? 0 : 1;
                    B[a * 2 + shift] += jac * cond->value;
                    B[b * 2 + shift] += jac * cond->value;
                }
                else if (type == NEUMANN_N || type == NEUMANN_T) { 
                    double* xy = type == NEUMANN_N ? N : T;
                    B[a * 2] += jac * cond->value * xy[e * 2];
                    B[a * 2 + 1] += jac * cond->value * xy[2 * j + 1];
                    B[b * 2] += jac * cond->value * xy[2 * j + 2];
                    B[b * 2 + 1] += jac * cond->value * xy[2 * j + 3];
                }
            }
        }
    }

    int *theConstrainedNodes = theProblem->constrainedNodes;     
    for (int i=0; i < theSystem->size; i++) {
        if (theConstrainedNodes[i] != -1) {
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femBoundaryType type = theProblem->conditions[theConstrainedNodes[i]]->type;

            if (type == DIRICHLET_X || type == DIRICHLET_Y) {
                femFullSystemConstrain(theSystem, i, value);
            }
        }
    }
                            
    double* B = femFullSystemEliminate(theSystem);

    free(x);
    free(y);
    free(phi);
    free(dphidxsi);
    free(dphideta);
    free(dphidx);
    free(dphidy);    

    return B;
}

double* femTension(femProblem* theProblem, double* UV) {
    femDiscrete* theSpace = theProblem->space;
    femGeo* theGeometry = theProblem->geometry;
    femNodes* theNodes = theGeometry->theNodes;
    femMesh* theMesh = theGeometry->theElements;
    int nLocal = theMesh->nLocalNode;
    double a = theProblem->A;
    double b = theProblem->B;
    double c = theProblem->C;
    
    int taille = 4 + (theProblem->planarStrainStress != AXISYM ? 0 : 1);

    double* U = &UV[0];
    double* V = &UV[1];
    double* xsi = malloc(sizeof(double) * nLocal);
    double* eta = malloc(sizeof(double) * nLocal);
    double* dphidxsi = malloc(sizeof(double) * nLocal);
    double* dphideta = malloc(sizeof(double) * nLocal);
    double* dphidx = malloc(sizeof(double) * nLocal);
    double* dphidy = malloc(sizeof(double) * nLocal);
    femDiscreteXsi2(theSpace, xsi, eta);

    int* NEXT = malloc(sizeof(int) * 10);  
    double* e = malloc(sizeof(double) * taille);
    double* sigma = malloc(sizeof(double) * theNodes->nNodes * taille);

    for (int i = 0; i < theNodes->nNodes; i++) {
        int count = -1;
        int c;
        for (int j = 0; j < theMesh->nElem * nLocal; j++) { 
            if (theMesh->elem[j] == i) {
                c = j % nLocal;
                NEXT[count++] = j - c;     
            }
        }
        double* dphidx = malloc(sizeof(double) * nLocal);
        double* dphidy = malloc(sizeof(double) * nLocal);
        double* dphidxsi = malloc(sizeof(double) * nLocal);
        double* dphideta = malloc(sizeof(double) * nLocal);
        for (int j = 0; j < count; j++) {  
            double dxdxsi = 0.0;
            double dydxsi = 0.0;
            double dydeta = 0.0;
            double dxdeta = 0.0;
            for (int k = 0; k < nLocal; k++) {
                femDiscreteDphi2(theSpace, xsi[c], eta[c], dphidxsi, dphideta);
                dxdxsi += theNodes->X[theMesh->elem[NEXT[j] + k]] * dphidxsi[k];
                dydxsi += theNodes->Y[theMesh->elem[NEXT[j] + k]] * dphidxsi[k];
                dydeta += theNodes->Y[theMesh->elem[NEXT[j] + k]] * dphideta[k];
                dxdeta += theNodes->X[theMesh->elem[NEXT[j] + k]] * dphideta[k];
            }
            double A = dydeta * dxdxsi;
            double B = dydxsi * dxdeta;
            double J = fabs(A - B);
            for (int l = 0; l < nLocal; l++) {
                dphidx[l] = (dphidxsi[l] * dydeta - dphideta[l] * dydxsi) / J;
                dphidy[l] = (dphideta[l] * dxdxsi - dphidxsi[l] * dxdeta) / J;
            }
            for (int m = 0; m < nLocal; m++) {
                int index = NEXT[j] + m;
                e[0] += (dphidx[m] * U[2 * theMesh->elem[index]]);
                e[1] += (dphidy[m] * V[2 * theMesh->elem[index]]);
                e[2] += (dphidx[m] * V[2 * theMesh->elem[index]] + dphidy[m] * U[2 * theMesh->elem[index]]) / 2;
                e[3] += (dphidy[m] * U[2 * theMesh->elem[index]] + dphidx[m] * V[2 * theMesh->elem[index]]) / 2;
            }
        }
        int sigmaIndex = i * taille;
        sigma[sigmaIndex] = a * e[0] + b * e[1];
        sigma[sigmaIndex + 1] = a * e[1] + b * e[0];
        sigma[sigmaIndex + 2] = 2 * c * e[2];
        sigma[sigmaIndex + 3] = 2 * c * e[3];
        if (theProblem->planarStrainStress == AXISYM) {
            e[4] += U[2 * theMesh->elem[i]] / theNodes->X[i];
            sigma[sigmaIndex] += b * e[4];
            sigma[sigmaIndex + 1] += b * e[4];
            sigma[sigmaIndex + 4] = b * (e[0] + e[1]) + a * e[4];}
        free(dphidx);
        free(dphidy);
        free(dphidxsi);
        free(dphideta);
        free(e);
        free(NEXT);
    }
    return sigma;
}

int femComputeBand(femMesh* theMesh) {
    int map[4] = {0, 0, 0, 0};
    int nLocal = theMesh->nLocalNode;
    int myBand = 0;

    for (int i = 0; i < theMesh->nElem; i++) {
        printf("number[%d] = %d\n", i, theMesh->number[i]);
        for (int j = 0; j < nLocal; ++j) 
            map[j] = theMesh->number[theMesh->elem[i * nLocal + j]];
       
        int myMin = map[0];
        int myMax = map[0];
        for (int j = 1; j < nLocal; j++) {
            myMax = fmax(map[j], myMax);
            myMin = fmin(map[j], myMin);
        }
        if (myBand < (myMax - myMin)) 
            myBand = myMax - myMin;
    }

    return myBand + 1; // largeur de la bande
}
