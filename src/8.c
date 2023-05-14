#include "fem.h"

/*
Une pression c'est une condition de Neuman normale.

*/


void geoMeshGenerate() {
    int ierr;
    femGeo* theGeometry = geoGetGeometry();

    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;
    double r = w * 0.5;

    double rectPos[] = { 0.0, h * 0.5 };
    double rectSize[] = { w, h };
    double topDiskPos[] = { w * 0.5, h * 1.5 };
    double bottomDiskPos[] = { w * 0.5, h * 0.5 };

    int rect[] = { 2, gmshModelOccAddRectangle(rectPos[0], rectPos[1], 0.0, rectSize[0], rectSize[1], -1, 0.0, &ierr)};
    int halfRect[] = { 2, gmshModelOccAddRectangle(rectPos[0], 0.0, 0.0, rectSize[0] * 0.5, rectSize[1] * 2.0, -1, 0, &ierr) };
    int diskTop[] = { 2, gmshModelOccAddDisk(topDiskPos[0], topDiskPos[1], 0.0, r, r, -1, NULL, 0, NULL, 0, &ierr)};
    int diskBottom[] = { 2, gmshModelOccAddDisk(bottomDiskPos[0], bottomDiskPos[1], 0.0, r, r, -1, NULL, 0, NULL, 0, &ierr) };

    gmshModelOccFuse(rect, rect[0], diskTop, diskTop[0], NULL,NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(rect, rect[0], diskBottom, diskBottom[0], NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccCut(rect, rect[0], halfRect, halfRect[0], NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);

    gmshModelOccSynchronize(&ierr);

    if (theGeometry->elementType == FEM_QUAD) {
        gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
        gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
        gmshOptionSetNumber("Mesh.Algorithm", 8, &ierr);
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);
        gmshModelGeoMeshSetRecombine(2, 1, 45, &ierr);
        gmshModelMeshGenerate(2, &ierr);
    }

    if (theGeometry->elementType == FEM_TRIANGLE) {
        gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
        gmshModelMeshGenerate(2, &ierr);
    }

    return;
}


double* femElasticitySolve(femProblem* theProblem)
{

    femFullSystem* theSystem = theProblem->system;
    femIntegration* theRule = theProblem->rule;
    femDiscrete* theSpace = theProblem->space;
    femGeo* theGeometry = theProblem->geometry;
    femNodes* theNodes = theGeometry->theNodes;
    femMesh* theMesh = theGeometry->theElements;


    double x[4], y[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4];
    int iElem, iInteg, iEdge, i, j, d, map[4], mapX[4], mapY[4];

    int nLocal = theMesh->nLocalNode;

    double a = theProblem->A;
    double b = theProblem->B;
    double c = theProblem->C;
    double rho = theProblem->rho;
    double g = theProblem->g;
    double** A = theSystem->A;
    double* B = theSystem->B;


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
            for (i = 0; i < theSpace->n; i++) {
                dxdxsi += x[i] * dphidxsi[i];
                dxdeta += x[i] * dphideta[i];
                dydxsi += y[i] * dphidxsi[i];
                dydeta += y[i] * dphideta[i];
            }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);

            for (i = 0; i < theSpace->n; i++) {
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
            }
            for (i = 0; i < theSpace->n; i++) {
                for (j = 0; j < theSpace->n; j++) {
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
            for (i = 0; i < theSpace->n; i++) {
                B[mapY[i]] -= phi[i] * g * rho * jac * weight;
            }
        }
    }

    int* theConstrainedNodes = theProblem->constrainedNodes;
    for (int i = 0; i < theSystem->size; i++) {
        if (theConstrainedNodes[i] != -1) {
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femFullSystemConstrain(theSystem, i, value);
        }
    }

    return femFullSystemEliminate(theSystem);
}