#include "fem.h"

void geoMeshGenerate() {
	femGeo* theGeometry = geoGetGeometry();

	double w = theGeometry->LxPlate;
	double h = theGeometry->LyPlate;

	int ierr;
	double r = w / 4;
	int idRect = gmshModelOccAddRectangle(0.0, 0.0, 0.0, w, h, -1, 0.0, &ierr);
	int idDisk = gmshModelOccAddDisk(w / 2.0, h / 2.0, 0.0, r, r, -1, NULL, 0, NULL, 0, &ierr);
	int idSlit = gmshModelOccAddRectangle(w / 2.0, h / 2.0 - r, 0.0, w, 2.0 * r, -1, 0.0, &ierr);
	int rect[] = { 2,idRect };
	int disk[] = { 2,idDisk };
	int slit[] = { 2,idSlit };

	gmshModelOccCut(rect, 2, disk, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
	gmshModelOccCut(rect, 2, slit, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
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

    int nLocal = theMesh->nLocalNode;

    double x[4], y[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4];
    int map[4];

    double a = theProblem->A;
    double b = theProblem->B;
    double c = theProblem->C;
    double rho = theProblem->rho;
    double g = theProblem->g;
    double** A = theSystem->A;
    double* B = theSystem->B;

	for (int Elem = 0; Elem < theMesh->nElem; Elem++){
		for (int i = 0; i < nLocal; i++) {
			map[i] = theMesh->elem[Elem * nLocal + i];
			x[i] = theMesh->nodes->X[map[i]];
			y[i] = theMesh->nodes->Y[map[i]];
		}
		for (int j = 0; j < theRule->n; j++) {
			double eta = theRule->eta[j];
			double xsi = theRule->xsi[j];
			double weight = theRule->weight[j];

			theSpace->phi2(xsi, eta, phi);
			theSpace->dphi2dx(xsi, eta, dphidxsi, dphideta);

			double dxdxsi = 0;
			double dydxsi = 0;
			double dxdeta = 0;
			double dydeta = 0;

			for (int i = 0; i < theSpace->n; i++) {
				dxdxsi += x[i] * dphidxsi[i];
				dydxsi += y[i] * dphidxsi[i];
				dxdeta += x[i] * dphideta[i];
				dydeta += y[i] * dphideta[i];
			}
			double J =  dydeta* dxdxsi - dydxsi*dxdeta ;
			for (int i = 0; i < theSpace->n; i++) {
				dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / J;
				dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / J;
			}

			for (int i = 0; i < theSpace->n; i++) {
				double JW = J * weight;
				theSystem->B[2 * map[i] + 1] -= phi[i] * rho * g * g * J * weight;
				theSystem->B[2 * map[i]] = 0;
				for (int j = 0; j < theSpace->n; j++) {
					theSystem->A[2 * map[i]][2 * map[j]] += (dphidx[i] * dphidx[j] * a + dphidy[i] * dphidy[j] * c) * JW;
					theSystem->A[2 * map[i]][2 * map[j]+1] += (dphidx[i] * dphidy[j] * b + dphidy[i] * dphidx[j] * c) * JW;
					theSystem->A[2 * map[i] + 1][2 * map[j]] += (dphidy[i] * dphidx[j] * b + dphidx[i] * dphidy[j] * c) * JW;
					theSystem->A[2 * map[i]+1][2 * map[j]+1] += (dphidy[i] * dphidy[j] * a + dphidx[i] * dphidx[j] * c) * JW;

				}
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
