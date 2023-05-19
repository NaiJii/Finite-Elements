#include "fem.h"
#include <math.h>

// Il faut un fifrelin generaliser ce code.....
//  (1) Ajouter l'axisymétrique !    (mandatory)
//  (2) Ajouter les conditions de Neumann !   (mandatory)
//  (3) Ajouter les conditions en normal et tangentiel !   (strongly advised)
//  (4) Et remplacer le solveur plein par un truc un fifrelin plus subtil  (mandatory)

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

	double a = theProblem->A;
	double b = theProblem->B;
	double c = theProblem->C;
	double rho = theProblem->rho;
	double g = theProblem->g;
	double** A = theSystem->A;
	double* B = theSystem->B;

	int bandWidth = femMeshComputeBand(theMesh);
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
			double r = 0.0;

			for (i = 0; i < theSpace->n; i++) {
				dxdxsi += x[i] * dphidxsi[i];
				dxdeta += x[i] * dphideta[i];
				dydxsi += y[i] * dphidxsi[i];
				dydeta += y[i] * dphideta[i];
				r += x[i] * phi[i];
			}

			double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);

			for (i = 0; i < theSpace->n; i++) {
				dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
				dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
			}

			for (i = 0; i < theSpace->n; i++) {
				for (j = 0; j < theSpace->n; j++) {
					if (problemCase == AXISYM) {
						A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] * r +
							dphidy[i] * c * dphidy[j] * r + phi[i] * (b * dphidx[j] + a * phi[j] / r) +
							dphidx[i] * b * phi[j]) * jac * weight;
						A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] * r +
							dphidy[i] * c * dphidx[j] * r +
							phi[i] * b * dphidy[j]) * jac * weight;
						A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] * r +
							dphidx[i] * c * dphidy[j] * r +
							dphidy[i] * b * phi[j]) * jac * weight;
						A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] * r +
							dphidx[i] * c * dphidx[j] * r) * jac * weight;
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
					B[mapY[i]] -= phi[i] * g * rho * jac * weight * r;
				else
					B[mapY[i]] -= phi[i] * g * rho * jac * weight;
			}
		}
	}

	int* theConstrainedNodes = theProblem->constrainedNodes;
	for (int i = 0; i < theSystem->size; i++) {
		if (theConstrainedNodes[i] != -1) {
			double value = theProblem->conditions[theConstrainedNodes[i]]->value;
			femBoundaryType type = theProblem->conditions[theConstrainedNodes[i]]->type;

			if (type == DIRICHLET_X || type == DIRICHLET_Y) {
				femFullSystemConstrain(theSystem, i, value);
			}
		}
	}

	B = femFullSystemEliminate(theSystem);

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

	int taille = 3 + (theProblem->planarStrainStress != AXISYM ? 0 : 1);

	double* U = &UV[0];
	double* V = &UV[1];
	double* xsi = malloc(sizeof(double) * nLocal);
	double* eta = malloc(sizeof(double) * nLocal);
	double* dphidxsi = malloc(sizeof(double) * nLocal);
	double* dphideta = malloc(sizeof(double) * nLocal);
	double* dphidx = malloc(sizeof(double) * nLocal);
	double* dphidy = malloc(sizeof(double) * nLocal);
	femDiscreteXsi2(theSpace, xsi, eta);

	int* next = malloc(sizeof(int) * 10);
	double* e = malloc(sizeof(double) * taille);
	double* sigma = malloc(sizeof(double) * theNodes->nNodes * taille);

	for (int i = 0; i < theNodes->nNodes; i++) {
		int count = -1;
		int c;
		for (int j = 0; j < theMesh->nElem * nLocal; j++) {
			if (theMesh->elem[j] == i) {
				c = j % nLocal;
				next[count++] = j - c;
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
				dxdxsi += theNodes->X[theMesh->elem[next[j] + k]] * dphidxsi[k];
				dydxsi += theNodes->Y[theMesh->elem[next[j] + k]] * dphidxsi[k];
				dydeta += theNodes->Y[theMesh->elem[next[j] + k]] * dphideta[k];
				dxdeta += theNodes->X[theMesh->elem[next[j] + k]] * dphideta[k];
			}
			double A = dydeta * dxdxsi;
			double B = dydxsi * dxdeta;
			double J = fabs(A - B);
			for (int l = 0; l < nLocal; l++) {
				dphidx[l] = (dphidxsi[l] * dydeta - dphideta[l] * dydxsi) / J;
				dphidy[l] = (dphideta[l] * dxdxsi - dphidxsi[l] * dxdeta) / J;
			}
			for (int m = 0; m < nLocal; m++) {
				int index = next[j] + m;
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
			sigma[sigmaIndex + 4] = b * (e[0] + e[1]) + a * e[4];
		}
	}

	free(e);
	free(next);
	free(xsi);
	free(eta);
	free(dphidx);
	free(dphidy);
	free(dphidxsi);
	free(dphideta);

	return sigma;
}