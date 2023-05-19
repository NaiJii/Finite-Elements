#include "fem.h"
#include <math.h>

// Il faut un fifrelin generaliser ce code.....
//  (1) Ajouter l'axisymétrique !    (mandatory)
//  (2) Ajouter les conditions de Neumann !   (mandatory)
//  (3) Ajouter les conditions en normal et tangentiel !   (strongly advised)
//  (4) Et remplacer le solveur plein par un truc un fifrelin plus subtil  (mandatory)

// Note:
// fem.c est légèrement modifié, pour ajouter les conditions nécessaires
// aussi, 
// femComputeBoundaryCondition(theProblem, theBoundary)
// est appelé dans femElasticityAddBoundaryCondition.

#define SOLVER_BAND

double* femTemp = NULL;

int femCompare(const void* a, const void* b) {
	const int ia = *(const int*)a;
	const int ib = *(const int*)b;
	const double d = femTemp[ia] - femTemp[ib];
	if (d < 0) return 1;
	if (d > 0) return -1;
	return 0;
}

int femMeshComputeBand(femMesh* theMesh) {
	int iElem, j, myMax, myMin, myBand, map[4];
	memset(map, 0, 4 * sizeof(int));
	int nLocal = theMesh->nLocalNode;
	myBand = 0;
	for (iElem = 0; iElem < theMesh->nElem; iElem++) {
		for (j = 0; j < nLocal; ++j)
			map[j] = theMesh->number[theMesh->elem[iElem * nLocal + j]];
		myMin = map[0];
		myMax = map[0];
		for (j = 1; j < nLocal; j++) {
			myMax = fmax(map[j], myMax);
			myMin = fmin(map[j], myMin);
		}
		if (myBand < (myMax - myMin)) myBand = myMax - myMin;
	}

	return (++myBand) * 2;
}

void femBandSystemInit(femBandSystem* myBandSystem)
{
	int i;
	int size = myBandSystem->size;
	int band = myBandSystem->band;
	for (i = 0; i < size * (band + 1); i++)
		myBandSystem->B[i] = 0;
}

femBandSystem* femBandSystemCreate(int size, int band)
{
	femBandSystem* myBandSystem = malloc(sizeof(femBandSystem));
	myBandSystem->B = malloc(sizeof(double) * size * (band + 1));
	myBandSystem->A = malloc(sizeof(double*) * size);
	myBandSystem->size = size;
	myBandSystem->band = band;
	myBandSystem->A[0] = myBandSystem->B + size;
	int i;
	for (i = 1; i < size; i++)
		myBandSystem->A[i] = myBandSystem->A[i - 1] + band - 1;
	femBandSystemInit(myBandSystem);
	return (myBandSystem);
}

double* femBandSystemEliminate(femBandSystem* myBand)
{
	double** A, * B, factor;
	int     i, j, k, jend, size, band;
	A = myBand->A;
	B = myBand->B;
	size = myBand->size;
	band = myBand->band;

	/* Incomplete Cholesky factorization */

	for (k = 0; k < size; k++) {
		if (fabs(A[k][k]) <= 1e-4) {
			Error("Cannot eleminate with such a pivot");
		}
		jend = fmin(k + band, size);
		for (i = k + 1; i < jend; i++) {
			factor = A[k][i] / A[k][k];
			for (j = i; j < jend; j++)
				A[i][j] = A[i][j] - A[k][j] * factor;
			B[i] = B[i] - B[k] * factor;
		}
	}

	/* Back-substitution */

	for (i = (size - 1); i >= 0; i--) {
		factor = 0;
		jend = fmin(i + band, size);
		for (j = i + 1; j < jend; j++)
			factor += A[i][j] * B[j];
		B[i] = (B[i] - factor) / A[i][i];
	}

	return(myBand->B);
}

void femMeshRenumber(femProblem* theProblem, femRenumType renumType) {
	int i, * inverse;
	femMesh* theMesh = theProblem->geometry->theElements;

	switch (renumType) {
	case FEM_NO:
		for (i = 0; i < theMesh->nodes->nNodes; i++)
			theMesh->number[i] = i;
		break;
	case FEM_XNUM:
		inverse = malloc(sizeof(int) * theMesh->nodes->nNodes);
		for (i = 0; i < theMesh->nodes->nNodes; i++)
			inverse[i] = i;
		femTemp = theMesh->nodes->X;
		qsort(inverse, theMesh->nodes->nNodes, sizeof(int), femCompare);
		for (i = 0; i < theMesh->nodes->nNodes; i++)
			theMesh->number[inverse[i]] = i;
		free(inverse);
		break;
	case FEM_YNUM:
		inverse = malloc(sizeof(int) * theMesh->nodes->nNodes);
		for (i = 0; i < theMesh->nodes->nNodes; i++)
			inverse[i] = i;
		femTemp = theMesh->nodes->Y;
		qsort(inverse, theMesh->nodes->nNodes, sizeof(int), femCompare);
		for (i = 0; i < theMesh->nodes->nNodes; i++)
			theMesh->number[inverse[i]] = i;
		free(inverse);
		break;
	default: Error("Unexpected renumbering option");
	}
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

	double a = theProblem->A;
	double b = theProblem->B;
	double c = theProblem->C;
	double rho = theProblem->rho;
	double g = theProblem->g;
	double** A = theSystem->A;
	double* B = theSystem->B;

	int bandWidth = femMeshComputeBand(theMesh);

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

			if (type >= DIRICHLET_X && type <= DIRICHLET_T) 
				femFullSystemConstrain(theSystem, i, value);
		}
	}

	double* renumB = malloc(theSystem->size * sizeof(double));
	if (!renumB)
		Error("Allocation error");

	// pivots
	for (int k = 0; k < theSystem->size / 2; k++) {
		if ((A[2 * k][2 * k] == 0) && (A[2 * k + 1][2 * k + 1] == 0)) {
			printf("Fixing disconnected node %d\n", k);
			A[2 * k][2 * k] = 1;
			A[2 * k + 1][2 * k + 1] = 1;
		}
	}

#ifdef SOLVER_BAND
	bandWidth = bandWidth * 2 + 1;
	femBandSystem* theBandSystem = femBandSystemCreate(theSystem->size, bandWidth);
	printf("Solving Ax=b using band solver\n");

	for (i = 0; i < theSystem->size; i++) {
		int jmin = (0 > i - bandWidth / 2) ? 0 : i - bandWidth / 2;
		int jmax = (theSystem->size < i + bandWidth / 2 + 1) ? theSystem->size : i + bandWidth / 2 + 1;
		for (j = jmin; j < jmax; j++) {
			theBandSystem->A[i][j] = A[i][j];
		}

		theBandSystem->B[i] = theSystem->B[i];
	}

	theBandSystem->B = femBandSystemEliminate(theBandSystem);
	memcpy(renumB, theBandSystem->B, theSystem->size * sizeof(double));
#else 
	printf("Solving Ax=b using full solver\n");
	B = femFullSystemEliminate(theSystem);
	memcpy(renumB, B, theSystem->size * sizeof(double));
#endif

	for (i = 0; i < theMesh->nodes->nNodes; i++) {
		B[2 * i] = renumB[2 * theMesh->number[i]];
		B[2 * i + 1] = renumB[2 * theMesh->number[i] + 1];
	}

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