/*
 *  main.c
 *  Projet 2022-2023
 *  Elasticite lineaire plane
 *
 *  Preprocesseur
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#include "glfem.h"

int main(void)
{
	printf("\n\n    V : Mesh and displacement norm \n");
	printf("    D : Domains \n");
	printf("    N : Next domain highlighted\n\n\n");

	//
	//  -1- Lecture des donnees
	//

	femGeo* theGeometry = geoGetGeometry();
	geoMeshRead("../../../../data/mesh.txt");

	femProblem* theProblem = femElasticityRead(theGeometry, "../../../../data/problem.txt");
	double* theSoluce = theProblem->system->B;
	int n = theGeometry->theNodes->nNodes;

	femFieldRead(&n, 2, &theSoluce[0], "../../../../data/U.txt");
	femFieldRead(&n, 2, &theSoluce[1], "../../../../data/V.txt");
	femElasticityPrint(theProblem);

	//
	//  -2- Deformation du maillage pour le plot final
	//      Creation du champ de la norme du deplacement
	//		Tout en gardant le maillage original
	//

	femNodes* theNodes = theGeometry->theNodes;
	// original
	double* X = malloc(n * sizeof(double));
	double* Y = malloc(n * sizeof(double));
	memcpy(X, theNodes->X, n * sizeof(double));
	memcpy(Y, theNodes->Y, n * sizeof(double));
	// deformed 
	double* Xdef = malloc(n * sizeof(double));
	double* Ydef = malloc(n * sizeof(double));

	double deformationFactor = 4e5;
	double* normDisplacement = malloc(n * sizeof(double));

	for (int i = 0; i < n; i++) {
		theNodes->X[i] += theSoluce[2 * i + 0] * deformationFactor;
		theNodes->Y[i] += theSoluce[2 * i + 1] * deformationFactor;
		normDisplacement[i] = sqrt(theSoluce[2 * i + 0] * theSoluce[2 * i + 0] +
			theSoluce[2 * i + 1] * theSoluce[2 * i + 1]);
	}

	double hMin = femMin(normDisplacement, n);
	double hMax = femMax(normDisplacement, n);
	memcpy(Xdef, theNodes->X, n * sizeof(double));
	memcpy(Ydef, theNodes->Y, n * sizeof(double));

	double* u = calloc(n, sizeof(double));

	printf(" ==== Minimum displacement          : %14.7e \n", hMin);
	printf(" ==== Maximum displacement          : %14.7e \n", hMax);

	//
	//  -3- Visualisation
	//

	int mode = 1;
	int domain = 0;
	int freezingButton = FALSE;
	double t, told = 0;
	char theMessage[MAXNAME];

	GLFWwindow* window = glfemInit("EPL1110 : Project 2022-23 ");
	glfwMakeContextCurrent(window);

	do {
		int w, h;
		glfwGetFramebufferSize(window, &w, &h);
		glfemReshapeWindows(theGeometry->theNodes, w, h);

		t = glfwGetTime();
		if (glfwGetKey(window, 'D') == GLFW_PRESS) { mode = 0; }
		if (glfwGetKey(window, 'V') == GLFW_PRESS) { mode = 1; }
		if (glfwGetKey(window, 'O') == GLFW_PRESS) { mode = 2; }
		if (glfwGetKey(window, 'S') == GLFW_PRESS) { mode = 3; }
		if (glfwGetKey(window, 'K') == GLFW_PRESS) { mode = 4; }
		if (glfwGetKey(window, 'L') == GLFW_PRESS) { mode = 5; }
		if (glfwGetKey(window, 'N') == GLFW_PRESS && freezingButton == FALSE) { domain++; freezingButton = TRUE; told = t; }

		if (t - told > 0.5) { freezingButton = FALSE; }
		if (mode == 4 || mode == 5) {
			const int frames = 40;
			double* delta = malloc(n * sizeof(double));
			double* field = normDisplacement;
			double min = hMin;
			double max = hMax;

			if (mode == 5) {
				/*
				* field = stress;
								max = hMax;
				min = hMin;
				*/
			}

			for (int i = 0; i < n; i++) {
				theGeometry->theNodes->X[i] = X[i];
				theGeometry->theNodes->Y[i] = Y[i];
				delta[i] = field[i] / frames;
				field[i] = 0;
			}

			const double field0 = field[0];
			for (int i = 0; i < frames; i++) {
				field[0] = max;

				glfemPlotField(theGeometry->theElements, field);
				glfemPlotMesh(theGeometry->theElements);
				glColor3f(1.0, 0.0, 0.0); glfemMessage(theMessage);
				_sleep(100);
				glfwSwapBuffers(window);
				glfwPollEvents();
				glfwGetFramebufferSize(window, &w, &h);
				glfemReshapeWindows(theGeometry->theNodes, w, h);

				for (int j = 0; j < n; j++) {
					theGeometry->theNodes->X[j] += theSoluce[2 * j + 0] * deformationFactor / frames;
					theGeometry->theNodes->Y[j] += theSoluce[2 * j + 1] * deformationFactor / frames;
					field[j] += delta[j];
				}
			}
			field[0] = field0;
			free(delta);

			mode -= 3;
		}
		if (mode == 3) {
			//glfemPlotField(theGeometry->theElements, stress);
			sprintf(theMessage, "Number of elements : %d ", theGeometry->theElements->nElem);
			glColor3f(1.0, 0.0, 0.0); glfemMessage(theMessage);
		}
		if (mode == 2) {
			memcpy(theGeometry->theElements->nodes->X, X, n * sizeof(double));
			memcpy(theGeometry->theElements->nodes->Y, Y, n * sizeof(double));
			glfemPlotField(theGeometry->theElements, u);
			glfemPlotMesh(theGeometry->theElements);
			memcpy(theGeometry->theNodes->X, Xdef, n * sizeof(double));
			memcpy(theGeometry->theNodes->Y, Ydef, n * sizeof(double));
			sprintf(theMessage, "Number of elements : %d ", theGeometry->theElements->nElem);
			glColor3f(1.0, 0.0, 0.0); glfemMessage(theMessage);
		}
		if (mode == 1) {
			glfemPlotField(theGeometry->theElements, normDisplacement);
			glfemPlotMesh(theGeometry->theElements);
			sprintf(theMessage, "Number of elements : %d ", theGeometry->theElements->nElem);
			glColor3f(1.0, 0.0, 0.0); glfemMessage(theMessage);
		}
		if (mode == 0) {
			domain = domain % theGeometry->nDomains;
			glfemPlotDomain(theGeometry->theDomains[domain]);
			sprintf(theMessage, "%s : %d ", theGeometry->theDomains[domain]->name, domain);
			glColor3f(1.0, 0.0, 0.0); glfemMessage(theMessage);
		}

		glfwSwapBuffers(window);
		glfwPollEvents();
	} while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS &&
		glfwWindowShouldClose(window) != 1);

	// Check if the ESC key was pressed or the window was closed

	free(normDisplacement);
	free(X);
	free(Y);
	free(Xdef);
	free(Ydef);
	free(u);
	femElasticityFree(theProblem);
	geoFree();
	glfwTerminate();

	exit(EXIT_SUCCESS);
	return 0;
}