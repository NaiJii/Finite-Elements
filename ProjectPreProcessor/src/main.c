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

#define OUR_GEO

int main(void)
{
	//
	//  -1- Construction de la geometrie
	//

	double Lx = 1.0;
	double Ly = 1.0;
	geoInitialize();
	femGeo* theGeometry = geoGetGeometry();

	theGeometry->LxPlate = Lx;
	theGeometry->LyPlate = Ly;
	theGeometry->h = Lx * 0.05;
	theGeometry->elementType = FEM_TRIANGLE;

	//    geoMeshGenerate();      // Utilisation de OpenCascade

	//  geoMeshGenerateGeo();   // Utilisation de outils de GMSH
								// Attention : les entités sont différentes !
								// On a aussi inversé la géomtrie pour rire !

#ifdef OUR_GEO
	geoMeshGenerateGeoFile("../../../data/our_mesh.geo");   // Lecture fichier geo
#else
	geoMeshGenerateGeoFile("../../../../data/mesh.geo");   // Lecture fichier geo
#endif

	geoMeshImport();
#ifdef OUR_GEO
	//geoMeshGenerateGeoFile("../../../../data/our_mesh.geo");   // Lecture fichier geo
	geoSetDomainName(4, "Bottom");
	geoSetDomainName(5, "Top");
	geoSetDomainName(6, "LeftIn");
	geoSetDomainName(7, "RightIn");
	geoSetDomainName(8, "BottomIn");
	geoSetDomainName(9, "TopIn");
	geoSetDomainName(0, "Left");
	geoSetDomainName(2, "Right");
	geoMeshWrite("../../../../data/mesh.txt");
#else
	//geoMeshGenerateGeoFile("../../../../data/mesh.geo");   // Lecture fichier geo
	geoSetDomainName(0, "Symetry");
	geoSetDomainName(1, "Bottom");
	geoMeshWrite("../../../../data/mesh.txt");
#endif

	//
	//  -2- Definition du probleme
	//

	double E = 220.e9;
	double nu = 0.3;
	double rho = 7.85e3;
	double g = 9.81;
	femProblem* theProblem = femElasticityCreate(theGeometry, E, nu, rho, g, PLANAR_STRAIN);

#ifdef OUR_GEO
	/* 
	 When should I use Dirichlet boundary conditions(displacements) ?
	 When should I use Neumann boundary conditions (forces) ?
	 What is the difference between Dirichlet and Neumann boundary conditions ?
	 The Dirichlet boundary conditions are used to impose the displacement of the structure on the boundary.
	 The Neumann boundary conditions are used to impose the forces on the boundary.
	 In the case of a gas bottle,
	 We want no displacement on the bottom of the bottle (Dirichlet boundary conditions)
	*/

	/// femElasticityAddBoundaryCondition(theProblem, "Bottom", DIRICHLET_X, 0.0);
	femElasticityAddBoundaryCondition(theProblem, "Bottom", DIRICHLET_Y, 0.0);

	/*
	Avg pressure on a gaz tank is 10 bar at around 30 - 35 °C
	(https://www.obelink.fr/conseils/pourquoi-un-detendeur-de-gaz-sur-une-bouteille-de-gaz.html#:~:text=En%20r%C3%A8gle%20g%C3%A9n%C3%A9rale%2C%20la%20pression,de%2030%20ou%2050%20millibars*.)
	10 bar == 1e6 N/m^2

	bottle:
	- width = 30.6 cm
	- height = 147 cm
	- thickness = 7.65 cm
	- radius = width / 2 = 15.3 cm

	for the top/bottom part (symmetrical): width - thickness * 2 = 15.3 cm 
	-> 1e6 * 0.153 = 153 [kN/m]  
	for the left/right part (symmetrical): height - radius * 2 = 101.1 cm
	-> 1e6 * 1.011 = 1011 [kN/m]
	Nos forces sont en /m car on est dans un problème 2D (on a pas de profondeur)
	*/
	
	double pressureSide = 1011e3;
	double pressureTop = 153e3;
	// The gas is pushing from within the bottle (Neumann boundary conditions)
	femElasticityAddBoundaryCondition(theProblem, "LeftIn", NEUMANN_N, -pressureSide); // Towards the inside.
	femElasticityAddBoundaryCondition(theProblem, "RightIn", NEUMANN_N, -pressureSide); // Towards the inside.
	femElasticityAddBoundaryCondition(theProblem, "TopIn", NEUMANN_N, -pressureTop); // Towards the inside.
	femElasticityAddBoundaryCondition(theProblem, "BottomIn", NEUMANN_N, -pressureTop); // Towards the inside.

	femElasticityPrint(theProblem);
	femElasticityWrite(theProblem, "../../../../data/problem.txt");
#else
	femElasticityAddBoundaryCondition(theProblem, "Symetry", DIRICHLET_X, 0.0);
	femElasticityAddBoundaryCondition(theProblem, "Bottom", DIRICHLET_Y, 0.0);
	femElasticityPrint(theProblem);
	femElasticityWrite(theProblem, "../../../../data/problem.txt");
#endif

	//
	//  -3- Champ de la taille de référence du maillage
	//

	double* meshSizeField = malloc(theGeometry->theNodes->nNodes * sizeof(double));
	femNodes* theNodes = theGeometry->theNodes;
	for (int i = 0; i < theNodes->nNodes; ++i)
		meshSizeField[i] = theGeometry->geoSize(theNodes->X[i], theNodes->Y[i]);
	double hMin = femMin(meshSizeField, theNodes->nNodes);
	double hMax = femMax(meshSizeField, theNodes->nNodes);
	printf(" ==== Global requested h : %14.7e \n", theGeometry->h);
	printf(" ==== Minimum h          : %14.7e \n", hMin);
	printf(" ==== Maximum h          : %14.7e \n", hMax);

	//
	//  -4- Visualisation
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
		if (glfwGetKey(window, 'N') == GLFW_PRESS && freezingButton == FALSE) { domain++; freezingButton = TRUE; told = t; }

		if (t - told > 0.5) { freezingButton = FALSE; }
		if (mode == 1) {
			glfemPlotField(theGeometry->theElements, meshSizeField);
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

	free(meshSizeField);
	femElasticityFree(theProblem);
	geoFree();
	glfwTerminate();

	exit(EXIT_SUCCESS);
	return 0;
}