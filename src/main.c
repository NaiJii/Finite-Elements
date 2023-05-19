/*
 *  main.c
 *  Projet 2022-2023
 *  Elasticite lineaire plane
 *
 *  Code de calcul
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#include "fem.h"

int main(char argc, char* argv[]) {
	if (argc != 1 && argc != 3)
		Error("Usage : main [meshfile problemfile]");

	femGeo* theGeometry = geoGetGeometry();
	femProblem* theProblem = NULL;

	if (argc == 1) {
		geoMeshRead("../../../data/mesh.txt");
		theProblem = femElasticityRead(theGeometry, "../../../data/problem.txt");
	}
	else {
		geoMeshRead(argv[1]);
		theProblem = femElasticityRead(theGeometry, argv[2]);
	}

	if (!theProblem)
		Error("No problem!(or maybe yes)");

	femElasticityPrint(theProblem);
	double* theSoluce = femElasticitySolve(theProblem);
	femNodes* theNodes = theGeometry->theNodes;
	femFieldWrite(theNodes->nNodes, 2, &theSoluce[0], "../../../data/U.txt");
	femFieldWrite(theNodes->nNodes, 2, &theSoluce[1], "../../../data/V.txt");
	femElasticityFree(theProblem);
	geoFree();
	return 0;
}