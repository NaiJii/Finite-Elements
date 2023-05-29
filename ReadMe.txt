Ce projet a été configuré sur Windows principalement, il est possible que vous ayez des problèmes de compilation sur Linux 
-bien que nous ayons fait de notre mieux pour assurer la cross-compatibility et cela ne devrait pas arriver.

1. Rajouter GMSH dans chaque projet (Project, Post/Pre-Processing) et copier le dll de gmsh dans le répertoire de sortie si nécessaire.

2. Générer et compiler ProjectPreProcessing, le fichier éxécutable devrait avoir le chemin: "Finite-Elements\ProjectPreProcessor\out\build\x64-Debug\myFem.exe"

3. Générer et compiler Project, le fichier éxécutable devrait avoir le chemin: "Finite-Elements\out\build\x64-Debug\myFem.exe"

4. Générer et compiler ProjectPostProcessing, le fichier éxécutable devrait avoir le chemin: "Finite-Elements\ProjectPostProcessor\out\build\x64-Debug\myFem.exe"

5. Différentes touches (D,V,K) peuvent être utilisées pour changer le mode d'affichage.

Si vous voulez juste observer le résultat final, seul les points 1., 4. et 5. sont nécessaires.

Notre simulation porte sur une bonbonne de gaz subissant une pression sur toutes ses parois internes (en négligeant la pression externe).
Vous pourrez trouver les détails de pression dans notre main.c de PreProcessing et dans le rapport.

Les différents détails d'implémentation sont disponibles dans homework.c du Projet.

Merci d'avoir lu ce README, nous espérons que vous apprécierez notre projet.