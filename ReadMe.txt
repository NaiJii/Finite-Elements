Ce projet a �t� configur� sur Windows principalement, il est possible que vous ayez des probl�mes de compilation sur Linux 
-bien que nous ayons fait de notre mieux pour assurer la cross-compatibility et cela ne devrait pas arriver.

1. Rajouter GMSH dans chaque projet (Project, Post/Pre-Processing) et copier le dll de gmsh dans le r�pertoire de sortie si n�cessaire.

2. G�n�rer et compiler ProjectPreProcessing, le fichier �x�cutable devrait avoir le chemin: "Finite-Elements\ProjectPreProcessor\out\build\x64-Debug\myFem.exe"

3. G�n�rer et compiler Project, le fichier �x�cutable devrait avoir le chemin: "Finite-Elements\out\build\x64-Debug\myFem.exe"

4. G�n�rer et compiler ProjectPostProcessing, le fichier �x�cutable devrait avoir le chemin: "Finite-Elements\ProjectPostProcessor\out\build\x64-Debug\myFem.exe"

5. Diff�rentes touches (D,V,K) peuvent �tre utilis�es pour changer le mode d'affichage.

Si vous voulez juste observer le r�sultat final, seul les points 1., 4. et 5. sont n�cessaires.

Notre simulation porte sur une bonbonne de gaz subissant une pression sur toutes ses parois internes (en n�gligeant la pression externe).
Vous pourrez trouver les d�tails de pression dans notre main.c de PreProcessing et dans le rapport.

Les diff�rents d�tails d'impl�mentation sont disponibles dans homework.c du Projet.

Merci d'avoir lu ce README, nous esp�rons que vous appr�cierez notre projet.