16/05/06 README

Utilisation du code.

Le code peut être utilisé soit à 1 échelle soit avec les deux échelles.

1 échelle:

	mpirun -np 1 ./build/MP domain_macro.msh domain_macro.phy

	Utiliser 1 process pour éviter les affichages multiples.

	conditions: marche en dirichlet et periodique. Il faut changer le type dans domain_macro.phy

	dirichlet: type="dirichlet" 
	periodique: type="periodic"

	Il faut ensuite changer les paramètres correspondants dans le .phy
	Regarder dans l'exemple "l.phy"

	!!! Si des conditions périodiques sont utilisées il faut un domaine macroscopique rectangulaire (ici le L ne marchera pas). !!!

2 échelles:

	mpirun -np 4 ./build/MP domain_macro.msh domain_macro.phy

	Peut utiliser plusieurs process (ici 4).

	conditions: marche en dirichlet. Il faut changer le type dans domain_macro.phy

	dirichlet: type="fe2D" 

	Il faut ensuite changer les paramètres correspondants dans le .phy
	Regarder dans l'exemple "l.phy"

	le domaine micro est compris dans micro.msh et micro.phy et est automatiquement lu par le programme 
	(les fichiers dovent se trouver dans le répertoire courrant)

	Pour changer les paramètres micro (toujours périodiques) modifier micro.phy

