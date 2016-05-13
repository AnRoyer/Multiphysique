#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <vector>
#include <map>
#include <mpi.h>
#include <string>
#include "gmshio.h"
#include "physicalio.h"
#include "fem.h"
#include "FE2.h"

using namespace std;

int main(int argc, char **argv)
{
    if(argc < 3)
    {
		// MPI initialization to avoid multiple displayong.
		MPI_Init(&argc, &argv);
		MPI_Status status;
		int nbproc, myrank ;
		MPI_Comm_rank( MPI_COMM_WORLD, &myrank);
		MPI_Comm_size( MPI_COMM_WORLD, &nbproc);

        if(myrank ==0) cout << "Usage: " << argv[0] << " file.msh file.phy" << endl;

		MPI_Finalize();		

        return 0;
    }

//----------------------------------------------------FLAG PERMETTANT DE CHOISIR ENTRE THERMIQUE OU ELECTRIQUE
	FemFlag thermalOrElectrical = THERMALFLAG;
//------------------------------------------------------------------------------------------------------------

//Il serait utile d'ajouter un paramètre dans le .phy (de type Nature par exemple) afin de pouvoir choisir entre thermique ou électrique.
//J'ai essayé mais je ne suis pas arriver à trouver comment faire pour lire un paramètre de type Nature.
    int natureFlag = 0;

    Type type;
    vector<Parameter*> parameters;
    Periodique conditions;
    Micro micro;
    double eps = 0;
	std::vector<int> methodFE2(2);

    Type type_micro;
    vector<Parameter*> parameters_micro;
    Periodique conditions_micro;
    Micro micro_micro;
    double eps_micro = 0;
	std::vector<int> methodFE2_micro(2); //no signification !!

    //lecture des PHYs
    readPHY(argv[2], parameters, conditions, micro, type, eps, methodFE2, natureFlag);

	if(type == DIRICHLET || type == PERIODIC || type == VONNEUMANN)
	{
		// MPI initialization to avoid multiple displayong.
		MPI_Init(&argc, &argv);
		MPI_Status status;
		int nbproc, myrank ;
		MPI_Comm_rank( MPI_COMM_WORLD, &myrank);
		MPI_Comm_size( MPI_COMM_WORLD, &nbproc);

		if(nbproc>1)
		{
			if(myrank == 0) cout << endl << "Error : FE1 method must be run with only 1 process !" << endl << endl;
			MPI_Finalize();
			return 0;
		}
		else
		{
			MPI_Finalize();
		}
		
	}

	int thermalOrElectricalFlag = 0;

	if(type == DIRICHLET || type == PERIODIC || type == VONNEUMANN)
	{
		while(1)
		{
		    int nbrchoix = 1;
		    std::map <int, FemFlag> choix;
			cout << endl;
		    cout << "What do you want:" << endl;

		    if((natureFlag & THERMALDATA) !=0)
		    {
		        cout << nbrchoix << ") \t Thermal computation" << endl;
		        choix[nbrchoix] = THERMALFLAG;
		        nbrchoix++;
		    }
		    if((natureFlag & ELECTRICALDATA) !=0)
		    {
		        cout << nbrchoix << ") \t Electrical computation" << endl;
		        choix[nbrchoix] = ELECTRICFLAG;
		        nbrchoix++;
		    }
		    /*if((natureFlag & THERMALDATA) !=0 && (natureFlag & ELECTRICALDATA) !=0 )
		    {
		        cout << nbrchoix << ") \t A coupling computation" << endl;
		        nbrchoix++;
		    }*/

		    cout << nbrchoix << ") \t Exit" << endl;
		    nbrchoix ++;

		    cout << "? >> ";
		    int userChoix;
		    cin >> userChoix;

		    if(userChoix == nbrchoix-1)
		    {
		        return 0;
		    }

		    if(userChoix <= 0 || userChoix >= nbrchoix)
		    {
		        cout << "Error : Invalid selection" << endl;
		    }
		    else
		    {
		        thermalOrElectrical = choix[userChoix];
				if (thermalOrElectrical == THERMALFLAG) thermalOrElectricalFlag = 1;
				else if (thermalOrElectrical == ELECTRICFLAG) thermalOrElectricalFlag = 2;
		        break;
		    }
		}//end while
	}// end if

	//if(thermalOrElectrical == THERMALFLAG) cout << "thermal !!!" << endl;
	//if(thermalOrElectrical == ELECTRICFLAG) cout << "electric !!!" << endl;

	if(type == DIRICHLET || type == PERIODIC || type == VONNEUMANN)
	{
		cout << endl;
		cout << "\t############################################################" << endl;
		cout << "\t############################################################" << endl;
		cout << "\t##                                                        ##" << endl;
		cout << "\t##                                                 11     ##" << endl;
		cout << "\t##                                               1111     ##" << endl;
		cout << "\t##     FFFFFFFFFFFFFFFF    EEEEEEEEEEEEEEEE     11 11     ##" << endl;
		cout << "\t##     FF                  EE                      11     ##" << endl;
		cout << "\t##     FF                  EE                      11     ##" << endl;
		cout << "\t##     FF                  EE                             ##" << endl;
		cout << "\t##     FF                  EE                             ##" << endl;
		cout << "\t##     FFFFFFFFFFFF        EEEEEEEEEEEE                   ##" << endl;
		cout << "\t##     FF                  EE                             ##" << endl;
		cout << "\t##     FF                  EE                             ##" << endl;
		cout << "\t##     FF                  EE                             ##" << endl;
		cout << "\t##     FF                  EE                             ##" << endl;
		cout << "\t##     FF                  EEEEEEEEEEEEEEEE               ##" << endl;
		cout << "\t##                                                        ##" << endl;
		cout << "\t############################################################" << endl;
		cout << "\t############################################################" << endl << endl << endl;
	}

    if(type == FE2withDIRICHLET || type == FE2withVONNEUMANN || type == FE2withPERIODIC)
    {
        int poubelle;
        readPHY(micro.filePhy.c_str(), parameters_micro, conditions_micro, micro_micro, type_micro, eps_micro, methodFE2_micro, poubelle);
    }

    vector<Node*> nodes;
    vector<Element*> elements;
    vector<Physical*> physicals;

    vector<Node*> nodes_micro;
    vector<Element*> elements_micro;
    vector<Physical*> physicals_micro;

    //lecture des MSHs
    readMSH(argv[1], nodes, elements, physicals);

    if(type == FE2withDIRICHLET || type == FE2withVONNEUMANN || type == FE2withPERIODIC)
    {
        readMSH(micro.fileMsh.c_str(), nodes_micro, elements_micro, physicals_micro);
    }

	if(type == DIRICHLET || type == PERIODIC || type == VONNEUMANN)
	{
    //Affichage des infos principales
		cout << endl;
		if(thermalOrElectrical == THERMALFLAG) cout << "SOLVING A THERMIC PROBLEM." << endl << endl;
		if(thermalOrElectrical == ELECTRICFLAG) cout << "SOLVING AN ELECTRIC PROBLEM." << endl << endl;
		//cout << "Read " << nodes.size() << " nodes and " << elements.size() << " elements." << endl << endl;
		if(type == DIRICHLET) cout << "Calling the FE1 method with DIRICHLET conditions." << endl;
		if(type == VONNEUMANN) cout << "Calling the FE1 method with VONNEUMANN conditions." << endl;
		if(type == PERIODIC) cout << "Calling the FE1 method with PERIODIC conditions." << endl;
		cout << endl;
		cout << "-----------------------------" << endl;
	}

	map<Node*, vector<double> > solutionTemperature;
	map<Node*, vector<double> > solutionFlux;
	map<Node*, vector<double> > solutionTemperature_micro;
	map<Node*, vector<double> > solutionFlux_micro;

	map<Node*, vector<double> > solutionPotential;
	map<Node*, vector<double> > solutionCurrent;
	map<Node*, vector<double> > solutionPotential_micro;
	map<Node*, vector<double> > solutionCurrent_micro;

	/*if(type == DIRICHLET || type == PERIODIC || type == VONNEUMANN)
	{
		for(unsigned int i = 0; i < parameters.size(); i++)
		{
		    if(parameters[i]->dim == 1 && (type == DIRICHLET || type == VONNEUMANN || type == FE2withDIRICHLET || type == FE2withVONNEUMANN))
		    {
		        if(parameters[i]->temperature != -1 && parameters[i]->voltage != -1)
		        {
		            cout << "Parameter " << parameters[i]->name << " has temperature of " << parameters[i]->temperature << " K and an electrical potential of " << parameters[i]->voltage << " V." << endl;
		        }
		        else if(parameters[i]->temperature != -1 && parameters[i]->voltage == -1)
		        {
		            cout << "Parameter " << parameters[i]->name << " has temperature of " << parameters[i]->temperature << " K." << endl;
		        }
		        else if(parameters[i]->temperature == -1 && parameters[i]->voltage != -1)
		        {
		            cout << "Parameter " << parameters[i]->name << " has an electrical potential of " << parameters[i]->voltage << " V." << endl;
		        }

		        if(parameters[i]->fluxTemperature != -1 && parameters[i]->fluxVoltage != -1)
		        {
		            cout << "Parameter " << parameters[i]->name << " has heat flux of " << parameters[i]->fluxTemperature << " W/m and an electrical flux of " << parameters[i]->fluxVoltage << "." << endl;
		        }
		        else if(parameters[i]->fluxTemperature != -1 && parameters[i]->fluxVoltage == -1)
		        {
		            cout << "Parameter " << parameters[i]->name << " has heat flux of " << parameters[i]->fluxTemperature << " W/m." << endl;
		        }
		        else if(parameters[i]->fluxTemperature == -1 && parameters[i]->fluxVoltage != -1)
		        {
		            cout << "Parameter " << parameters[i]->name << " has an electrical potential of " << parameters[i]->fluxVoltage << "." << endl;
		        }
		    }
		    else if(parameters[i]->dim == 2)
		    {
		        cout << "Parameter " << parameters[i]->name << " has:" << endl;

		        cout << "\t - A thermal conductivity [W/mK] of:" << endl;
		        for(unsigned int j = 0; j < parameters[i]->thermalConductivity.size(); j++)
		        {
		            cout << "\t\t" << parameters[i]->thermalConductivity[j]->name << " = " << endl;
		            cout << "\t\t|" << parameters[i]->thermalConductivity[j]->conductivity[0][0] << "\t" << parameters[i]->thermalConductivity[j]->conductivity[0][1] << "|" << endl;
		            cout << "\t\t|" << parameters[i]->thermalConductivity[j]->conductivity[1][0] << "\t" << parameters[i]->thermalConductivity[j]->conductivity[1][1] << "|" << endl << endl;
		        }

		        cout << "\t - A heat production of " << parameters[i]->thermalGeneration << " W/m^2" << endl;

		        cout << "\t - A electrical conductivity of :" << endl;
		        for(unsigned int j = 0; j < parameters[i]->electricalConductivity.size(); j++)
		        {
		            cout << "\t\t" << parameters[i]->electricalConductivity[j]->name << " = " << endl;
		            cout << "\t\t|" << parameters[i]->electricalConductivity[j]->conductivity[0][0] << "\t" << parameters[i]->electricalConductivity[j]->conductivity[0][1] << "|" << endl;
		            cout << "\t\t|" << parameters[i]->electricalConductivity[j]->conductivity[1][0] << "\t" << parameters[i]->electricalConductivity[j]->conductivity[1][1] << "|" << endl << endl;
		        }
		        cout << "\t - A change in charge density of " << parameters[i]->electricalGeneration << " A/m^3." << endl;
		    }
		}
	}

    if(type == PERIODIC)
    {
        cout << "There is periodic conditions with:" << endl;
        cout << "\t - A mean temperature: " << conditions.meanTemperature << endl;
        cout << "\t - A mean x gradient: " << conditions.xGradient << endl;
        cout << "\t - A mean y gradient: " << conditions.yGradient << endl;
    	cout << endl << endl;
    }*/

    //Initial guess of the temperature field, reading macro.msh and macro.phy. The initial guess will correspond to solutionTemperature_macro.pos and will be run in DIRICHLET (see the macro.phy)


	// Pour amméliorer => utilise MPI !!!!!!!!!!!!!!!!!!!!!!!!!!

	if(thermalOrElectrical == THERMALFLAG)
	{
		if(type == PERIODIC || type == FE2withPERIODIC)
		{
		    fem(nodes, elements, physicals, parameters, solutionTemperature, solutionFlux, THERMALFLAG, PERIODICFLAG, conditions, eps, type);
		}
		else if(type == DIRICHLET || type == FE2withDIRICHLET)
		{
		    fem(nodes, elements, physicals, parameters, solutionTemperature, solutionFlux, THERMALFLAG, DIRICHLETFLAG, conditions, eps, type);
		}
		else if(type == VONNEUMANN || type == FE2withVONNEUMANN)
		{
		    fem(nodes, elements, physicals, parameters, solutionTemperature, solutionFlux, THERMALFLAG, VONNEUMANNFLAG, conditions, eps, type);
		}
	}
	else if(thermalOrElectrical == ELECTRICFLAG)
	{
		if(type == PERIODIC || type == FE2withPERIODIC)
		{
		    fem(nodes, elements, physicals, parameters, solutionPotential, solutionCurrent, ELECTRICFLAG, PERIODICFLAG, conditions, eps, type);
		}
		else if(type == DIRICHLET || type == FE2withDIRICHLET)
		{
		    fem(nodes, elements, physicals, parameters, solutionPotential, solutionCurrent, ELECTRICFLAG, DIRICHLETFLAG, conditions, eps, type);
		}
		else if(type == VONNEUMANN || type == FE2withVONNEUMANN)
		{
		    fem(nodes, elements, physicals, parameters, solutionPotential, solutionCurrent, ELECTRICFLAG, VONNEUMANNFLAG, conditions, eps, type);
		}
	}

    //FE2 method.
    if(type == FE2withDIRICHLET || type == FE2withVONNEUMANN || type == FE2withPERIODIC)
    {
        FE2(nodes_micro, elements_micro, physicals_micro, parameters_micro, solutionTemperature_micro, solutionFlux_micro,
            conditions_micro, nodes, elements, physicals, parameters, solutionTemperature, solutionFlux, eps,
			methodFE2, THERMALFLAG, type, argc, argv, natureFlag);
    }
    /*else if(thermalOrElectrical == ELECTRICFLAG && (type == FE2withDIRICHLET || type == FE2withVONNEUMANN || type == FE2withPERIODIC))
    {
        FE2(nodes_micro, elements_micro, physicals_micro, parameters_micro, solutionPotential_micro, solutionCurrent_micro,
            conditions_micro, nodes, elements, physicals, parameters, solutionPotential, solutionCurrent, eps,
			methodFE2, ELECTRICFLAG, type, argc, argv);
    }*/

	if(type == DIRICHLET || type == VONNEUMANN || type == PERIODIC)
	{
		if(thermalOrElectrical == THERMALFLAG)
		{
			writeMSH((char*)"solutionTemperature.pos", solutionTemperature);
			writeMSH((char*)"solutionFlux.pos", solutionFlux);
		}
		else if(thermalOrElectrical == ELECTRICFLAG)
		{
			writeMSH((char*)"solutionPotential.pos", solutionPotential);
			writeMSH((char*)"solutionCurrent.pos", solutionCurrent);
		}
		cout << endl;
		cout << "-----------------------------" << endl;
		if(type == DIRICHLET && thermalOrElectrical == THERMALFLAG) 
		cout << "The THERMIC problem has been solved in ONE SCALE with DIRICHLET conditions." << endl;
		if(type == VONNEUMANN && thermalOrElectrical == THERMALFLAG) 
		cout <<"The THERMIC problem has been solved in ONE SCALE with VON NEUMANN conditions." << endl;
		if(type == PERIODIC && thermalOrElectrical == THERMALFLAG) 
		cout <<  "The THERMIC problem has been solved in ONE SCALE with PERIODIC conditions."  << endl;
		if(type == DIRICHLET && thermalOrElectrical == ELECTRICFLAG) 
		cout << "The ELECTRIC problem has been solved in ONE SCALE with DIRICHLET conditions." << endl;
		if(type == VONNEUMANN && thermalOrElectrical == ELECTRICFLAG) 
		cout <<"The ELECTRIC problem has been solved in ONE SCALE with VON NEUMANN conditions." << endl;
		if(type == PERIODIC && thermalOrElectrical == ELECTRICFLAG) 
		cout <<  "The ELECTRIC problem has been solved in ONE SCALE with PERIODIC conditions."  << endl;
		cout << endl;
		/*cout << "Press Enter to show the solution." << endl;
		std::cin.ignore();
		system("gmsh l.msh solutionTemperature.pos &");	*/
	}
    return 0;
}



