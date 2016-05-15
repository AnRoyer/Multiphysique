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
		// MPI initialization to avoid multiple displaying.
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
	FemFlag thermalOrElectrical = DEFAULTFLAG;
//------------------------------------------------------------------------------------------------------------

    int natureFlag = 0;
	Type type = DEFAULTTYPE;
	int nature;
	int nature_micro;

    Type type_thermic = DEFAULTTYPE;
    Type type_electric = DEFAULTTYPE;
    vector<Parameter*> parameters;
    Periodique conditions;
    Micro micro;
    double eps = 0;
	std::vector<int> methodFE2(2);

    Type type_micro_thermic;
    Type type_micro_electric;
    vector<Parameter*> parameters_micro;
    Periodique conditions_micro;
    Micro micro_micro;
    double eps_micro = 0;
	std::vector<int> methodFE2_micro(2); //no signification !!

    //lecture des PHYs
    readPHY(argv[2], parameters, conditions, micro, type_thermic, type_electric, nature, eps, methodFE2, natureFlag);
	/*for(int i=0;i<9;i++)
	{
		cout << parameters[i]->name << endl;
	}	*/
	int FlagFE2 = 0;

	//cout << "nature " << nature;

	if(nature == 1)//thermic
	{
		thermalOrElectrical = THERMALFLAG;
		type = type_thermic;
		if(type == FE2withDIRICHLET) FlagFE2 = 1;
	}
	else if(nature == 2)//electric
	{
		thermalOrElectrical = ELECTRICFLAG;
		type = type_electric;
		if(type == FE2withDIRICHLET) FlagFE2 = 1;
	}
	else if(nature == 3)//couppled
	{
		if(type_electric == FE2withDIRICHLET || type_thermic == FE2withDIRICHLET) FlagFE2 = 1;
	}

	//cout << "FlagFE2 " << FlagFE2;
	
	if(FlagFE2 ==0)
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


	if(FlagFE2 ==0)
	{
		if(thermalOrElectrical == ELECTRICFLAG)
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
			cout << "\t##     FFFFFFFFFFFF        EEEEEEEEEEEE             ,/    ##" << endl;
			cout << "\t##     FF                  EE                     ,'/     ##" << endl;
			cout << "\t##     FF                  EE                   ,' /      ##" << endl;
			cout << "\t##     FF                  EE                 ,'  /_____, ##" << endl;
			cout << "\t##     FF                  EE               .'____    ,'  ##" << endl;
			cout << "\t##     FF                  EEEEEEEEEEEEEEEE      /  ,'    ##" << endl;
			cout << "\t##                                              / ,'      ##" << endl;
			cout << "\t###############################################/,' #########" << endl;
			cout << "\t##############################################/'############" << endl << endl << endl;
		}
		else if(thermalOrElectrical == THERMALFLAG)
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
			cout << "\t##     FF                  EE                     (       ##" << endl;
			cout << "\t##     FFFFFFFFFFFF        EEEEEEEEEEEE           ))      ##" << endl;
			cout << "\t##     FF                  EE                     {_}     ##" << endl;
			cout << "\t##     FF                  EE                    .-;-.    ##" << endl;
			cout << "\t##     FF                  EE                   |'-=-'|   ##" << endl;
			cout << "\t##     FF                  EE                   |     |   ##" << endl;
			cout << "\t##     FF                  EEEEEEEEEEEEEEEE     |     |   ##" << endl;
			cout << "\t##                                              |     |   ##" << endl;
			cout << "\t################################################|     |#####" << endl;
			cout << "\t################################################'.___.'#####" << endl << endl << endl; 
    	} 
	}

    if(FlagFE2 ==1)
    {
        int poubelle;
		int nature_micro;
        readPHY(micro.filePhy.c_str(), parameters_micro, conditions_micro, micro_micro, type_micro_thermic, type_micro_electric, nature_micro, eps_micro, methodFE2_micro, poubelle);
    }

    vector<Node*> nodes;
    vector<Element*> elements;
    vector<Physical*> physicals;

    vector<Node*> nodes_micro;
    vector<Element*> elements_micro;
    vector<Physical*> physicals_micro;

    //lecture des MSHs
    readMSH(argv[1], nodes, elements, physicals);

    if(FlagFE2 ==1)
    {
        readMSH(micro.fileMsh.c_str(), nodes_micro, elements_micro, physicals_micro);
    }

	if(FlagFE2 ==0)
	{
    //Affichage des infos principales
		cout << endl;
		//if(thermalOrElectrical == THERMALFLAG) cout << "SOLVING A THERMIC PROBLEM." << endl << endl;
		//else if(thermalOrElectrical == ELECTRICFLAG) cout << "SOLVING AN ELECTRIC PROBLEM." << endl << endl;
		//cout << "Read " << nodes.size() << " nodes and " << elements.size() << " elements." << endl << endl;
		if(type == DIRICHLET) cout << "Calling the FE1 method with DIRICHLET conditions." << endl;
		else if(type == VONNEUMANN) cout << "Calling the FE1 method with VONNEUMANN conditions." << endl;
		else if(type == PERIODIC) cout << "Calling the FE1 method with PERIODIC conditions." << endl;
		cout << endl;
		cout << "-----------------------------" << endl;
	}

	map<Node*, vector<double> > solutionTemperature;
	map<Node*, vector<double> > solutionFlux;

	map<Node*, vector<double> > solutionPotential;
	map<Node*, vector<double> > solutionCurrent;

	// Pour amméliorer => utilise MPI !!!!!!!!!!!!!!!!!!!!!!!!!!
	if(FlagFE2 == 0 && nature !=3)
	{
		if(thermalOrElectrical == THERMALFLAG)
		{
			if(type == PERIODIC)
			{
				fem(nodes, elements, physicals, parameters, solutionTemperature, solutionFlux, THERMALFLAG, PERIODICFLAG, conditions, eps, type, 0);
			}
			else if(type == DIRICHLET)
			{
				fem(nodes, elements, physicals, parameters, solutionTemperature, solutionFlux, THERMALFLAG, DIRICHLETFLAG, conditions, eps, type, 0);
			}
			else if(type == VONNEUMANN)
			{
				fem(nodes, elements, physicals, parameters, solutionTemperature, solutionFlux, THERMALFLAG, VONNEUMANNFLAG, conditions, eps, type, 0);
			}
		}
		else if(thermalOrElectrical == ELECTRICFLAG)
		{
			if(type == PERIODIC)
			{
				fem(nodes, elements, physicals, parameters, solutionPotential, solutionCurrent, ELECTRICFLAG, PERIODICFLAG, conditions, eps, type, 0);
			}
			else if(type == DIRICHLET)
			{
				fem(nodes, elements, physicals, parameters, solutionPotential, solutionCurrent, ELECTRICFLAG, DIRICHLETFLAG, conditions, eps, type, 0);
			}
			else if(type == VONNEUMANN)
			{
				fem(nodes, elements, physicals, parameters, solutionPotential, solutionCurrent, ELECTRICFLAG, VONNEUMANNFLAG, conditions, eps, type, 0);
			}
		}
	}

	else if(FlagFE2 == 0 && nature == 3) //Coupled FE1 FE1
	{
		type = type_electric;
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
		cout << "\t##     FFFFFFFFFFFF        EEEEEEEEEEEE             ,/    ##" << endl;
		cout << "\t##     FF                  EE                     ,'/     ##" << endl;
		cout << "\t##     FF                  EE                   ,' /      ##" << endl;
		cout << "\t##     FF                  EE                 ,'  /_____, ##" << endl;
		cout << "\t##     FF                  EE               .'____    ,'  ##" << endl;
		cout << "\t##     FF                  EEEEEEEEEEEEEEEE      /  ,'    ##" << endl;
		cout << "\t##                                              / ,'      ##" << endl;
		cout << "\t###############################################/,' #########" << endl;
		cout << "\t##############################################/'############" << endl << endl << endl;

		if(type == PERIODIC)
		{
			fem(nodes, elements, physicals, parameters, solutionPotential, solutionCurrent, ELECTRICFLAG, PERIODICFLAG, conditions, eps, type, 0);
		}
		else if(type == DIRICHLET)
		{
			fem(nodes, elements, physicals, parameters, solutionPotential, solutionCurrent, ELECTRICFLAG, DIRICHLETFLAG, conditions, eps, type, 0);
		}
		else if(type == VONNEUMANN)
		{
			fem(nodes, elements, physicals, parameters, solutionPotential, solutionCurrent, ELECTRICFLAG, VONNEUMANNFLAG, conditions, eps, type, 0);
		}

		writeMSH((char*)"solutionPotential.pos", solutionPotential);
		writeMSH((char*)"solutionCurrent.pos", solutionCurrent);

		//Write in .dat file
		FILE *fp = fopen("dataMatlabV.dat", "w");
		std::map<Node*, std::vector<double> >::iterator itT = solutionPotential.begin();
		for(itT = solutionPotential.begin(); itT != solutionPotential.end(); itT++)
		fprintf(fp, "%.15f \t %.15f \t %.15f \n", itT->first->x, itT->first->y, itT->second[0]);

		fclose(fp);

		map<int, Parameter*> region;
		for(unsigned int i = 0; i < physicals.size(); i++)
		{
			for(unsigned int j = 0; j < parameters.size(); j++)
			{
				if(parameters[j]->name == physicals[i]->name)
				{
				    if(parameters[j]->dim != physicals[i]->dim)//Verification si erreurs entre les deux fichiers
				    {
				        cout << "Error: file.phy and file.msh do not correspond" << endl;
				    }
				    else
				    {
				        region[physicals[i]->num] = parameters[j];
				    }
				}
			}
		}
		Average_Joule_FE1(solutionPotential,elements,region);
		/*for(int e=0;e<elements.size();e++)
		{
			cout << elements[e]->s_Me;
		}*/
		cout<< endl;
		cout << "ELECTRICAL COMPUTATION FINISHED, STARTING THERMAL COMPUTATION..." << endl;
		cout<< endl<< endl;
		type = type_thermic;

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
		cout << "\t##     FF                  EE                     (       ##" << endl;
		cout << "\t##     FFFFFFFFFFFF        EEEEEEEEEEEE           ))      ##" << endl;
		cout << "\t##     FF                  EE                     {_}     ##" << endl;
		cout << "\t##     FF                  EE                    .-;-.    ##" << endl;
		cout << "\t##     FF                  EE                   |'-=-'|   ##" << endl;
		cout << "\t##     FF                  EE                   |     |   ##" << endl;
		cout << "\t##     FF                  EEEEEEEEEEEEEEEE     |     |   ##" << endl;
		cout << "\t##                                              |     |   ##" << endl;
		cout << "\t################################################|     |#####" << endl;
		cout << "\t################################################'.___.'#####" << endl << endl << endl;

		if(type == PERIODIC)
		{
			fem(nodes, elements, physicals, parameters, solutionTemperature, solutionFlux, THERMALFLAG, PERIODICFLAG, conditions, eps, type, 1);
		}
		else if(type == DIRICHLET)
		{
			fem(nodes, elements, physicals, parameters, solutionTemperature, solutionFlux, THERMALFLAG, DIRICHLETFLAG, conditions, eps, type, 1);
		}
		else if(type == VONNEUMANN)
		{
			fem(nodes, elements, physicals, parameters, solutionTemperature, solutionFlux, THERMALFLAG, VONNEUMANNFLAG, conditions, eps, type, 1);
		}

		writeMSH((char*)"solutionTemperature.pos", solutionTemperature);
		writeMSH((char*)"solutionFlux.pos", solutionFlux);

		//Write in .dat file
		fp = fopen("dataMatlabT.dat", "w");
		itT = solutionTemperature.begin();
		for(itT = solutionTemperature.begin(); itT != solutionTemperature.end(); itT++)
		fprintf(fp, "%.15f \t %.15f \t %.15f \n", itT->first->x, itT->first->y, itT->second[0]);

		fclose(fp);
	}

    //FE2 method.
    if(FlagFE2 ==1)
    {

		map<Node*, vector<double> > solutionScalars;
		map<Node*, vector<double> > solutionVectors;
		map<Node*, vector<double> > solutionScalars_micro;
		map<Node*, vector<double> > solutionVectors_micro;

        FE2(nodes_micro, elements_micro, physicals_micro, parameters_micro, solutionScalars_micro, solutionVectors_micro,
            conditions_micro, nodes, elements, physicals, parameters, solutionScalars, solutionVectors, conditions, eps,
			methodFE2, type_thermic, type_electric, argc, argv, nature);
    }

	if(FlagFE2 ==0 && nature !=3)
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

		//Write in .dat file
		FILE *fp = fopen("dataMatlabT.dat", "w");

		std::map<Node*, std::vector<double> >::iterator itT = solutionTemperature.begin();

		for(itT = solutionTemperature.begin(); itT != solutionTemperature.end(); itT++)
		{
		fprintf(fp, "%.15f \t %.15f \t %.15f \n", itT->first->x, itT->first->y, itT->second[0]);
		}

		fclose(fp);

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
	}
    return 0;
}



