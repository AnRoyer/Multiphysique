#include <iostream>
#include <cstdio>
#include <vector>
#include <map>
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
        cout << "Usage: " << argv[0] << " file.msh file.phy" << endl;
        return 1;
    }

    Type type;
    vector<Parameter*> parameters;
    Periodique conditions;
    Micro micro;
    double eps = 0;
	int methodFE2 = 0;

    Type type_micro;
    vector<Parameter*> parameters_micro;
    Periodique conditions_micro;
    Micro micro_micro;
    double eps_micro = 0;
	int methodFE2_micro = 0; //no signification !!

    //lecture des PHYs
    readPHY(argv[2], parameters, conditions, micro, type, eps, methodFE2);

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
        readPHY(micro.filePhy.c_str(), parameters_micro, conditions_micro, micro_micro, type_micro, eps_micro, methodFE2_micro);
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
		cout << "Read " << nodes.size() << " nodes and " << elements.size() << " elements." << endl << endl;
		if(type == DIRICHLET) cout << "Calling the FE1 method with DIRICHLET conditions." << endl;
		if(type == VONNEUMANN) cout << "Calling the FE1 method with VONNEUMANN conditions." << endl;
		if(type == PERIODIC) cout << "Calling the FE1 method with PERIODIC conditions." << endl;
		cout << endl;
	}

    map<Node*, vector<double> > solutionTemperature;
    map<Node*, vector<double> > solutionFlux;
    map<Node*, vector<double> > solutionTemperature_micro;
    map<Node*, vector<double> > solutionFlux_micro;

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

    //FE2 method.
    if(type == FE2withDIRICHLET || type == FE2withVONNEUMANN || type == FE2withPERIODIC)
    {
        FE2(nodes_micro, elements_micro, physicals_micro, parameters_micro, solutionTemperature_micro, solutionFlux_micro,
            conditions_micro, nodes, elements, physicals, parameters, solutionTemperature, solutionFlux, eps, argc, argv,
			methodFE2, THERMALFLAG, type);
    }

	if(type == DIRICHLET || type == VONNEUMANN || type == PERIODIC)
	{
		writeMSH((char*)"solutionTemperature.pos", solutionTemperature);
		writeMSH((char*)"solutionFlux.pos", solutionFlux);
		cout << endl;
		if(type == DIRICHLET) cout << "The problem has been solved in one scale using DIRICHLET conditions." << endl;
		if(type == VONNEUMANN) cout << "The problem has been solved in one scale using VONNEUMANN conditions." << endl;
		if(type == PERIODIC) cout << "The problem has been solved in one scale using PERIODIC conditions." << endl;
		cout << endl;
	}

    return 0;
}
