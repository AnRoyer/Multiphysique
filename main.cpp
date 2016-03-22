#include <iostream>
#include <cstdio>
#include <vector>
#include <map>
#include <string>
#include "gmshio.h"
#include "physicalio.h"
#include "fem.h"

using namespace std;

int main(int argc, char **argv)
{
    cout << "\t############################################################" << endl;
    cout << "\t############################################################" << endl;
    cout << "\t##                                                        ##" << endl;
    cout << "\t##                                             222222     ##" << endl;
    cout << "\t##                                                 22     ##" << endl;
    cout << "\t##     FFFFFFFFFFFFFFFF    EEEEEEEEEEEEEEEE    222222     ##" << endl;
    cout << "\t##     FF                  EE                  22         ##" << endl;
    cout << "\t##     FF                  EE                  222222     ##" << endl;
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

    if(argc < 3)
    {
        cout << "Usage: " << argv[0] << " file.msh file.phy" << endl;
        return 1;
    }

    vector<Node*> nodes;
    vector<Element*> elements;
    vector<Physical*> physicals;
    //lecture du MSH
    readMSH(argv[1], nodes, elements, physicals);

    //Affichage des infos principales
    cout << "Read " << nodes.size() << " nodes and " << elements.size() << " elements" << endl << endl;
    map<Node*, vector<double> > solutionTemperature;
    map<Node*, vector<double> > solutionFlux;


    vector<Parameter*> parameters;
    Periodique conditions;
    readPHY(argv[2], parameters, conditions);


    for(unsigned int i = 0; i < parameters.size(); i++)
    {
        if(parameters[i]->dim == 1 && conditions.exist == false)
        {
            cout << "Parameter " << parameters[i]->name << " has temperature of " << parameters[i]->temperature << " K and a electrical potential of " << parameters[i]->voltage << " V." << endl;
        }
        else if(parameters[i]->dim == 2)
        {
            cout << "Parameter " << parameters[i]->name << " has:" << endl;
            cout << "\t - A thermal conductivity of :" << endl;
            for(unsigned int j = 0; j < parameters[i]->thermalConductivity.size(); j++)
            {
                cout << "\t\t" << parameters[i]->thermalConductivity[j]->name << " = " << endl;
                cout << "\t\t|" << parameters[i]->thermalConductivity[j]->conductivity[0][0] << "\t" << parameters[i]->thermalConductivity[j]->conductivity[0][1] << "|" << endl;
                cout << "\t\t|" << parameters[i]->thermalConductivity[j]->conductivity[1][0] << "\t" << parameters[i]->thermalConductivity[j]->conductivity[1][1] << "|" << endl << endl;
            }
            cout << "\t - A heat production of " << parameters[i]->thermalGeneration << " K/s;" << endl;
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

    if(conditions.exist == true)
    {
        cout << "There is periodic conditions with:" << endl;
        cout << "\t - A mean temperature: " << conditions.meanTemperature << endl;
        cout << "\t - A mean x gradient: " << conditions.xGradient << endl;
        cout << "\t - A mean y gradient: " << conditions.yGradient << endl;
    }

    cout << endl << endl;

    if(conditions.exist == true)
    {
        fem(nodes, elements, physicals, parameters, solutionTemperature, solutionFlux, THERMALFLAG, PERIODICFLAG, conditions);
    }
    else
    {
        fem(nodes, elements, physicals, parameters, solutionTemperature, solutionFlux, THERMALFLAG, DIRICHLETFLAG, conditions);
    }

    writeMSH((char*)"solutionTemperature.pos", solutionTemperature);
    writeMSH((char*)"solutionFlux.pos", solutionFlux);

    FILE *fp = fopen("dataMatlabT.dat", "w");

    std::map<Node*, std::vector<double> >::iterator itT = solutionT.begin();

    for(itT = solutionT.begin(); itT != solutionT.end(); itT++)
    {
        fprintf(fp, "%.15f \t %.15f \t %.15f \n", itT->first->x, itT->first->y, itT->second[0]);
    }

    fclose(fp);

	/*

    fp = fopen("dataMatlabE.dat", "w");

    std::map<Node*, std::vector<double> >::iterator itE = solutionE.begin();

    for(itE = solutionE.begin(); itE != solutionE.end(); itE++)
    {
        fprintf(fp, "%.15f \t %.15f \t %.15f \n", itE->first->x, itE->first->y, itE->second[0]);
    }

    fclose(fp);*/

    return 0;
}
