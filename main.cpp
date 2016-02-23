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
    map<Node*, vector<double> > solutionT;
    map<Node*, vector<double> > solutionE;


    vector<Parameter*> parameters;
    Periodique conditions;
    readPHY(argv[2], parameters, conditions);

    for(unsigned int i = 0; i < parameters.size(); i++)
    {
        if(parameters[i]->dim == 1 && conditions.exist == false)
        {
            cout << "Parameter " << parameters[i]->name << " has temperature of " << parameters[i]->value[0] << " K and a electrical potential of " << parameters[i]->value[1] << " V." << endl;
        }
        else if(parameters[i]->dim == 2)
        {
            cout << "Parameter " << parameters[i]->name << " has:" << endl;
            cout << "\t - A thermal diffusivity of " << parameters[i]->value[0] << " m^2/s with a heat production of " << parameters[i]->value[1] << " K/s;" << endl;
            cout << "\t - A electrical conductivity of " << parameters[i]->value[2] << " S/m with a change in charge density of " << parameters[i]->value[3] << " A/m^3." << endl;
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
        fem(nodes, elements, physicals, parameters, solutionT, THERMAL, PERIODIC, conditions);
    }
    else
    {
        fem(nodes, elements, physicals, parameters, solutionT, THERMAL, DIRICHLET, conditions);
    }

    writeMSH((char*)"solT.pos", 0, 1, solutionT);

    /*FILE *fp = fopen("dataMatlabT.dat", "w");

    std::map<Node*, std::vector<double> >::iterator itT = solutionT.begin();

    for(itT = solutionT.begin(); itT != solutionT.end(); itT++)
    {
        fprintf(fp, "%.15f \t %.15f \t %.15f \n", itT->first->x, itT->first->y, itT->second[0]);
    }

    fclose(fp);


    fp = fopen("dataMatlabE.dat", "w");

    std::map<Node*, std::vector<double> >::iterator itE = solutionE.begin();

    for(itE = solutionE.begin(); itE != solutionE.end(); itE++)
    {
        fprintf(fp, "%.15f \t %.15f \t %.15f \n", itE->first->x, itE->first->y, itE->second[0]);
    }

    fclose(fp);*/

    return 0;
}
