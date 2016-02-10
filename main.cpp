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
    readPHY(argv[2], parameters);

    for(unsigned int i = 0; i < parameters.size(); i++)
    {
        if(parameters[i]->dim == 1)
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

    //FEM thermique
    fem(nodes, elements, physicals, parameters, solutionT, 0);
    //FEM électrique
    fem(nodes, elements, physicals, parameters, solutionE, 1);

    //Ecriture MSH thermique
    writeMSH((char*)"solT.pos", 0, 1, solutionT);
    //Ecriture MSH électrique
    writeMSH((char*)"solE.pos", 0, 1, solutionE);

    //fichier Matlab
    FILE *fp = fopen("dataMatlabT.dat", "w");

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

    fclose(fp);

    return 0;
}
