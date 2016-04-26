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

    if(argc < 5)
    {
        cout << "Usage: " << argv[0] << " macro.msh macro.phy micro.msh micro.phy" << endl;
        return 1;
    }

    vector<Node*> nodes_macro;
    vector<Node*> nodes_micro;
    vector<Element*> elements_macro;
    vector<Element*> elements_micro;
    vector<Physical*> physicals_macro;
    vector<Physical*> physicals_micro;
    //lecture des MSHs
    readMSH(argv[1], nodes_macro, elements_macro, physicals_macro);
    readMSH(argv[3], nodes_micro, elements_micro, physicals_micro);

    //Affichage des infos principales
    cout << "Read " << nodes_macro.size() << " nodes and " << elements_macro.size() << " elements for the macroscopic domain" << endl << endl;
    cout << "Read " << nodes_micro.size() << " nodes and " << elements_micro.size() << " elements for the microscopic domain" << endl << endl;
    map<Node*, vector<double> > solutionTemperature_macro;
    map<Node*, vector<double> > solutionFlux_macro;
    map<Node*, vector<double> > solutionTemperature_micro;
    map<Node*, vector<double> > solutionFlux_micro;


    vector<Parameter*> parameters_macro;
    vector<Parameter*> parameters_micro;
    Periodique conditions_macro;
    Periodique conditions_micro;
    readPHY(argv[2], parameters_macro, conditions_macro);
    readPHY(argv[4], parameters_micro, conditions_micro);


    cout << "For the macroscopic domain:" << endl;

    for(unsigned int i = 0; i < parameters_macro.size(); i++)
    {
        if(parameters_macro[i]->dim == 1 && conditions_macro.exist == false)
        {
            cout << "Parameter " << parameters_macro[i]->name << " has temperature of " << parameters_macro[i]->temperature << " K and a electrical potential of " << parameters_macro[i]->voltage << " V." << endl;
        }
        else if(parameters_macro[i]->dim == 2)
        {
            cout << "Parameter " << parameters_macro[i]->name << " has:" << endl;
            cout << "\t - A thermal conductivity of :" << endl;
            for(unsigned int j = 0; j < parameters_macro[i]->thermalConductivity.size(); j++)
            {
                cout << "\t\t" << parameters_macro[i]->thermalConductivity[j]->name << " = " << endl;
                cout << "\t\t|" << parameters_macro[i]->thermalConductivity[j]->conductivity[0][0] << "\t" << parameters_macro[i]->thermalConductivity[j]->conductivity[0][1] << "|" << endl;
                cout << "\t\t|" << parameters_macro[i]->thermalConductivity[j]->conductivity[1][0] << "\t" << parameters_macro[i]->thermalConductivity[j]->conductivity[1][1] << "|" << endl << endl;
            }
            cout << "\t - A heat production of " << parameters_macro[i]->thermalGeneration << " K/s;" << endl;
            cout << "\t - A electrical conductivity of :" << endl;
            for(unsigned int j = 0; j < parameters_macro[i]->electricalConductivity.size(); j++)
            {
                cout << "\t\t" << parameters_macro[i]->electricalConductivity[j]->name << " = " << endl;
                cout << "\t\t|" << parameters_macro[i]->electricalConductivity[j]->conductivity[0][0] << "\t" << parameters_macro[i]->electricalConductivity[j]->conductivity[0][1] << "|" << endl;
                cout << "\t\t|" << parameters_macro[i]->electricalConductivity[j]->conductivity[1][0] << "\t" << parameters_macro[i]->electricalConductivity[j]->conductivity[1][1] << "|" << endl << endl;
            }
            cout << "\t - A change in charge density of " << parameters_macro[i]->electricalGeneration << " A/m^3." << endl;
        }
    }

    if(conditions_macro.exist == true)
    {
        cout << "There is periodic conditions_macro with:" << endl;
        cout << "\t - A mean temperature: " << conditions_macro.meanTemperature << endl;
        cout << "\t - A mean x gradient: " << conditions_macro.xGradient << endl;
        cout << "\t - A mean y gradient: " << conditions_macro.yGradient << endl;
    }

    cout << endl << endl;

    cout << "For the microscopic domain:" << endl;

    for(unsigned int i = 0; i < parameters_micro.size(); i++)
    {
        if(parameters_micro[i]->dim == 1 && conditions_micro.exist == false)
        {
            cout << "Parameter " << parameters_micro[i]->name << " has temperature of " << parameters_micro[i]->temperature << " K and a electrical potential of " << parameters_micro[i]->voltage << " V." << endl;
        }
        else if(parameters_micro[i]->dim == 2)
        {
            cout << "Parameter " << parameters_micro[i]->name << " has:" << endl;
            cout << "\t - A thermal conductivity of :" << endl;
            for(unsigned int j = 0; j < parameters_micro[i]->thermalConductivity.size(); j++)
            {
                cout << "\t\t" << parameters_micro[i]->thermalConductivity[j]->name << " = " << endl;
                cout << "\t\t|" << parameters_micro[i]->thermalConductivity[j]->conductivity[0][0] << "\t" << parameters_micro[i]->thermalConductivity[j]->conductivity[0][1] << "|" << endl;
                cout << "\t\t|" << parameters_micro[i]->thermalConductivity[j]->conductivity[1][0] << "\t" << parameters_micro[i]->thermalConductivity[j]->conductivity[1][1] << "|" << endl << endl;
            }
            cout << "\t - A heat production of " << parameters_micro[i]->thermalGeneration << " K/s;" << endl;
            cout << "\t - A electrical conductivity of :" << endl;
            for(unsigned int j = 0; j < parameters_micro[i]->electricalConductivity.size(); j++)
            {
                cout << "\t\t" << parameters_micro[i]->electricalConductivity[j]->name << " = " << endl;
                cout << "\t\t|" << parameters_micro[i]->electricalConductivity[j]->conductivity[0][0] << "\t" << parameters_micro[i]->electricalConductivity[j]->conductivity[0][1] << "|" << endl;
                cout << "\t\t|" << parameters_micro[i]->electricalConductivity[j]->conductivity[1][0] << "\t" << parameters_micro[i]->electricalConductivity[j]->conductivity[1][1] << "|" << endl << endl;
            }
            cout << "\t - A change in charge density of " << parameters_micro[i]->electricalGeneration << " A/m^3." << endl;
        }
    }

    if(conditions_micro.exist == true)
    {
        cout << "There is periodic conditions_micro with:" << endl;
        cout << "\t - A mean temperature: " << conditions_micro.meanTemperature << endl;
        cout << "\t - A mean x gradient: " << conditions_micro.xGradient << endl;
        cout << "\t - A mean y gradient: " << conditions_micro.yGradient << endl;
    }

    cout << endl << endl;



    //Initial guess of the temperature field, reading macro.msh and macro.phy. The initial guess will correspond to solutionTemperature_macro.pos and will be run in DIRICHLET (see the macro.phy)
    double vol_macro;
    if(conditions_macro.exist == true)
    {
        fem(nodes_macro, elements_macro, physicals_macro, parameters_macro, solutionTemperature_macro, solutionFlux_macro, THERMALFLAG, PERIODICFLAG, conditions_macro);
    }
    else
    {
        fem(nodes_macro, elements_macro, physicals_macro, parameters_macro, solutionTemperature_macro, solutionFlux_macro, THERMALFLAG, DIRICHLETFLAG, conditions_macro);
    }

    //FE2 method.
    FE2(nodes_micro,elements_micro,physicals_micro,parameters_micro,solutionTemperature_micro,solutionFlux_micro,conditions_micro,nodes_macro, elements_macro,physicals_macro,parameters_macro,solutionTemperature_macro);

    writeMSH((char*)"solutionTemperature.pos", solutionTemperature_macro);
    writeMSH((char*)"solutionFlux.pos", solutionFlux_macro);

    FILE *fp = fopen("dataMatlabT.dat", "w");

    std::map<Node*, std::vector<double> >::iterator itT = solutionTemperature_macro.begin();

    for(itT = solutionTemperature_macro.begin(); itT != solutionTemperature_macro.end(); itT++)
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
