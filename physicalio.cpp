#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <string>
#include "stack.h"
#include "physicalio.h"
#include <mpi.h>

using namespace std;


//fonction lisant le fichier PHY en tranférant les infos qu'il contient dans parameters
void readPHY(const char *fileName, std::vector<Parameter*> &parameters, Periodique &conditions, Micro &micro,
             Type &typeUsed_thermic, Type &typeUsed_electric, int &nature, double &eps, std::vector<int> &methodFE2, int &natureFlag)
{
    ifstream fp(fileName);
    if(!fp.is_open())//On verifie que le fichier soit bien ouvert
    {
        cout << "Error: cannot open file " << fileName << endl;
        return;
    }

    string word;
    char c;
    Stack pile;
    Parameter* p;
    Type type_thermic = DEFAULTTYPE;
    Type type_electric = DEFAULTTYPE;
    Nature currentNature = DEFAULTNATURE;
    Dim currentDim = DEFAULTDIM;
    double epsRead = 0;
	int methodFE2Read_thermic;
	int methodFE2Read_electric;

    while(true)
    {
        fp.get(c);
        if(c == '<')
        {
            fp.get(c);
            if(c == '/')
            {
                word.clear();
                fp.get(c);
                while(c != '>')
                {
                    word += c;
                    fp.get(c);
                }

                if(word != pile.pop())
                {
                    cout << "Error: corrupted file .phy" << endl;
                }

                if(word == "phy")
                {
                    break;
                }
                else if(word == "line")
                {
                    parameters.push_back(p);
                    p = NULL;
                }
                else if(word == "surface")
                {
                    parameters.push_back(p);
                    p = NULL;
                }
            }
            else
            {
                word.clear();
                while(c != '>' && c != ' ')
                {
                    word += c;
                    fp.get(c);
                }

                pile.push(word);

                if(word == "phy")
                {
                    XMLparam param = readParam(fp);

                    while(param.name != "NULL")
                    {
                        if(param.name == "type_thermic")
                        {
                            if(param.value == "dirichlet")
                            {
                                //cout << "Reading a dirichlet phy file ..." << endl;
                                type_thermic = DIRICHLET;
                                conditions.exist = false;
                            }
                            else if(param.value == "periodic")
                            {
                                //cout << "Reading a periodic phy file ..." << endl;
                                type_thermic = PERIODIC;
                                conditions.exist = true;
                            }
                            else if(param.value == "vonNeumann")
                            {
                                //cout << "Reading a von Neumann phy file ..." << endl;
                                type_thermic = VONNEUMANN;
                                conditions.exist = false;
                            }
                            else if(param.value == "fe2D")
                            {
                                //cout << "Reading a FE2 with dirichlet phy file ..." << endl;
                                type_thermic = FE2withDIRICHLET;
                                conditions.exist = false;
                            }
                            else if(param.value == "fe2V")
                            {
                                //cout << "Reading a FE2 with von Neumann phy file ..." << endl;
                                type_thermic = FE2withVONNEUMANN;
                                conditions.exist = false;
                            }
                            else if(param.value == "fe2P")
                            {
                                //cout << "Reading a FE2 with periodic phy file ..." << endl;
                                type_thermic = FE2withPERIODIC;
                                conditions.exist = true;
                            }
                            else
                            {
                                cout << "Error: unknown type" << endl;
                            }
						}
						else if(param.name == "type_electric")
                        {
                            if(param.value == "dirichlet")
                            {
                                //cout << "Reading a dirichlet phy file ..." << endl;
                                type_electric = DIRICHLET;
                                conditions.exist = false;
                            }
                            else if(param.value == "periodic")
                            {
                                //cout << "Reading a periodic phy file ..." << endl;
                                type_electric = PERIODIC;
                                conditions.exist = true;
                            }
                            else if(param.value == "vonNeumann")
                            {
                                //cout << "Reading a von Neumann phy file ..." << endl;
                                type_electric = VONNEUMANN;
                                conditions.exist = false;
                            }
                            else if(param.value == "fe2D")
                            {
                                //cout << "Reading a FE2 with dirichlet phy file ..." << endl;
                                type_electric = FE2withDIRICHLET;
                                conditions.exist = false;
                            }
                            else if(param.value == "fe2V")
                            {
                                //cout << "Reading a FE2 with von Neumann phy file ..." << endl;
                                type_electric = FE2withVONNEUMANN;
                                conditions.exist = false;
                            }
                            else if(param.value == "fe2P")
                            {
                                //cout << "Reading a FE2 with periodic phy file ..." << endl;
                                type_electric = FE2withPERIODIC;
                                conditions.exist = true;
                            }
                            else
                            {
                                cout << "Error: unknown type" << endl;
                            }
						}
                        else if(param.name == "microMSH")
                        {
                            micro.fileMsh = param.value;
                        }
                        else if(param.name == "microPHY")
                        {
                            micro.filePhy = param.value;
                        }
                        else if(param.name == "res")
                        {
                            sscanf(param.value.c_str(), "%lf", &epsRead);
                        }
                        else if(param.name == "methodFE2_thermic")
                        {
                            sscanf(param.value.c_str(), "%d", &methodFE2Read_thermic);
                        }
                        else if(param.name == "methodFE2_electric")
                        {
                            sscanf(param.value.c_str(), "%d", &methodFE2Read_electric);
                        }
                        else if(param.name == "nature")
                        {
                            if(param.value == "thermic") nature = 1;
                            if(param.value == "electric") nature = 2;
                            if(param.value == "coupled") nature = 3;
                        }
                        else
                        {
                            cout << "Error: unknown parameter " << param.name << " for markup <" << pile.peek() << ">" << endl;
                        }

                        param = readParam(fp);
                    }

                    if(type_electric == DEFAULTTYPE && type_thermic == DEFAULTTYPE)
                    {
                        cout << "Error: no type specified" << endl;
                        return;
                    }
                }
                else if(word == "line")
                {
                    currentDim = LINE;

                    XMLparam param = readParam(fp);

                    while(param.name != "NULL")
                    {
                        if(param.name == "name")
                        {
                            p = new Parameter;
                            p->name = param.value;
                            p->dim = 1;
                        }
                        else
                        {
                            cout << "Error: unknown parameter " << param.name << " for markup <" << pile.peek() << ">" << endl;
                        }

                        param = readParam(fp);
                    }

                    if(p == NULL)
                    {
                        cout << "Error: no name specified for line" << endl;
                        return;
                    }
                }
                else if(word == "surface")
                {
                    currentDim = SURFACE;

                    XMLparam param = readParam(fp);

                    while(param.name != "NULL")
                    {
                        if(param.name == "name")
                        {
                            p = new Parameter;
                            p->name = param.value;
                            p->dim = 2;
                        }
                        else
                        {
                            cout << "Error: unknown parameter " << param.name << " for markup <" << pile.peek() << ">" << endl;
                        }

                        param = readParam(fp);
                    }

                    if(p == NULL)
                    {
                        cout << "Error: no name specified for surface" << endl;
                        return;
                    }
                }
                else if(word == "global")
                {
                    currentDim = GLOBAL;
                }
                else if(word == "thermal")
                {
                    currentNature = THERMAL;
                    natureFlag = (natureFlag|THERMALDATA);
                }
                else if(word == "electrical")
                {
                    currentNature = ELECTRICAL;
                    natureFlag = (natureFlag|ELECTRICALDATA);
                }
                else if(word == "value")
                {
                    if(currentDim == LINE)
                    {
                        XMLparam param = readParam(fp);

                        int typeValue = 0;
                        std::string name;

                        while(param.name != "NULL")
                        {
                            if(param.value == "reference")
                            {
                                typeValue = 1;
                            }
                            else if(param.value == "flux")
                            {
                                typeValue = 2;
                            }

                            param = readParam(fp);
                        }

                        if(typeValue == 0)
                        {
                            cout << "Error: no type value specified" << endl;
                        }
                        else if(typeValue == 1)
                        {
                            if(currentNature == THERMAL)
                            {
                                p->temperature = readValue(fp);
                            }
                            else if(currentNature == ELECTRICAL)
                            {
                                p->voltage = readValue(fp);
                            }
                            else
                            {
                               cout << "Error: unknown current nature" << endl;
                            }
                        }
                        else if(typeValue == 2)
                        {
                            if(currentNature == THERMAL)
                            {
                                p->fluxTemperature = readValue(fp);
                            }
                            else if(currentNature == ELECTRICAL)
                            {
                                p->fluxVoltage = readValue(fp);
                            }
                            else
                            {
                               cout << "Error: unknown current nature" << endl;
                            }
                        }
                    }
                    else if(currentDim == SURFACE)
                    {
                        XMLparam param = readParam(fp);

                        int typeValue = 0;
                        std::string name;

                        while(param.name != "NULL")
                        {
                            if(param.name == "type")
                            {
                                if(param.value == "conductivity")
                                {
                                    typeValue = 1;
                                }
                                else if(param.value == "generation")
                                {
                                    typeValue = 2;
                                }
                            }
                            else if(param.name == "name")
                            {
                                name = param.value;
                            }

                            param = readParam(fp);
                        }

                        if(typeValue == 0)
                        {
                            cout << "Error: no type value specified" << endl;
                        }
                        else if(typeValue == 1)
                        {
                            Conductivity* newCond = new Conductivity;
                            newCond->name = name;

                            newCond->conductivity[0][0] = readValue(fp);
                            newCond->conductivity[0][1] = readValue(fp);
                            newCond->conductivity[1][0] = readValue(fp);
                            newCond->conductivity[1][1] = readValue(fp);

                            if(currentNature == THERMAL)
                            {
                                p->thermalConductivity.push_back(newCond);
                            }
                            else if(currentNature == ELECTRICAL)
                            {
                                p->electricalConductivity.push_back(newCond);
                            }
                            else
                            {
                               cout << "Error: unknown current nature" << endl;
                            }
                        }
                        else if(typeValue == 2)
                        {
                            if(currentNature == THERMAL)
                            {
                                p->thermalGeneration = readValue(fp);
                            }
                            else if(currentNature == ELECTRICAL)
                            {
                                p->electricalGeneration = readValue(fp);
                            }
                        }
                    }
                }
                else if(word == "mean")
                {
                    if(currentNature == THERMAL)
                    {
                        conditions.meanTemperature = readValue(fp);
                    }
                    else if(currentNature == ELECTRICAL)
                    {
                        conditions.meanVoltage = readValue(fp);
                    }
                    else
                    {
                        cout << "Error: unknown current nature" << endl;
                    }
                }
                else if(word == "xgradient")
                {
                    if(currentNature == THERMAL)
                    {
                        conditions.xGradient = readValue(fp);
                    }
                    else if(currentNature == ELECTRICAL)
                    {
                        conditions.xGradient = readValue(fp);
                    }
                    else
                    {
                        cout << "Error: unknown current nature" << endl;
                    }
                }
                else if(word == "ygradient")
                {
                    if(currentNature == THERMAL)
                    {
                        conditions.yGradient = readValue(fp);
                    }
                    else if(currentNature == ELECTRICAL)
                    {
                        conditions.yGradient = readValue(fp);
                    }
                    else
                    {
                        cout << "Error: unknown current nature" << endl;
                    }
                }
            }
        }
    }

	typeUsed_thermic = type_thermic;
	typeUsed_electric = type_electric;

    if(epsRead == 0)
    {
        eps = 1e-5;
    }
    else
    {
        eps = epsRead;
    }

    if(methodFE2Read_thermic == 0)
    {
        methodFE2[0] = 2;
    }
    else
    {
        methodFE2[0] = methodFE2Read_thermic;
    }

    if(methodFE2Read_electric == 0)
    {
        methodFE2[1] = 2;
    }
    else
    {
        methodFE2[1] = methodFE2Read_electric;
    }

    //cout << "End of the file reaches" << endl;

    fp.close();
}

XMLparam readParam(ifstream& fp)
{
    XMLparam exit;

    string word;
    char c;

    while(c != '=')
    {
        fp.get(c);
        if(c != ' ' && c != '=' && c != '>')
        {
            word += c;
        }
        else if(c == '>')
        {
            exit.name = "NULL";
            return exit;
        }
    }

    exit.name = word;

    word.clear();
    while(true)
    {
        fp.get(c);

        if(c == '"')
        {
            while(true)
            {
                fp.get(c);
                if(c == '"')
                {
                    break;
                }
                else
                {
                    word += c;
                }
            }

            exit.value = word;
            break;
        }
    }

    return exit;
}

double readValue(ifstream& fp)
{
    string word;
    char c;
    double value = -1;

    while(true)
    {
        fp.get(c);
        if(c == '$')
        {
            fp.get(c);
            while(c != '$')
            {
                word += c;
                fp.get(c);
            }

            if(word == "NULL")
            {
                value = -1;
            }
            else
            {
                value = atof(word.c_str());
            }

            break;
        }
        else if(c == '<')
        {
            cout << "Error: missing a parameter" << endl;
        }
    }

    return value;
}
