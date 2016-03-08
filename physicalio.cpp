#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <string>
#include "stack.h"
#include "physicalio.h"

using namespace std;


//fonction lisant le fichier PHY en tranférant les infos qu'il contient dans parameters
void readPHY(const char *fileName, std::vector<Parameter*> &parameters, Periodique &conditions)
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
    Type type = DEFAULTTYPE;
    Nature currentNature = DEFAULTNATURE;
    Dim currentDim = DEFAULTDIM;

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
                    if(c == ' ')
                    {
                        XMLparam param = readParam(fp);
                        if(param.name == "NULL")
                        {
                            cout << "Error: no type specified" << endl;
                        }

                        if(param.name == "type")
                        {
                            if(param.value == "dirichlet")
                            {
                                cout << "Reading a dirichlet phy file ..." << endl;
                                type = DIRICHLET;
                                conditions.exist = false;
                            }
                            else if(param.value == "periodic")
                            {
                                cout << "Reading a periodic phy file ..." << endl;
                                type = PERIODIC;
                                conditions.exist = true;
                            }
                            else
                            {
                                cout << "Error: unknown type" << endl;
                            }
                        }
                        else
                        {
                            cout << "Error: unknown parameter " << param.name << " for markup <" << pile.peek() << ">" << endl;
                        }
                    }
                    else if(c == '>')
                    {
                        cout << "Error: no type specified" << endl;
                    }
                }
                else if(word == "line")
                {
                    currentDim = LINE;
                    if(c == ' ')
                    {
                        XMLparam param = readParam(fp);
                        if(param.name == "NULL")
                        {
                            cout << "Error: no name specified" << endl;
                        }

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
                    }
                    else if(c == '>')
                    {
                        cout << "Error: no name specified" << endl;
                    }
                }
                else if(word == "surface")
                {
                    currentDim = SURFACE;
                    if(c == ' ')
                    {
                        XMLparam param = readParam(fp);
                        if(param.name == "NULL")
                        {
                            cout << "Error: no name specified" << endl;
                        }

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
                    }
                    else if(c == '>')
                    {
                        cout << "Error: no name specified" << endl;
                    }
                }
                else if(word == "global")
                {
                    currentDim = GLOBAL;
                }
                else if(word == "thermal")
                {
                    currentNature = THERMAL;
                }
                else if(word == "electrical")
                {
                    currentNature = ELECTRICAL;
                }
                else if(word == "value" && currentDim == LINE)
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
                else if(word == "value" && currentDim == SURFACE)
                {
                    XMLparam param1 = readParam(fp);
                    XMLparam param2 = readParam(fp);

                    XMLparam param;
                    XMLparam nameParam;

                    if(param1.name == "type")
                    {
                        param = param1;
                    }
                    else if(param2.name == "type")
                    {
                        param = param2;
                    }

                    if(param1.name == "name")
                    {
                        nameParam = param1;
                    }
                    else if(param2.name == "name")
                    {
                        nameParam = param2;
                    }

                    if(param.name == "type")
                    {
                        if(param.value == "conductivity")
                        {
                            if(currentNature == THERMAL)
                            {
                                if(nameParam.name == "NULL")
                                {
                                    cout << "Error: no name specified" << endl;
                                }

                                Conductivity* newCond = new Conductivity;
                                newCond->name = nameParam.value;

                                cout << nameParam.value << endl;

                                newCond->conductivity[0][0] = readValue(fp);
                                newCond->conductivity[0][1] = readValue(fp);
                                newCond->conductivity[1][0] = readValue(fp);
                                newCond->conductivity[1][1] = readValue(fp);

                                p->thermalConductivity.push_back(newCond);
                            }
                            else if(currentNature == ELECTRICAL)
                            {
                                if(nameParam.name == "NULL")
                                {
                                    cout << "Error: no name specified" << endl;
                                }

                                Conductivity* newCond = new Conductivity;
                                newCond->name = nameParam.value;

                                newCond->conductivity[0][0] = readValue(fp);
                                newCond->conductivity[0][1] = readValue(fp);
                                newCond->conductivity[1][0] = readValue(fp);
                                newCond->conductivity[1][1] = readValue(fp);

                                p->electricalConductivity.push_back(newCond);
                            }
                            else
                            {
                               cout << "Error: unknown current nature" << endl;
                            }
                        }
                        else if(param.value == "generation")
                        {
                            if(currentNature == THERMAL)
                            {
                                p->thermalGeneration = readValue(fp);
                            }
                            else if(currentNature == ELECTRICAL)
                            {
                                p->electricalGeneration = readValue(fp);
                            }
                            else
                            {
                               cout << "Error: unknown current nature" << endl;
                            }
                        }
                        else
                        {
                            cout << "Error: unknown value " << param.value << " in parameter " << param.name << " for markup <" << pile.peek() << ">" << endl;
                        }
                    }
                    else
                    {
                        cout << "Error: unknown parameter " << param.name << " for markup <" << pile.peek() << ">" << endl;
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
                        conditions.meanTemperature = readValue(fp);
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

    cout << "End of the file reaches" << endl;

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
    double value = 0;

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
