#include <iostream>
#include <cstring>
#include <cstdio>
#include <vector>
#include <string>
#include "physicalio.h"

using namespace std;


//fonction lisant le fichier PHY en tranférant les infos qu'il contient dans parameters
void readPHY(const char *fileName, std::vector<Parameter*> &parameters, std::vector<double> &cond_periodiques)
{
    FILE *fp = fopen(fileName, "r");
    if(!fp)//On verifie que le fichier soit bien ouvert
    {
        cout << "Error: cannot open file " << fileName << endl;
        return;
    }

    char string[256] = "";

    while(1)//Boucle de lecture
    {
        fscanf(fp, "%s", string);
        if(strcmp(string, "END") == 0)//On verifie si on est pas à la fin du fichier
        {
            break;
        }
        else if(strcmp(string, "PC") == 0)//On verifie si on est pas à la fin du fichier
        {
            fscanf(fp, "%lf", &cond_periodiques[0]);
            fscanf(fp, "%lf", &cond_periodiques[1]);
            fscanf(fp, "%lf", &cond_periodiques[2]);
        }
        else//Si pas "END" alors c'est un parametre
        {
            Parameter *p = new Parameter;
            p->name = string;//nom du parametre (= nom des physicals créer dans le .geo)

            if(fscanf(fp, "%d", &p->dim) != 1)//Si la lecture échoué
            {
                fclose(fp);
                return ;
            }

            double v;

            if(p->dim == 1)//Si c'est un parametre de dimension 1
            {
                for(unsigned int i = 0; i < 2; i++)
                {
                    if(fscanf(fp, "%lf", &v) != 1)
                    {
                        fclose(fp);
                        cout << "Error: corrupted file " << fileName << "." << endl;
                        return ;
                    }
                    p->value.push_back(v);//On stock la valeur des parametres
                }
            }
            else if(p->dim == 2)//Si c'est un parametre de dimension 2
            {
                for(unsigned int i = 0; i < 4; i++)
                {
                    if(fscanf(fp, "%lf", &v) != 1)
                    {
                        fclose(fp);
                        cout << "Error: corrupted file " << fileName << "." << endl;
                        return ;
                    }
                    p->value.push_back(v);//On stock la valeur des parametres
                }
            }

            parameters.push_back(p);//On stock le nouveau parametre
        }
    }

    fclose(fp);
}
