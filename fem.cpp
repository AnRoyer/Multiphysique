#include <iostream>
#include <cstdio>
#include <vector>
#include <map>
#include <string>
#include "gmm/gmm.h"
#include "fem.h"

using namespace std;

/*
Code de calcul FEM
*/

//type = 0 : thermique, type = 1 = électrique
void fem(std::vector<Node*> &nodes, std::vector<Element*> &elements, std::vector<Physical*> &physicals, std::vector<Parameter*> &parameters, std::map<Node*, std::vector<double> > &solution, int type)
{
    //Boundaries
    map<int, double> linesRegion;//Stock le lien entre le numéro du physical de msh (stocker dans "physicals") et la valeur du parametre de "parametres" pour les régions de dimension 1 (ligne)
    map<int, std::vector<double> > surfaceRegion;//Stock le lien entre le numéro du physical de msh (stocker dans "physicals") et la valeur du parametre de "parametres" pour les régions de dimension 2 (surface)

    /*Chargement des paramètres de la struc Parameter contenant les paramètres d'entrées.
    A chaque éléments de physicals, on lui associe l'élément de "parameters" correspondant. La correspondance est mappé dans linesRegion et surfaceRegion en fonction du type de paramètre
    */
    for(unsigned int i = 0; i < physicals.size(); i++)
    {
        //Si paramètre correspondant à une ligne (ex : condition de bord)
        if(physicals[i]->dim == 1)
        {
            for(unsigned int j = 0; j < parameters.size(); j++)
            {
                if(parameters[j]->name == physicals[i]->name)
                {
                    if(parameters[j]->dim != physicals[i]->dim)//Verification si erreurs entre les deux fichiers
                    {
                        cout << "Error: file.phy and file.msh do not correspond" << endl;
                    }

                    if(type == 0)
                    {
                        if(parameters[j]->value[0] != -1)
                        {
                            linesRegion[physicals[i]->num] = parameters[j]->value[0];
                        }
                    }
                    else if(type == 1)
                    {
                        if(parameters[j]->value[1] != -1)
                        {
                            linesRegion[physicals[i]->num] = parameters[j]->value[1];
                        }
                    }
                }
            }
        }
        else if(physicals[i]->dim == 2)
        {
            for(unsigned int j = 0; j < parameters.size(); j++)
            {
                if(parameters[j]->name == physicals[i]->name)
                {
                    if(parameters[j]->dim != physicals[i]->dim)
                    {
                        cout << "Error: file.phy and file.msh do not correspond" << endl;
                    }

                    if(type == 0)
                    {
                        std::vector<double> v(2);

                        v[0] = parameters[j]->value[0];
                        v[1] = parameters[j]->value[1];

                        surfaceRegion[physicals[i]->num] = v;
                    }
                    else if(type == 1)
                    {
                        std::vector<double> v(2);

                        v[0] = parameters[j]->value[2];
                        v[1] = parameters[j]->value[3];

                        surfaceRegion[physicals[i]->num] = v;
                    }
                }
            }
        }
    }

    //K matrix
    gmm::row_matrix< gmm::wsvector<double> > Tmp(nodes.size(), nodes.size());

    for(unsigned int i = 0; i < elements.size(); i++)
    {
        if(elements[i]->type == 2)//If triangle
        {
            Node *n1 = elements[i]->nodes[0];
            Node *n2 = elements[i]->nodes[1];
            Node *n3 = elements[i]->nodes[2];
            double x1 = n1->x;
            double y1 = n1->y;
            double x2 = n2->x;
            double y2 = n2->y;
            double x3 = n3->x;
            double y3 = n3->y;
            double cons = surfaceRegion[elements[i]->region][0];

            gmm::dense_matrix<double> Ke(3,3);
            double J = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);

            Ke(0,0) = ((y2-y3)*(y2-y3)+(x3-x2)*(x3-x2))/J;
            Ke(0,1) = (- (x2-x3)*(x1-x3) - (y2-y3)*(y1-y3))/J;
            Ke(0,2) = ((x2-x3)*(x1-x2) + (y2-y3)*(y1-y2))/J;
            Ke(1,0) = ((x1-x3)*(x3-x2) + (y1-y3)*(y3-y2))/J;
            Ke(1,1) = ((x1-x3)*(x1-x3) + (y1-y3)*(y1-y3))/J;
            Ke(1,2) = (- (x1-x3)*(x1-x2) - (y1-y3)*(y1-y2))/J;
            Ke(2,0) = ((x1-x2)*(x2-x3) + (y1-y2)*(y2-y3))/J;
            Ke(2,1) = (- (x1-x2)*(x1-x3) - (y1-y2)*(y1-y3))/J;
            Ke(2,2) = ((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2))/J;

            Tmp(n1->num-1, n1->num-1) += cons/2*Ke(0,0);
            Tmp(n1->num-1, n2->num-1) += cons/2*Ke(0,1);
            Tmp(n1->num-1, n3->num-1) += cons/2*Ke(0,2);
            Tmp(n2->num-1, n1->num-1) += cons/2*Ke(1,0);
            Tmp(n2->num-1, n2->num-1) += cons/2*Ke(1,1);
            Tmp(n2->num-1, n3->num-1) += cons/2*Ke(1,2);
            Tmp(n3->num-1, n1->num-1) += cons/2*Ke(2,0);
            Tmp(n3->num-1, n2->num-1) += cons/2*Ke(2,1);
            Tmp(n3->num-1, n3->num-1) += cons/2*Ke(2,2);
        }
    }
    //Dirichlet sur K
    for(unsigned int i = 0; i < elements.size(); i++)
    {
        if(elements[i]->type == 1)//If line
        {
            if(linesRegion.count(elements[i]->region) == 1)//If linesRegion contains elements[i]->region
            {
                for(unsigned int j = 0; j < elements[i]->nodes.size(); j++)
                {
                    for(unsigned int k = 0; k < nodes.size(); k++)
                    {
                        if(k == elements[i]->nodes[j]->num-1)
                        {
                            Tmp(elements[i]->nodes[j]->num-1, k) = 1;
                        }
                        else
                        {
                            Tmp(elements[i]->nodes[j]->num-1, k) = 0;
                        }
                    }
                }
            }
        }
    }

    //f vector
    vector<double> f(nodes.size());

    for(unsigned int i = 0; i < elements.size(); i++)
    {
        if(elements[i]->type == 2)//If triangle
        {
            Node *n1 = elements[i]->nodes[0];
            Node *n2 = elements[i]->nodes[1];
            Node *n3 = elements[i]->nodes[2];
            double x1 = n1->x;
            double y1 = n1->y;
            double x2 = n2->x;
            double y2 = n2->y;
            double x3 = n3->x;
            double y3 = n3->y;
            double cons = surfaceRegion[elements[i]->region][1];

            double J = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);

            f[n1->num-1] += cons*J/3;
            f[n2->num-1] += cons*J/6;
            f[n3->num-1] += cons*J/6;
        }
    }
    //Dirichlet sur f
    for(unsigned int i = 0; i < elements.size(); i++)
    {
        if(elements[i]->type == 1)//If line
        {
            if(linesRegion.count(elements[i]->region) == 1)//If linesRegion contains elements[i]->region
            {
                for(unsigned int j = 0; j < elements[i]->nodes.size(); j++)
                {
                    f[elements[i]->nodes[j]->num-1] = linesRegion[elements[i]->region];
                }
            }
        }
    }

    //System
    vector<double> x(nodes.size());
    gmm::csr_matrix<double> K;
    gmm::copy(Tmp,K);
    gmm::lu_solve(K, x, f);

    //Solution (écriture)
    for(unsigned int i = 0; i < nodes.size(); i++)
    {
        std::vector<double> val(1, x[i]);
        solution[nodes[i]] = val;
    }
}
