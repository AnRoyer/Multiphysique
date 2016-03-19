#include <iostream>
#include <cstdio>
#include <vector>
#include <map>
#include <string>
#include "gmm/gmm.h"
#include "fem.h"
#include "NewtonRaphson.h"

using namespace std;

using namespace std;

/*---------------------------------------------------Code de calcul FEM-------------------------------------------------------*/


//type = 0 : thermique, type = 1 = électrique
void fem(std::vector<Node*> &nodes, std::vector<Element*> &elements, std::vector<Physical*> &physicals, std::vector<Parameter*> &parameters, std::map<Node*, std::vector<double> > &solutionTemperature, std::map<Node*, std::vector<double> > &solutionFlux, FemFlag type, FemFlag method, Periodique &conditions)
{
   // cout << "In FEM.cpp " << endl;
    //Boundaries
    map<int, Parameter*> region;//Stock le lien entre le numéro du physical de msh (stocker dans "physicals") et la valeur du parametre de "parametres" pour les régions de dimension 1 (ligne)

    /*Chargement des paramètres de la struc Parameter contenant les paramètres d'entrées.
    A chaque éléments de physicals, on lui associe l'élément de "parameters" correspondant. La correspondance est mappé dans linesRegion et surfaceRegion en fonction du type de paramètre
    */
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

    //Tri des noeuds
    Node *C1 = NULL,*C2 = NULL,*C3 = NULL,*C4 = NULL;//pointeurs vers les noeuds des coins
    vector<Node*> LeftNodes,RightNodes,TopNodes,BottomNodes;//vecteurs contenant les noeuds des bords
    map<Node*, Node*> NodesCorresp;//Vecteur de correspondance qui servira à additionner les lignes des noeuds en vis-a-vis

    double xmin = nodes[0]->x;
    double ymin = nodes[0]->y;
    double xmax = nodes[0]->x;
    double ymax = nodes[0]->y;
    for(unsigned int i = 0; i < nodes.size(); i++)
    {
        if(nodes[i]->x < xmin)
        {
            xmin = nodes[i]->x;
        }
        else if(nodes[i]->x > xmax)
        {
            xmax = nodes[i]->x;
        }
        else if(nodes[i]->y < ymin)
        {
            ymin = nodes[i]->y;
        }
        else if(nodes[i]->y > ymax)
        {
        ymax = nodes[i]->y;
        }
    }

    for(unsigned int i = 0; i < nodes.size(); i++)//Classement
    {
        NodesCorresp[nodes[i]] = nodes[i];//Initialisation du mapping en faisant pointer tous les noeuds vers eux-mêmes
        double x = nodes[i]->x;
        double y = nodes[i]->y;

        if(x == xmin && y != ymin && y != ymax)
        {
            LeftNodes.push_back(nodes[i]);
        }
        else if(x == xmax && y != ymin && y != ymax)
        {
            RightNodes.push_back(nodes[i]);
        }
        else if(y == ymax && x != xmin && x != xmax)
        {
            TopNodes.push_back(nodes[i]);
        }
        else if(y == ymin && x != xmin && x != xmax)
        {
            BottomNodes.push_back(nodes[i]);
        }
        else if(x == xmin && y == ymin)
        {
            C1 = nodes[i];
        }
        else if(x == xmax && y == ymin)
        {
            C2 = nodes[i];
        }
        else if(x == xmax && y == ymax)
        {
            C3 = nodes[i];
        }
        else if(x == xmin && y == ymax)
        {
            C4 = nodes[i];
        }
    }

    if(method == PERIODICFLAG)/*Si conditions periodiques alors le mapping "NodesCorresp" associe les noeuds en vis-a-vis. Sinon, le mapping est quand même utilisé mais tous les noeuds sont mappés vers eux-mêmes.*/
    {
        for(unsigned int i = 0; i < LeftNodes.size(); i++)
        {
            for(unsigned int j = 0; j < RightNodes.size(); j++)
            {
                if(abs(LeftNodes[i]->y - RightNodes[j]->y)<1e-5)
                {
                    NodesCorresp[RightNodes[j]] = LeftNodes[i];
                }
           }
        }

        for(unsigned int i = 0; i < BottomNodes.size(); i++)
        {
            for(unsigned int j = 0; j < TopNodes.size(); j++)
            {
                if(abs(BottomNodes[i]->x - TopNodes[j]->x)<1e-5)
                {
                    NodesCorresp[TopNodes[j]] = BottomNodes[i];
                }
            }
        }
    }





    //f vector
    vector<double> f(nodes.size());
    f_function(f, nodes, elements, region, type, 0); //dernier paramètre de la fonction f nul =>

    //Theta_K and delta_theta_k vector
    std::vector<double> theta_k(nodes.size(),0);
    std::vector<double> delta_theta_k(nodes.size());


    //Including the Dirichlet condition on theta_k
	std::vector<double> flag_theta(nodes.size()); //Vector aimed at keeping in mind the allocation
	double flag = 145;

    for(unsigned int l = 0; l < elements.size(); l++)
    {
        if(elements[l]->type == 1)//If line
        {
            if(region.count(elements[l]->region) == 1)//If linesRegion contains elements[i]->region
            {
                for(unsigned int j = 0; j < elements[l]->nodes.size(); j++)
                {
					if(flag_theta[elements[l]->nodes[j]->num-1]!= flag)
					{
					    if(region[elements[l]->region]->temperature != -1)
                        {
                            theta_k[elements[l]->nodes[j]->num-1] = region[elements[l]->region]->temperature;
                            flag_theta[elements[l]->nodes[j]->num-1]= flag;
                        }

					}
                }
            }
        }

    }

    std::vector<double> RHS(nodes.size());//Right Hand Side in Newton Raphson algorithm
    std::vector<double> q_m_x(nodes.size());
    std::vector<double> q_m_y(nodes.size());

    bool Criterion = false;
    int iter = 1; //counter of the number of iterations
    double normRHS0;//norm of initial RHS

    //Loop as long as criterion is not respected
    while(Criterion == false)
    {
        //Newton Raphson routine
        NewtonRaphson(nodes, elements, physicals, theta_k, f, method, NodesCorresp, delta_theta_k, region, RHS);
        //Initialize value for normRHS0
        if(iter==1)
            normRHS0 = gmm::vect_norm2(RHS);

        //Check the convergence criterion
        Criterion = End_Criterion(RHS,normRHS0);
        cout<<"Iteration number "<<iter<<endl;
        iter++;
    }

    Internal_flux(theta_k, region, elements, f, q_m_x, q_m_y);




    /*


     A UTILISER POUR L'INCLUSION DES CONDITIONS PERIODIQUES

    //Conditions periodiques sur K et f.
    double lx = abs(C2->x - C1->x);
    double ly = abs(C4->y - C1->y);
    double vol= lx * ly;//Volume du domaine

    if(method == PERIODICFLAG)
    {
        double Tavg = conditions.meanTemperature;
        double gradAvg_x = conditions.xGradient;
        double gradAvg_y = conditions.yGradient;

        //Condition correspondant au noeud 1: temperature moyenne.
        unsigned int numC1= C1->num-1;
        vector<double> c(nodes.size());//coefficients

        f_function(c,nodes,elements,surfaceRegion,1);

        vol = 0.0;
        for(unsigned int j=0; j<nodes.size(); j++)
        {
            Tmp(numC1,j) = c[j];
            vol += c[j];
        }

        f[numC1] = Tavg*vol;//condition

        //Conditions sur les noeuds 2 3 et 4.
        unsigned int numC2 = C2->num-1;
        unsigned int numC3 = C3->num-1;
        unsigned int numC4 = C4->num-1;

        for(unsigned int j=0; j<nodes.size(); j++)
        {
            Tmp(numC2,j) = 0;
            Tmp(numC3,j) = 0;
            Tmp(numC4,j) = 0;

            if(j == numC1)
            {
                Tmp(numC2,j) = -1;
                Tmp(numC3,j) = -1;
                Tmp(numC4,j) = -1;
            }

            if(j == numC2)
            {
                Tmp(numC2,j) = 1;
            }
            else if(j == numC3)
            {
                Tmp(numC3,j) = 1;
            }
            else if(j == numC4)
            {
                Tmp(numC4,j) = 1;
            }

        }

        //conditions sur f correspondants aux noeuds des coins
        f[numC2] = gradAvg_x * lx;
        f[numC3] = gradAvg_x * lx + gradAvg_y * ly;
        f[numC4] = gradAvg_y * ly;


        //conditions sur f correspondant aux noeuds à droite et en haut
        for(unsigned int i=0;i<RightNodes.size();i++)
        {
            unsigned int numNode = RightNodes[i]->num-1;
            for(unsigned int j=0; j<nodes.size(); j++)
            {
                Tmp(numNode,j) = 0;

                if(j == NodesCorresp[RightNodes[i]]->num-1)
                {
                    Tmp(numNode,j) = -1;
                }

                if(j == numNode)
                {
                    Tmp(numNode,j) = 1;
                }

            }

            f[numNode] = gradAvg_x * lx;
        }

        for(unsigned int i=0;i<TopNodes.size();i++)
        {
            unsigned int numNode = TopNodes[i]->num-1;

            for(unsigned int j=0; j<nodes.size(); j++)
            {
                Tmp(numNode,j) = 0;

                if(j == NodesCorresp[TopNodes[i]]->num-1)
                {
                    Tmp(numNode,j) = -1;
                }

                if(j == numNode)
                {
                    Tmp(numNode,j) = 1;
                }

            }
            f[numNode] = gradAvg_y * ly;
        }
    }
    */


    //Solution (écriture)
    for(unsigned int i = 0; i < nodes.size(); i++)
    {
        std::vector<double> val(1, theta_k[i]);
        solutionTemperature[nodes[i]] = val;
        std::vector<double> val2(3);
        val2[0] = q_m_x[i];
        val2[1] = q_m_y[i];
        val2[2] = 0;
        solutionFlux[nodes[i]] = val2;
    }
}//end of FEM routine



/* fonction permettant de calculer le vecteur f ainsi que la condition periodique sur le noeud 1 (condition sur la temperature moyenne)*/
void f_function(std::vector<double> &f, std::vector<Node*> &nodes, std::vector<Element*> &elements, std::map<int,Parameter*> &region, FemFlag type, int constantProperty)
{
    double cons;
    for(unsigned int i = 0; i < elements.size(); i++)
    {
        if(elements[i]->type == 2)//If triangle
        {
            if(constantProperty != 0)
            {
                cons = constantProperty;
            }
            else
            {
                if(type == THERMALFLAG)
                {
                    cons = region[elements[i]->region]->thermalGeneration;
                }
                else if(type == ELECTRICFLAG)
                {
                    cons = region[elements[i]->region]->electricalGeneration;
                }
                Node *n1 = elements[i]->nodes[0];
                Node *n2 = elements[i]->nodes[1];
                Node *n3 = elements[i]->nodes[2];
                double x1 = n1->x;
                double y1 = n1->y;
                double x2 = n2->x;
                double y2 = n2->y;
                double x3 = n3->x;
                double y3 = n3->y;
                double J = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);

                f[n1->num-1] += cons*J/6;
                f[n2->num-1] += cons*J/6;
                f[n3->num-1] += cons*J/6;
            }
        }
    }
}
