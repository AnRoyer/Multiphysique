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
void fem(std::vector<Node*> &nodes, std::vector<Element*> &elements, std::vector<Physical*> &physicals, std::vector<Parameter*> &parameters, std::map<Node*, std::vector<double> > &solution, FemFlag type, FemFlag method, Periodique &conditions)
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
        if(physicals[i]->dim == 1 && method == DIRICHLET)
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

    //Tri des noeuds
    Node *C1,*C2,*C3,*C4;//pointeurs vers les noeuds des coins
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
        int num = nodes[i]->num;
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

    if(method == PERIODIC)/*Si conditions periodiques alors le mapping "NodesCorresp" associe les noeuds en vis-a-vis. Sinon, le mapping est quand même utilisé mais tous les noeuds sont mappés vers eux-mêmes.*/
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
            double cons = surfaceRegion[elements[i]->region][0];// alpha

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


            /*Utilisation du mapping NodesCorresp.
            Si une valeur doit être ajoutée à la ligne d'un noeud de droite, celle-ci est directement ajoutée à la ligne correspondant au noeud de gauche en vis-a-vis*/
            int num1 = NodesCorresp[n1]->num-1;
            int num2 = NodesCorresp[n2]->num-1;
            int num3 = NodesCorresp[n3]->num-1;

            Tmp(num1, n1->num-1) += cons/2*Ke(0,0);
            Tmp(num1, n2->num-1) += cons/2*Ke(0,1);
            Tmp(num1, n3->num-1) += cons/2*Ke(0,2);
            Tmp(num2, n1->num-1) += cons/2*Ke(1,0);
            Tmp(num2, n2->num-1) += cons/2*Ke(1,1);
            Tmp(num2, n3->num-1) += cons/2*Ke(1,2);
            Tmp(num3, n1->num-1) += cons/2*Ke(2,0);
            Tmp(num3, n2->num-1) += cons/2*Ke(2,1);
            Tmp(num3, n3->num-1) += cons/2*Ke(2,2);
        }
    }

    //f vector
    vector<double> f(nodes.size());
    f_function(f,nodes,elements,surfaceRegion,0); //dernier paramètre de la fonction f nul =>

    //Dirichlet sur K
    if(method == DIRICHLET)
    {
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
    }

    //Conditions periodiques sur K et f.
    double lx = abs(C2->x - C1->x);
    double ly = abs(C4->y - C1->y);
    double vol= lx * ly;//Volume du domaine
    if(method == PERIODIC)
    {
        double Tavg = conditions.meanTemperature;
        double gradAvg_x = conditions.xGradient;
        double gradAvg_y = conditions.yGradient;

        //Condition correspondant au noeud 1: temperature moyenne.
        int numC1= C1->num-1;
        vector<double> c(nodes.size());//coefficients

        f_function(c,nodes,elements,surfaceRegion,1);

        vol = 0.0;
        for(int j=0; j<nodes.size(); j++)
        {
            Tmp(numC1,j) = c[j];
            vol += c[j];
        }

        f[numC1] = Tavg*vol;//condition

        //Conditions sur les noeuds 2 3 et 4.
        int numC2 = C2->num-1;
        int numC3 = C3->num-1;
        int numC4 = C4->num-1;

        for(int j=0; j<nodes.size(); j++)
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
        for(int i=0;i<RightNodes.size();i++)
        {
            int numNode = RightNodes[i]->num-1;
            for(int j=0; j<nodes.size(); j++)
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

        for(int i=0;i<TopNodes.size();i++)
        {
            int numNode = TopNodes[i]->num-1;

            for(int j=0; j<nodes.size(); j++)
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

    //System
    vector<double> x(nodes.size());
    gmm::csr_matrix<double> K;
    gmm::copy(Tmp,K);
    gmm::lu_solve(K, x, f);

    //gradiant moyen x
    double fluxx = 0.0;
    vector<double> g(nodes.size());
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

            g[n1->num-1] += cons*(1./2)*(y2-y3);
            g[n2->num-1] += cons*(1./2)*(y3-y1);
            g[n3->num-1] += cons*(1./2)*(y1-y2);
        }
    }
    for(unsigned int i = 0; i < nodes.size(); i++)
    {
        fluxx += g[i] * x[i];
    }
    fluxx = fluxx/vol;

    //gradiant moyen y
    double fluxy = 0.0;
    vector<double> gy(nodes.size());
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

            gy[n1->num-1] += cons*(1./2)*(x3-x2);
            gy[n2->num-1] += cons*(1./2)*(x1-x3);
            gy[n3->num-1] += cons*(1./2)*(x2-x1);
        }
    }
    for(unsigned int i = 0; i < nodes.size(); i++)
    {
        fluxy += gy[i] * x[i];
    }
    fluxy = fluxy/vol;

    cout << "Average heat flux along x: q_x = " << fluxx << " [m K/s]" << endl;
    cout << "Average heat flux along y: q_y = " << fluxy << " [m K/s]" << endl;


    //Solution (écriture)
    for(unsigned int i = 0; i < nodes.size(); i++)
    {
        std::vector<double> val(1, x[i]);
        solution[nodes[i]] = val;
    }
}

/* fonction permettant de calculer le vecteur f ainsi que la condition periodique sur le noeud 1 (condition sur la temperature moyenne)*/
void f_function(std::vector<double> &f, std::vector<Node*> &nodes, std::vector<Element*> &elements, std::map<int,std::vector<double> > &surfaceRegion, int constantProperty)
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
                cons = surfaceRegion[elements[i]->region][1];// Dans le cas du calcul du vecteur f.
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
