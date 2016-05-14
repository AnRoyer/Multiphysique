#include <iostream>
#include <cstdio>
#include <vector>
#include <map>
#include <string>
#include "gmm/gmm.h"
#include "fem.h"
#include "NewtonRaphson.h"
#include <mpi.h>

using namespace std;

/*---------------------------------------------------Code de calcul FEM-------------------------------------------------------*/


void fem(std::vector<Node*> &nodes, std::vector<Element*> &elements, std::vector<Physical*> &physicals, std::vector<Parameter*> &parameters,
         std::map<Node*, std::vector<double> > &solutionTemperature, std::map<Node*, std::vector<double> > &solutionFlux,
		 FemFlag thermalOrElectrical,
         FemFlag method, Periodique &conditions, double eps, Type type, int JouleEffect)
{
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
    NodeCorner corner;
    NodeBorder border;//vecteurs contenant les noeuds des bords
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
            border.LeftNodes.push_back(nodes[i]);
        }
        else if(x == xmax && y != ymin && y != ymax)
        {
            border.RightNodes.push_back(nodes[i]);
        }
        else if(y == ymax && x != xmin && x != xmax)
        {
            border.TopNodes.push_back(nodes[i]);
        }
        else if(y == ymin && x != xmin && x != xmax)
        {
            border.BottomNodes.push_back(nodes[i]);
        }
        else if(x == xmin && y == ymin)
        {
            corner.C1 = nodes[i];
        }
        else if(x == xmax && y == ymin)
        {
            corner.C2 = nodes[i];
        }
        else if(x == xmax && y == ymax)
        {
            corner.C3 = nodes[i];
        }
        else if(x == xmin && y == ymax)
        {
            corner.C4 = nodes[i];
        }
    }

    if(method == PERIODICFLAG)/*Si conditions periodiques alors le mapping "NodesCorresp" associe les noeuds en vis-a-vis. Sinon, le mapping est quand même utilisé mais tous les noeuds sont mappés vers eux-mêmes.*/
    {
        for(unsigned int i = 0; i < border.LeftNodes.size(); i++)
        {
            for(unsigned int j = 0; j < border.RightNodes.size(); j++)
            {
                if(abs(border.LeftNodes[i]->y - border.RightNodes[j]->y) < 1e-5)
                {
                    NodesCorresp[border.RightNodes[j]] = border.LeftNodes[i];
                }
           }
        }

        for(unsigned int i = 0; i < border.BottomNodes.size(); i++)
        {
            for(unsigned int j = 0; j < border.TopNodes.size(); j++)
            {
                if(abs(border.BottomNodes[i]->x - border.TopNodes[j]->x) < 1e-5)
                {
                    NodesCorresp[border.TopNodes[j]] = border.BottomNodes[i];
                }
            }
        }
    }

    //f vector
    vector<double> f(nodes.size());
    f_function(f, nodes, elements, region, thermalOrElectrical, 0, JouleEffect); //dernier paramètre de la fonction f nul =>
    
	//Theta_K and delta_theta_k vector
    std::vector<double> theta_k(nodes.size(),1);
    std::vector<double> delta_theta_k(nodes.size());

    if(method == VONNEUMANNFLAG)
    {
        for(unsigned int l = 0; l < elements.size(); l++)
        {
            if(elements[l]->type == 1)//If line
            {
                if(region.count(elements[l]->region) == 1)
                {
                    if(thermalOrElectrical == THERMALFLAG)
                    {
                        if(region[elements[l]->region]->fluxTemperature != -1)
                        {
                            Node *n1 = elements[l]->nodes[0];
                            Node *n2 = elements[l]->nodes[1];

                            double x1 = n1->x;
                            double y1 = n1->y;
                            double x2 = n2->x;
                            double y2 = n2->y;

                            double longueur = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));

                            f[elements[l]->nodes[0]->num-1] -= longueur*region[elements[l]->region]->fluxTemperature/2;
                            f[elements[l]->nodes[1]->num-1] -= longueur*region[elements[l]->region]->fluxTemperature/2;
                        }

                        if(region[elements[l]->region]->temperature != -1)
                        {
                            theta_k[elements[l]->nodes[0]->num-1] = region[elements[l]->region]->temperature;
                            theta_k[elements[l]->nodes[1]->num-1] = region[elements[l]->region]->temperature;
                        }
                    }
                    else if(thermalOrElectrical == ELECTRICFLAG)
                    {
                        if(region[elements[l]->region]->fluxVoltage != -1)
                        {
                            Node *n1 = elements[l]->nodes[0];
                            Node *n2 = elements[l]->nodes[1];

                            double x1 = n1->x;
                            double y1 = n1->y;
                            double x2 = n2->x;
                            double y2 = n2->y;

                            double longueur = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));

                            f[elements[l]->nodes[0]->num-1] -= longueur*region[elements[l]->region]->fluxVoltage/2;
                            f[elements[l]->nodes[1]->num-1] -= longueur*region[elements[l]->region]->fluxVoltage/2;
                        }

                        if(region[elements[l]->region]->voltage != -1)
                        {
                            theta_k[elements[l]->nodes[0]->num-1] = region[elements[l]->region]->voltage;
                            theta_k[elements[l]->nodes[1]->num-1] = region[elements[l]->region]->voltage;
                        }
                    }
                }
            }
        }
    }

    if(method == DIRICHLETFLAG)//Including the Dirichlet condition on theta_k
    {
        for(unsigned int l = 0; l < elements.size(); l++)
        {
            if(elements[l]->type == 1)//If line
            {
                if(region.count(elements[l]->region) == 1)
                {
                    for(unsigned int j = 0; j < elements[l]->nodes.size(); j++)
                    {
                        if(thermalOrElectrical == THERMALFLAG)
                        {
                            if(region[elements[l]->region]->temperature != -1)
                            {
                                theta_k[elements[l]->nodes[j]->num-1] = region[elements[l]->region]->temperature;
                            }
                        }
                        else if(thermalOrElectrical == ELECTRICFLAG)
                        {
                            if(region[elements[l]->region]->voltage != -1)
                            {
                                theta_k[elements[l]->nodes[j]->num-1] = region[elements[l]->region]->voltage;
                            }
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
        NewtonRaphson(nodes, elements, physicals, theta_k, f, method, NodesCorresp, delta_theta_k, region, RHS, corner, border, conditions, thermalOrElectrical);

        //Initialize value for normRHS0
        if(iter==1)
            normRHS0 = gmm::vect_norm2(RHS);

        //Check the convergence criterion
        Criterion = End_Criterion(RHS, normRHS0, eps, type);
        //cout << "Iteration number " << iter << endl;
        iter++;
    }

    std::vector<double> qint(nodes.size());
    Internal_flux(theta_k, region, elements, qint, q_m_x, q_m_y, thermalOrElectrical);//Chargement de la solution pour la solution du flux

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
void f_function(std::vector<double> &f, std::vector<Node*> &nodes, std::vector<Element*> &elements, std::map<int,Parameter*> &region, FemFlag thermalOrElectrical, int constantProperty, int JouleEffect)
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
                if(thermalOrElectrical == THERMALFLAG)
                {
                    cons = region[elements[i]->region]->thermalGeneration;
					if(JouleEffect == 1) cons = cons + elements[i]->s_Me;
                }
                else if(thermalOrElectrical == ELECTRICFLAG)
                {
                    cons = region[elements[i]->region]->electricalGeneration;
                }
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
        else if(elements[i]->type == 3)//If quad
        {
            if(constantProperty != 0)
            {
                cons = constantProperty;
            }
            else
            {
                if(thermalOrElectrical == THERMALFLAG)
                {
                    cons = region[elements[i]->region]->thermalGeneration;
                }
                else if(thermalOrElectrical == ELECTRICFLAG)
                {
                    cons = region[elements[i]->region]->electricalGeneration;
                }
            }

            Node *n1 = elements[i]->nodes[0];
            Node *n2 = elements[i]->nodes[1];
            Node *n3 = elements[i]->nodes[2];
            Node *n4 = elements[i]->nodes[3];
            double x1 = n1->x;
            double y1 = n1->y;
            double x2 = n2->x;
            double y2 = n2->y;
            double x3 = n3->x;
            double y3 = n3->y;
            double x4 = n4->x;
            double y4 = n4->y;
            double intDetJ = (x1*y2)/2 - (x2*y1)/2 - (x1*y4)/2 + (x2*y3)/2 - (x3*y2)/2 + (x4*y1)/2 + (x3*y4)/2 - (x4*y3)/2;

            f[n1->num-1] += cons*intDetJ/4;
            f[n2->num-1] += cons*intDetJ/4;
            f[n3->num-1] += cons*intDetJ/4;
            f[n4->num-1] += cons*intDetJ/4;
        }
    }
}

void Average_flux(std::map<Node*, std::vector<double> > &solutionTemperature, std::map<int,Parameter*> &region, std::vector<Element*> &elements, std::vector<double> &q_Me, double vol, FemFlag thermalOrElectrical)
{
    gmm::dense_matrix<double> alpha(2,2); // matrice alpha
    gmm::dense_matrix<double> beta(2,2); // matrice beta pour la conductivité
    q_Me[0] = 0.0;
    q_Me[1] = 0.0;


    for(unsigned int i = 0; i < elements.size(); i++)//loop over the elements
    {
        if(elements[i]->type == 2)//If triangle
        {
            //cout << "In LOOP 1 !" << endl;
            Node *n1 = elements[i]->nodes[0];
            Node *n2 = elements[i]->nodes[1];
            Node *n3 = elements[i]->nodes[2];
            double x1 = n1->x;
            double y1 = n1->y;
            double x2 = n2->x;
            double y2 = n2->y;
            double x3 = n3->x;
            double y3 = n3->y;

			if(thermalOrElectrical == THERMALFLAG)
			{
		        if(region[elements[i]->region]->thermalConductivity[0]->name == "alpha")
		        {
		            alpha(0,0) = region[elements[i]->region]->thermalConductivity[0]->conductivity[0][0];
		            alpha(0,1) = region[elements[i]->region]->thermalConductivity[0]->conductivity[0][1];
		            alpha(1,0) = region[elements[i]->region]->thermalConductivity[0]->conductivity[1][0];
		            alpha(1,1) = region[elements[i]->region]->thermalConductivity[0]->conductivity[1][1];
		        }
		        else if(region[elements[i]->region]->thermalConductivity[1]->name == "alpha")
		        {
		            alpha(0,0) = region[elements[i]->region]->thermalConductivity[1]->conductivity[0][0];
		            alpha(0,1) = region[elements[i]->region]->thermalConductivity[1]->conductivity[0][1];
		            alpha(1,0) = region[elements[i]->region]->thermalConductivity[1]->conductivity[1][0];
		            alpha(1,1) = region[elements[i]->region]->thermalConductivity[1]->conductivity[1][1];
		        }

		        if(region[elements[i]->region]->thermalConductivity[0]->name == "beta")
		        {
		            beta(0,0) = region[elements[i]->region]->thermalConductivity[0]->conductivity[0][0];
		            beta(0,1) = region[elements[i]->region]->thermalConductivity[0]->conductivity[0][1];
		            beta(1,0) = region[elements[i]->region]->thermalConductivity[0]->conductivity[1][0];
		            beta(1,1) = region[elements[i]->region]->thermalConductivity[0]->conductivity[1][1];
		        }
		        else if(region[elements[i]->region]->thermalConductivity[1]->name == "beta")
		        {
		            beta(0,0) = region[elements[i]->region]->thermalConductivity[1]->conductivity[0][0];
		            beta(0,1) = region[elements[i]->region]->thermalConductivity[1]->conductivity[0][1];
		            beta(1,0) = region[elements[i]->region]->thermalConductivity[1]->conductivity[1][0];
		            beta(1,1) = region[elements[i]->region]->thermalConductivity[1]->conductivity[1][1];
		        }
			}//end if thermalflag

			else if(thermalOrElectrical == ELECTRICFLAG)
			{
		        if(region[elements[i]->region]->electricalConductivity[0]->name == "alpha")
		        {
		            alpha(0,0) = region[elements[i]->region]->electricalConductivity[0]->conductivity[0][0];
		            alpha(0,1) = region[elements[i]->region]->electricalConductivity[0]->conductivity[0][1];
		            alpha(1,0) = region[elements[i]->region]->electricalConductivity[0]->conductivity[1][0];
		            alpha(1,1) = region[elements[i]->region]->electricalConductivity[0]->conductivity[1][1];
		        }
		        else if(region[elements[i]->region]->electricalConductivity[1]->name == "alpha")
		        {
		            alpha(0,0) = region[elements[i]->region]->electricalConductivity[1]->conductivity[0][0];
		            alpha(0,1) = region[elements[i]->region]->electricalConductivity[1]->conductivity[0][1];
		            alpha(1,0) = region[elements[i]->region]->electricalConductivity[1]->conductivity[1][0];
		            alpha(1,1) = region[elements[i]->region]->electricalConductivity[1]->conductivity[1][1];
		        }

		        if(region[elements[i]->region]->electricalConductivity[0]->name == "beta")
		        {
		            beta(0,0) = region[elements[i]->region]->electricalConductivity[0]->conductivity[0][0];
		            beta(0,1) = region[elements[i]->region]->electricalConductivity[0]->conductivity[0][1];
		            beta(1,0) = region[elements[i]->region]->electricalConductivity[0]->conductivity[1][0];
		            beta(1,1) = region[elements[i]->region]->electricalConductivity[0]->conductivity[1][1];
		        }
		        else if(region[elements[i]->region]->electricalConductivity[1]->name == "beta")
		        {
		            beta(0,0) = region[elements[i]->region]->electricalConductivity[1]->conductivity[0][0];
		            beta(0,1) = region[elements[i]->region]->electricalConductivity[1]->conductivity[0][1];
		            beta(1,0) = region[elements[i]->region]->electricalConductivity[1]->conductivity[1][0];
		            beta(1,1) = region[elements[i]->region]->electricalConductivity[1]->conductivity[1][1];
		        }
			}//end if electricalflag

            std::vector<double> theta_nodes(3);//Temperature at nodes
            std::vector<double> sol;
            sol = solutionTemperature[n1];
            theta_nodes[0] = sol[0];
            sol = solutionTemperature[n2];
            theta_nodes[1] = sol[0];
            sol = solutionTemperature[n3];
            theta_nodes[2] = sol[0];
            double thetam = (theta_nodes[0]+theta_nodes[1]+theta_nodes[2])/3;//Average temperature on the element
            double detJ = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);//Determinant of the Jacobian matrix

            gmm::dense_matrix<double> kappa(2,2); //Conductivity tensor at given temperature

            //Assembly of the conductivity tensor kappa(k,j)
            for (unsigned int k = 0; k<2;k++)
            {
            	for(unsigned int j = 0 ; j<2; j++)
            	{
            		kappa(k,j) = alpha (k,j)*thetam+beta(k,j);
            	}
            }


            //Inverse of the Jacobian Matrix
            gmm::dense_matrix<double> inv_J(2,2);
            inv_J(0,0) = 1/detJ*(y3-y1);
            inv_J(0,1) = 1/detJ*(y1-y2);
            inv_J(1,0) = 1/detJ*(x1-x3);
            inv_J(1,1) = 1/detJ*(x2-x1);


            //Matrix containing the gradients of the shape functions in isoparametric coordinates
            gmm::dense_matrix<double> grad_phi(2,3);
            grad_phi(0,0) = -1;
            grad_phi(1,0) = -1;
            grad_phi(0,1) = 1;
            grad_phi(1,1) = 0;
            grad_phi(0,2) = 0;
            grad_phi(1,2) = 1;


            //Computation of q_Me
            gmm::dense_matrix<double> a1(2,2);
            gmm::mult(kappa, inv_J, a1);
            gmm::scale(a1, 0.5*detJ);
            std::vector<double> grad_phi_current(2);
            std::vector<double> a2(2);
            for(unsigned int k = 0; k<3; k++)
            {
                grad_phi_current[0] = grad_phi(0,k);
                grad_phi_current[1] = grad_phi(1,k);
                gmm::mult(a1, grad_phi_current, a2);
                gmm::scale(a2, theta_nodes[k]);
                gmm::add(a2,q_Me);
            }

        }//end if element = triangle


    }//end loop over the elements
    gmm::scale(q_Me, -1./vol);

}

void Average_Joule(std::map<Node*, std::vector<double> > &solutionTemperature, std::vector<Element*> &elements, std::map<int, Parameter*> &region, double &s_Me, double vol)
{
    gmm::dense_matrix<double> alpha(2,2); // matrice alpha
    gmm::dense_matrix<double> beta(2,2); // matrice beta pour la conductivité
	std::vector<double> q_Me(2);
	std::vector<double> E(2);

    for(unsigned int i = 0; i < elements.size(); i++)//loop over the elements
    {
        if(elements[i]->type == 2)//If triangle
        {
            //cout << "In LOOP 1 !" << endl;
            Node *n1 = elements[i]->nodes[0];
            Node *n2 = elements[i]->nodes[1];
            Node *n3 = elements[i]->nodes[2];
            double x1 = n1->x;
            double y1 = n1->y;
            double x2 = n2->x;
            double y2 = n2->y;
            double x3 = n3->x;
            double y3 = n3->y;

	        if(region[elements[i]->region]->electricalConductivity[0]->name == "alpha")
	        {
	            alpha(0,0) = region[elements[i]->region]->electricalConductivity[0]->conductivity[0][0];
	            alpha(0,1) = region[elements[i]->region]->electricalConductivity[0]->conductivity[0][1];
	            alpha(1,0) = region[elements[i]->region]->electricalConductivity[0]->conductivity[1][0];
	            alpha(1,1) = region[elements[i]->region]->electricalConductivity[0]->conductivity[1][1];
	        }
	        else if(region[elements[i]->region]->electricalConductivity[1]->name == "alpha")
	        {
	            alpha(0,0) = region[elements[i]->region]->electricalConductivity[1]->conductivity[0][0];
	            alpha(0,1) = region[elements[i]->region]->electricalConductivity[1]->conductivity[0][1];
	            alpha(1,0) = region[elements[i]->region]->electricalConductivity[1]->conductivity[1][0];
	            alpha(1,1) = region[elements[i]->region]->electricalConductivity[1]->conductivity[1][1];
	        }

	        if(region[elements[i]->region]->electricalConductivity[0]->name == "beta")
	        {
	            beta(0,0) = region[elements[i]->region]->electricalConductivity[0]->conductivity[0][0];
	            beta(0,1) = region[elements[i]->region]->electricalConductivity[0]->conductivity[0][1];
	            beta(1,0) = region[elements[i]->region]->electricalConductivity[0]->conductivity[1][0];
	            beta(1,1) = region[elements[i]->region]->electricalConductivity[0]->conductivity[1][1];
	        }
	        else if(region[elements[i]->region]->electricalConductivity[1]->name == "beta")
	        {
	            beta(0,0) = region[elements[i]->region]->electricalConductivity[1]->conductivity[0][0];
	            beta(0,1) = region[elements[i]->region]->electricalConductivity[1]->conductivity[0][1];
	            beta(1,0) = region[elements[i]->region]->electricalConductivity[1]->conductivity[1][0];
	            beta(1,1) = region[elements[i]->region]->electricalConductivity[1]->conductivity[1][1];
	        }

            std::vector<double> theta_nodes(3);//Temperature at nodes
            std::vector<double> sol;
            sol = solutionTemperature[n1];
            theta_nodes[0] = sol[0];
            sol = solutionTemperature[n2];
            theta_nodes[1] = sol[0];
            sol = solutionTemperature[n3];
            theta_nodes[2] = sol[0];
            double thetam = (theta_nodes[0]+theta_nodes[1]+theta_nodes[2])/3;//Average temperature on the element
            double detJ = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);//Determinant of the Jacobian matrix

            gmm::dense_matrix<double> kappa(2,2); //Conductivity tensor at given temperature

            //Assembly of the conductivity tensor kappa(k,j)
            for (unsigned int k = 0; k<2;k++)
            {
            	for(unsigned int j = 0 ; j<2; j++)
            	{
            		kappa(k,j) = alpha (k,j)*thetam+beta(k,j);
            	}
            }


            //Inverse of the Jacobian Matrix
            gmm::dense_matrix<double> inv_J(2,2);
            inv_J(0,0) = 1/detJ*(y3-y1);
            inv_J(0,1) = 1/detJ*(y1-y2);
            inv_J(1,0) = 1/detJ*(x1-x3);
            inv_J(1,1) = 1/detJ*(x2-x1);


            //Matrix containing the gradients of the shape functions in isoparametric coordinates
            gmm::dense_matrix<double> grad_phi(2,3);
            grad_phi(0,0) = -1;
            grad_phi(1,0) = -1;
            grad_phi(0,1) = 1;
            grad_phi(1,1) = 0;
            grad_phi(0,2) = 0;
            grad_phi(1,2) = 1;


            //Computation of q_Me
            gmm::dense_matrix<double> a1(2,2);
            gmm::dense_matrix<double> a1_E(2,2);
            gmm::mult(kappa, inv_J, a1);
            gmm::scale(a1, 0.5*detJ);
            gmm::copy(inv_J, a1_E);
            gmm::scale(a1_E, 0.5*detJ);
            std::vector<double> grad_phi_current(2);
            std::vector<double> a2(2);
            std::vector<double> a2_E(2);
            for(unsigned int k = 0; k<3; k++)
            {
                grad_phi_current[0] = grad_phi(0,k);
                grad_phi_current[1] = grad_phi(1,k);
                gmm::mult(a1, grad_phi_current, a2);
                gmm::scale(a2, theta_nodes[k]);
                gmm::add(a2,q_Me);
				gmm::mult(a1_E, grad_phi_current, a2_E);
                gmm::scale(a2_E, theta_nodes[k]);
                gmm::add(a2_E,E);

            }
        }//end if element = triangle


    }//end loop over the elements
    gmm::scale(q_Me, -1./vol);
	//cout << "q_Me " << q_Me << endl;
    gmm::scale(E, -1./vol);
	//cout << "E " << E << endl;
	s_Me = gmm::vect_sp(E, q_Me);

}

void Average_Joule_FE1(std::map<Node*, std::vector<double> > &solutionTemperature, std::vector<Element*> &elements, std::map<int, Parameter*> &region)
{
    gmm::dense_matrix<double> alpha(2,2); // matrice alpha
    gmm::dense_matrix<double> beta(2,2); // matrice beta pour la conductivité
	std::vector<double> J(2);
	std::vector<double> E(2);
	std::vector<double> gradT(2);

    for(unsigned int i = 0; i < elements.size(); i++)//loop over the elements
    {
        if(elements[i]->type == 2)//If triangle
        {
            //cout << "In LOOP 1 !" << endl;
            Node *n1 = elements[i]->nodes[0];
            Node *n2 = elements[i]->nodes[1];
            Node *n3 = elements[i]->nodes[2];
            double x1 = n1->x;
            double y1 = n1->y;
            double x2 = n2->x;
            double y2 = n2->y;
            double x3 = n3->x;
            double y3 = n3->y;

	        if(region[elements[i]->region]->electricalConductivity[0]->name == "alpha")
	        {
	            alpha(0,0) = region[elements[i]->region]->electricalConductivity[0]->conductivity[0][0];
	            alpha(0,1) = region[elements[i]->region]->electricalConductivity[0]->conductivity[0][1];
	            alpha(1,0) = region[elements[i]->region]->electricalConductivity[0]->conductivity[1][0];
	            alpha(1,1) = region[elements[i]->region]->electricalConductivity[0]->conductivity[1][1];
	        }
	        else if(region[elements[i]->region]->electricalConductivity[1]->name == "alpha")
	        {
	            alpha(0,0) = region[elements[i]->region]->electricalConductivity[1]->conductivity[0][0];
	            alpha(0,1) = region[elements[i]->region]->electricalConductivity[1]->conductivity[0][1];
	            alpha(1,0) = region[elements[i]->region]->electricalConductivity[1]->conductivity[1][0];
	            alpha(1,1) = region[elements[i]->region]->electricalConductivity[1]->conductivity[1][1];
	        }

	        if(region[elements[i]->region]->electricalConductivity[0]->name == "beta")
	        {
	            beta(0,0) = region[elements[i]->region]->electricalConductivity[0]->conductivity[0][0];
	            beta(0,1) = region[elements[i]->region]->electricalConductivity[0]->conductivity[0][1];
	            beta(1,0) = region[elements[i]->region]->electricalConductivity[0]->conductivity[1][0];
	            beta(1,1) = region[elements[i]->region]->electricalConductivity[0]->conductivity[1][1];
	        }
	        else if(region[elements[i]->region]->electricalConductivity[1]->name == "beta")
	        {
	            beta(0,0) = region[elements[i]->region]->electricalConductivity[1]->conductivity[0][0];
	            beta(0,1) = region[elements[i]->region]->electricalConductivity[1]->conductivity[0][1];
	            beta(1,0) = region[elements[i]->region]->electricalConductivity[1]->conductivity[1][0];
	            beta(1,1) = region[elements[i]->region]->electricalConductivity[1]->conductivity[1][1];
	        }

            std::vector<double> theta_nodes(3);//Temperature at nodes
            std::vector<double> sol;
            sol = solutionTemperature[n1];
            theta_nodes[0] = sol[0];
            sol = solutionTemperature[n2];
            theta_nodes[1] = sol[0];
            sol = solutionTemperature[n3];
            theta_nodes[2] = sol[0];
			double u1 = theta_nodes[0];
			double u2 = theta_nodes[1];
			double u3 = theta_nodes[2];

            double thetam = (theta_nodes[0]+theta_nodes[1]+theta_nodes[2])/3;//Average temperature on the element
            double detJ = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);//Determinant of the Jacobian matrix

            gmm::dense_matrix<double> kappa(2,2); //Conductivity tensor at given temperature

            //Assembly of the conductivity tensor kappa(k,j)
            for (unsigned int k = 0; k<2;k++)
            {
            	for(unsigned int j = 0 ; j<2; j++)
            	{
            		kappa(k,j) = alpha (k,j)*thetam+beta(k,j);
            	}
            }


            //Inverse of the Jacobian Matrix
            gmm::dense_matrix<double> inv_J(2,2);
            inv_J(0,0) = 1/detJ*(y3-y1);
            inv_J(0,1) = 1/detJ*(y1-y2);
            inv_J(1,0) = 1/detJ*(x1-x3);
            inv_J(1,1) = 1/detJ*(x2-x1);


            //Matrix containing the gradients of the shape functions in isoparametric coordinates
			std::vector<double> gradPhi1_red(2);
			std::vector<double> gradPhi2_red(2);
			std::vector<double> gradPhi3_red(2);
			std::vector<double> gradPhi1(2);
			std::vector<double> gradPhi2(2);
			std::vector<double> gradPhi3(2);

			gradPhi1_red[0] = -1.0;
			gradPhi1_red[1] = -1.0;
			gradPhi2_red[0] = 1.0;
			gradPhi2_red[1] = 0.0;
			gradPhi3_red[0] = 0.0;
			gradPhi3_red[1] = 1.0;

		    gmm::mult(inv_J, gradPhi1_red, gradPhi1);
		    gmm::mult(inv_J, gradPhi2_red, gradPhi2);
		    gmm::mult(inv_J, gradPhi3_red, gradPhi3);

			gradT[0] = u1*gradPhi1[0] +  u2*gradPhi2[0] + u3*gradPhi3[0];
	    	gradT[1] = u1*gradPhi1[1] +  u2*gradPhi2[1] + u3*gradPhi3[1];

			gmm::copy(gradT,E);
			gmm::scale(E,-1);
			gmm::mult(kappa,E,J);
			elements[i]->s_Me = gmm::vect_sp(J,E);

        }//end if element = triangle


    }//end loop over the elements

}
