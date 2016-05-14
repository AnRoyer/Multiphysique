#include <iostream>
#include <cstdio>
#include <vector>
#include <map>
#include <mpi.h>
#include <string>
#include <cmath>
#include "gmm/gmm.h"
#include "gmshio.h"
#include "physicalio.h"
#include "NewtonRaphson.h"
#include "fem.h"
#ifdef GMM_USES_MUMPS
#include "gmm/gmm_MUMPS_interface.h"
#endif
using namespace std;

/*---------------------------------------------NEWTON RAPHSON ROUTINE--------------------------------------------------------------*/

void NewtonRaphson(std::vector<Node*> &nodes, std::vector<Element*> &elements, std::vector<Physical*> &physicals, std::vector<double> &theta_k,
					std::vector<double> &qext, FemFlag method, std::map<Node*, Node*> &NodesCorresp,
					std::vector<double> &delta_theta_k, std::map<int,Parameter*> &region, std::vector<double> &RHS,
					NodeCorner &corner, NodeBorder &border, Periodique &conditions, FemFlag thermalOrElectrical)
{
    std::vector<double> qint(nodes.size()); //Internal flux vector
	gmm::row_matrix< gmm::wsvector<double> > KT_tmp(nodes.size(), nodes.size());//Tangent matrix used for build-in
    gmm::csr_matrix<double> KT;//Tangent matrix used in system solving

    //Internal Flux : qint
    std::vector<double> q_m_x(0); //Internal flux vector
    std::vector<double> q_m_y(0); //Internal flux vector

	Internal_flux(theta_k, region, elements, qint, q_m_x, q_m_y, thermalOrElectrical);

    //TANGENT STIFFNESS MATRIX :  KT
    Tangent_Stiffness_Matrix(theta_k, region, elements, KT_tmp, NodesCorresp, nodes, thermalOrElectrical);

    //Including the Dirichlet boundary conditions
    if(method == DIRICHLETFLAG)
    {
        //Dirichlet conditions in KT_tmp
        for(unsigned int i = 0; i < elements.size(); i++)
        {
            if(elements[i]->type == 1)//If line
            {
                if(thermalOrElectrical == THERMALFLAG)
                {
                    if(region.count(elements[i]->region) == 1 && region[elements[i]->region]->temperature != -1)//If linesRegion contains elements[i]->region
                    {
                        for(unsigned int j = 0; j < elements[i]->nodes.size(); j++)
                        {
                            for(unsigned int k = 0; k < nodes.size(); k++)
                            {
                                if(k == elements[i]->nodes[j]->num-1)
                                {
                                    KT_tmp(elements[i]->nodes[j]->num-1, k) = 1;
                                }
                                else
                                {
                                    KT_tmp(elements[i]->nodes[j]->num-1, k) = 0;
                                }
                            }
                        }
                    }
                }
                else if(thermalOrElectrical == ELECTRICFLAG)//If line
                {
                    if(region.count(elements[i]->region) == 1 && region[elements[i]->region]->voltage != -1)//If linesRegion contains elements[i]->region
                    {
                        for(unsigned int j = 0; j < elements[i]->nodes.size(); j++)
                        {
                            for(unsigned int k = 0; k < nodes.size(); k++)
                            {
                                if(k == elements[i]->nodes[j]->num-1)
                                {
                                    KT_tmp(elements[i]->nodes[j]->num-1, k) = 1;
                                }
                                else
                                {
                                    KT_tmp(elements[i]->nodes[j]->num-1, k) = 0;
                                }
                            }
                        }
                    }
                }
            }
        }

        //Building of the Right Hand Side
        for (unsigned int i=0 ; i<nodes.size(); i++)
        {
            RHS[i] = qint[i]+qext[i];
        }

        //Including the Dirichlet boundary conditions in RHS
        for(unsigned int i = 0; i < elements.size(); i++)
        {
            if(elements[i]->type == 1)//If line
            {
                if(thermalOrElectrical == THERMALFLAG)
                {
                    if(region.count(elements[i]->region) == 1 && region[elements[i]->region]->temperature != -1)//If linesRegion contains elements[i]->region
                    {
                        for(unsigned int j = 0; j < elements[i]->nodes.size(); j++)
                        {
                            RHS[elements[i]->nodes[j]->num-1] = 0;
                        }
                    }
                }
                else if(thermalOrElectrical == ELECTRICFLAG)
                {
                    if(region.count(elements[i]->region) == 1 && region[elements[i]->region]->voltage != -1)//If linesRegion contains elements[i]->region
                    {
                        for(unsigned int j = 0; j < elements[i]->nodes.size(); j++)
                        {
                            RHS[elements[i]->nodes[j]->num-1] = 0;
                        }
                    }
                }
            }
        }
    }
    else if(method == VONNEUMANNFLAG)
    {
        //Dirichlet conditions in KT_tmp
        for(unsigned int i = 0; i < elements.size(); i++)
        {
            if(elements[i]->type == 1)//If line
            {
                if(thermalOrElectrical == THERMALFLAG)
                {
                    if(region.count(elements[i]->region) == 1 && region[elements[i]->region]->temperature != -1)//If linesRegion contains elements[i]->region
                    {
                        for(unsigned int j = 0; j < elements[i]->nodes.size(); j++)
                        {
                            for(unsigned int k = 0; k < nodes.size(); k++)
                            {
                                if(k == elements[i]->nodes[j]->num-1)
                                {
                                    KT_tmp(elements[i]->nodes[j]->num-1, k) = 1;
                                }
                                else
                                {
                                    KT_tmp(elements[i]->nodes[j]->num-1, k) = 0;
                                }
                            }
                        }
                    }
                }
                else if(thermalOrElectrical == ELECTRICFLAG)
                {
                    if(region.count(elements[i]->region) == 1 && region[elements[i]->region]->voltage != -1)//If linesRegion contains elements[i]->region
                    {
                        for(unsigned int j = 0; j < elements[i]->nodes.size(); j++)
                        {
                            for(unsigned int k = 0; k < nodes.size(); k++)
                            {
                                if(k == elements[i]->nodes[j]->num-1)
                                {
                                    KT_tmp(elements[i]->nodes[j]->num-1, k) = 1;
                                }
                                else
                                {
                                    KT_tmp(elements[i]->nodes[j]->num-1, k) = 0;
                                }
                            }
                        }
                    }
                }
            }
        }

        //Building of the Right Hand Side
        for (unsigned int i=0 ; i<nodes.size(); i++)
        {
            RHS[i] = qint[i]+qext[i];
        }

        //Including the Dirichlet boundary conditions in RHS
        for(unsigned int i = 0; i < elements.size(); i++)
        {
            if(elements[i]->type == 1)//If line
            {
                if(thermalOrElectrical == THERMALFLAG)
                {
                    if(region.count(elements[i]->region) == 1 && region[elements[i]->region]->temperature != -1)//If linesRegion contains elements[i]->region
                    {
                        for(unsigned int j = 0; j < elements[i]->nodes.size(); j++)
                        {
                            RHS[elements[i]->nodes[j]->num-1] = 0;
                        }
                    }
                }
                else if(thermalOrElectrical == ELECTRICFLAG)
                {
                    if(region.count(elements[i]->region) == 1 && region[elements[i]->region]->voltage != -1)//If linesRegion contains elements[i]->region
                    {
                        for(unsigned int j = 0; j < elements[i]->nodes.size(); j++)
                        {
                            RHS[elements[i]->nodes[j]->num-1] = 0;
                        }
                    }
                }
            }
        }
    }
    else if(method == PERIODICFLAG)
    {
        double vol=0;
        double lx = abs(corner.C2->x - corner.C1->x);
        double ly = abs(corner.C4->y - corner.C1->y);


        double Tavg = 0;

        if(thermalOrElectrical == THERMALFLAG)
        {
            Tavg = conditions.meanTemperature;
        }
        else if(thermalOrElectrical == ELECTRICFLAG)
        {
            Tavg = conditions.meanVoltage;
        }

        double gradAvg_x = conditions.xGradient;
        double gradAvg_y = conditions.yGradient;

        //Condition correspondant au noeud 1: temperature moyenne.
        unsigned int numC1 = corner.C1->num-1;
        vector<double> c(nodes.size());//coefficients

        if(thermalOrElectrical == THERMALFLAG)
        {
            f_function(c, nodes, elements, region, THERMALFLAG, 1, 0);
        }
        else if(thermalOrElectrical == ELECTRICFLAG)
        {
            f_function(c, nodes, elements, region, ELECTRICFLAG, 1, 0);
        }

        for(unsigned int j=0; j<nodes.size(); j++)
        {
            KT_tmp(numC1,j) = c[j];
            vol += c[j]; //faux
        }

        double somTheta = 0;
        for(unsigned int j=0; j<nodes.size(); j++)
        {
            somTheta += c[j]*theta_k[j];
        }

        qext[numC1] = vol*Tavg - somTheta;//condition
        qint[numC1] = 0;

        //Conditions sur les noeuds 2 3 et 4.
        unsigned int numC2 = corner.C2->num-1;
        unsigned int numC3 = corner.C3->num-1;
        unsigned int numC4 = corner.C4->num-1;

        for(unsigned int j=0; j<nodes.size(); j++)
        {
            KT_tmp(numC2,j) = 0;
            KT_tmp(numC3,j) = 0;
            KT_tmp(numC4,j) = 0;

            if(j == numC1)
            {
                KT_tmp(numC2,j) = -1;
                KT_tmp(numC3,j) = -1;
                KT_tmp(numC4,j) = -1;
            }

            if(j == numC2)
            {
                KT_tmp(numC2,j) = 1;
            }
            else if(j == numC3)
            {
                KT_tmp(numC3,j) = 1;
            }
            else if(j == numC4)
            {
                KT_tmp(numC4,j) = 1;
            }

        }

        //conditions sur qext correspondants aux noeuds des coins.
        qext[numC2] = gradAvg_x*lx + theta_k[numC1] - theta_k[numC2];
        qext[numC3] = (gradAvg_x*lx + gradAvg_y*ly) + theta_k[numC1] - theta_k[numC3];
        qext[numC4] = gradAvg_y*ly + theta_k[numC1] - theta_k[numC4];

        qint[numC2] = 0;
        qint[numC3] = 0;
        qint[numC4] = 0;

        //conditions sur qext correspondant aux noeuds à droite et en haut
        for(unsigned int i=0;i<border.RightNodes.size();i++)
        {
            unsigned int numNode = border.RightNodes[i]->num-1;
            for(unsigned int j=0; j<nodes.size(); j++)
            {
                KT_tmp(numNode,j) = 0;

                if(j == NodesCorresp[border.RightNodes[i]]->num-1)
                {
                    KT_tmp(numNode,j) = -1;
                }

                if(j == numNode)
                {
                    KT_tmp(numNode,j) = 1;
                }

            }

            qint[NodesCorresp[border.RightNodes[i]]->num-1] += qint[numNode];
            qext[NodesCorresp[border.RightNodes[i]]->num-1] = 0;
            qext[numNode] = gradAvg_x*lx + theta_k[NodesCorresp[border.RightNodes[i]]->num-1] - theta_k[numNode];
            qint[numNode] = 0;
        }

        for(unsigned int i=0;i<border.TopNodes.size();i++)
        {
            unsigned int numNode = border.TopNodes[i]->num-1;

            for(unsigned int j=0; j<nodes.size(); j++)
            {
                KT_tmp(numNode,j) = 0;

                if(j == NodesCorresp[border.TopNodes[i]]->num-1)
                {
                    KT_tmp(numNode,j) = -1;
                }

                if(j == numNode)
                {
                    KT_tmp(numNode,j) = 1;
                }

            }

            qint[NodesCorresp[border.TopNodes[i]]->num-1] += qint[numNode];
            qext[NodesCorresp[border.TopNodes[i]]->num-1] = 0;
            qext[numNode] = gradAvg_y*ly + theta_k[NodesCorresp[border.TopNodes[i]]->num-1] - theta_k[numNode];
            qint[numNode] = 0;
        }

        //Building of the Right Hand Side
        for (unsigned int i=0 ; i<nodes.size(); i++)
        {
            RHS[i] = qint[i]+qext[i];
        }
    }

    //Solving the system
    gmm::copy(KT_tmp, KT);
#ifdef GMM_USES_MUMPS
    //std::cout << "solving linear system with MUMPS\n";
    gmm::MUMPS_solve(KT, delta_theta_k, RHS);
#else
    //std::cout << "solving linear system with gmm::lu_solve\n";
    gmm::lu_solve(KT, delta_theta_k, RHS);
#endif

    if(method == DIRICHLETFLAG)
    {
        //Reincluding the Dirichlet boundary conditions in delta_theta_k to avoid little error in system resolution
        for(unsigned int i = 0; i < elements.size(); i++)
        {
            if(elements[i]->type == 1)//If line
            {
                if(thermalOrElectrical == THERMALFLAG)
                {
                    if(region.count(elements[i]->region) == 1 && region[elements[i]->region]->temperature != -1)//If linesRegion contains elements[i]->region
                    {
                        for(unsigned int j = 0; j < elements[i]->nodes.size(); j++)
                        {
                            RHS[elements[i]->nodes[j]->num-1] = 0;
                        }
                    }
                }
                else if(thermalOrElectrical == ELECTRICFLAG)
                {
                    if(region.count(elements[i]->region) == 1 && region[elements[i]->region]->voltage != -1)//If linesRegion contains elements[i]->region
                    {
                        for(unsigned int j = 0; j < elements[i]->nodes.size(); j++)
                        {
                            RHS[elements[i]->nodes[j]->num-1] = 0;
                        }
                    }
                }
            }
        }
    }

    for (unsigned int i=0; i<nodes.size(); i++)
	{
		theta_k[i] += delta_theta_k[i];
    }

}//end of Newton Raphson


/*----------------------------------------------------------INTERNAL FLUX ROUTINE----------------------------------------------------------------------*/
void Internal_flux(std::vector<double> &theta_k, std::map<int,Parameter*> &region, std::vector<Element*> &elements,
                   std::vector<double> &qint, std::vector<double> &q_m_x, std::vector<double> &q_m_y, FemFlag thermalOrElectrical)
{
    gmm::dense_matrix<double> alpha(2,2); // matrice alpha
    gmm::dense_matrix<double> beta(2,2); // matrice beta pour la conductivité

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
            double qint_I; //Contribution of a given node to qint

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
            }
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
            }

            std::vector<double> theta_nodes(3);//Temperature at nodes
            theta_nodes[0] = theta_k[n1->num-1];
            theta_nodes[1] = theta_k[n2->num-1];
            theta_nodes[2] = theta_k[n3->num-1];
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

            //Computation of qint

            std::vector<double> tmp_vec(2,0);//Temporary vector in the computation of qint
            std::vector<double> tmp_vec_2(2,0);//Temporary vector in the computation of qint


            for(unsigned int j = 0; j<3 ; j++)
            {
                tmp_vec[0] += theta_nodes[j]*grad_phi(0,j);
                tmp_vec[1] += theta_nodes[j]*grad_phi(1,j);
            }

            for(unsigned int j = 0; j<2; j++)
            {
                for(unsigned int k = 0; k<2 ; k++)
                {
                    tmp_vec_2[j] += inv_J(j,k)*tmp_vec[k];
                }
            }

            tmp_vec[0]=0;
            tmp_vec[1]=0;
            for(unsigned int j = 0; j<2; j++)
            {
                for(unsigned int k = 0; k<2 ; k++)
                {
                    tmp_vec[j]+= kappa(j,k)*tmp_vec_2[k];
                }
            }

            if(q_m_x.size() != 0 && q_m_y.size() != 0)
            {
                q_m_x[n1->num-1] = -tmp_vec[0];
                q_m_x[n2->num-1] = -tmp_vec[0];
                q_m_x[n3->num-1] = -tmp_vec[0];

                q_m_y[n1->num-1] = -tmp_vec[1];
                q_m_y[n2->num-1] = -tmp_vec[1];
                q_m_y[n3->num-1] = -tmp_vec[1];
            }


            for(unsigned int l = 0; l<3 ; l++)//For each element, there are three contributions to qint at distinct nodes
            {
                tmp_vec_2[0]=0;
                tmp_vec_2[1]=0;
                //cout << "LOOP 3 : l = "<< l << endl;
            	for(unsigned int j = 0; j<2; j++)
                {
                     //cout << "LOOP 4 : j = "<< j << endl;
            		for(unsigned int k = 0; k<2 ; k++)
            		{
            		    //cout << "LOOP 5 : k = "<< k << endl;
            			tmp_vec_2[j] += inv_J(j,k)*grad_phi(k,l);
            		}
            	}



            	qint_I = -detJ/2*(tmp_vec[0]*tmp_vec_2[0]+tmp_vec[1]*tmp_vec_2[1]);

                //qint_i assigned to a different node in the qint vector
            	if(l+1 == 1)
                    qint[n1->num-1] += qint_I;
                if(l+1 == 2)
                    qint[n2->num-1] += qint_I;
                if(l+1 == 3)
                    qint[n3->num-1] += qint_I;
            }//end for
        }//end if element = triangle
        else if(elements[i]->type == 3)//If quad
        {
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
            }
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
            }

            std::vector<double> theta_nodes(4);//Temperature at nodes
            theta_nodes[0] = theta_k[n1->num-1];
            theta_nodes[1] = theta_k[n2->num-1];
            theta_nodes[2] = theta_k[n3->num-1];
            theta_nodes[3] = theta_k[n4->num-1];

            //Average temperature on the element
            double thetam = (theta_nodes[1]*((x1*y2)/4 - (x2*y1)/4 - (x1*y3)/4 + (x3*y1)/4 + (x2*y3)/4 - (x3*y2)/4)
                             + theta_nodes[0]*((x1*y2)/4 - (x2*y1)/4 - (x1*y4)/4 + (x4*y1)/4 + (x2*y4)/4 - (x4*y2)/4)
                             + theta_nodes[2]*((x1*y3)/4 - (x3*y1)/4 - (x1*y4)/4 + (x4*y1)/4 + (x3*y4)/4 - (x4*y3)/4)
                             + theta_nodes[3]*((x2*y3)/4 - (x3*y2)/4 - (x2*y4)/4 + (x4*y2)/4 + (x3*y4)/4 - (x4*y3)/4))
                             /((x1*y2)/2 - (x2*y1)/2 - (x1*y4)/2 + (x2*y3)/2 - (x3*y2)/2 + (x4*y1)/2 + (x3*y4)/2 - (x4*y3)/2);


            //Assembly of the conductivity tensor on the element
            gmm::dense_matrix<double> kappa(2,2);

            for (unsigned int k = 0; k<2;k++)
            {
                for(unsigned int j = 0 ; j<2; j++)
                {
                    kappa(k,j) = alpha(k,j)*thetam+beta(k,j);
                }
            }

            double k11 = kappa(0,0);
            double k12 = kappa(0,1);
            double k21 = kappa(1,0);
            double k22 = kappa(1,1);

            double t1 = theta_nodes[0];
            double t2 = theta_nodes[1];
            double t3 = theta_nodes[2];
            double t4 = theta_nodes[3];

            if(q_m_x.size() != 0 && q_m_y.size() != 0)
            {
                q_m_x[n1->num-1] = (k12*t1*x2 - k12*t2*x1 - k12*t1*x4 + k12*t4*x1 + k12*t2*x4 - k12*t4*x2 - k11*t1*y2 + k11*t2*y1 + k11*t1*y4 - k11*t4*y1 - k11*t2*y4 + k11*t4*y2)/(x1*y2 - x2*y1 - x1*y4 + x4*y1 + x2*y4 - x4*y2);
                q_m_x[n2->num-1] = (k12*t1*x2 - k12*t2*x1 - k12*t1*x3 + k12*t3*x1 + k12*t2*x3 - k12*t3*x2 - k11*t1*y2 + k11*t2*y1 + k11*t1*y3 - k11*t3*y1 - k11*t2*y3 + k11*t3*y2)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2);
                q_m_x[n3->num-1] = (k12*t2*x3 - k12*t3*x2 - k12*t2*x4 + k12*t4*x2 + k12*t3*x4 - k12*t4*x3 - k11*t2*y3 + k11*t3*y2 + k11*t2*y4 - k11*t4*y2 - k11*t3*y4 + k11*t4*y3)/(x2*y3 - x3*y2 - x2*y4 + x4*y2 + x3*y4 - x4*y3);
                q_m_x[n4->num-1] = (k12*t1*x3 - k12*t3*x1 - k12*t1*x4 + k12*t4*x1 + k12*t3*x4 - k12*t4*x3 - k11*t1*y3 + k11*t3*y1 + k11*t1*y4 - k11*t4*y1 - k11*t3*y4 + k11*t4*y3)/(x1*y3 - x3*y1 - x1*y4 + x4*y1 + x3*y4 - x4*y3);

                q_m_y[n1->num-1] = (k22*t1*x2 - k22*t2*x1 - k22*t1*x4 + k22*t4*x1 + k22*t2*x4 - k22*t4*x2 - k21*t1*y2 + k21*t2*y1 + k21*t1*y4 - k21*t4*y1 - k21*t2*y4 + k21*t4*y2)/(x1*y2 - x2*y1 - x1*y4 + x4*y1 + x2*y4 - x4*y2);
                q_m_y[n2->num-1] = (k22*t1*x2 - k22*t2*x1 - k22*t1*x3 + k22*t3*x1 + k22*t2*x3 - k22*t3*x2 - k21*t1*y2 + k21*t2*y1 + k21*t1*y3 - k21*t3*y1 - k21*t2*y3 + k21*t3*y2)/(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2);
                q_m_y[n3->num-1] = (k22*t2*x3 - k22*t3*x2 - k22*t2*x4 + k22*t4*x2 + k22*t3*x4 - k22*t4*x3 - k21*t2*y3 + k21*t3*y2 + k21*t2*y4 - k21*t4*y2 - k21*t3*y4 + k21*t4*y3)/(x2*y3 - x3*y2 - x2*y4 + x4*y2 + x3*y4 - x4*y3);
                q_m_y[n4->num-1] = (k22*t1*x3 - k22*t3*x1 - k22*t1*x4 + k22*t4*x1 + k22*t3*x4 - k22*t4*x3 - k21*t1*y3 + k21*t3*y1 + k21*t1*y4 - k21*t4*y1 - k21*t3*y4 + k21*t4*y3)/(x1*y3 - x3*y1 - x1*y4 + x4*y1 + x3*y4 - x4*y3);
            }

            qint[n1->num-1] += -(k22*t1*x2*x2 + k22*t1*x4*x4 - k22*t3*x2*x2 - k22*t3*x4*x4 + k11*t1*y2*y2 + k11*t1*y4*y4 - k11*t3*y2*y2 - k11*t3*y4*y4 - k22*t2*x1*x2 - 2*k22*t1*x2*x4 + k22*t2*x1*x4 + k22*t2*x2*x3 + k22*t4*x1*x2 - k22*t2*x3*x4 + 2*k22*t3*x2*x4 - k22*t4*x1*x4 - k22*t4*x2*x3 + k22*t4*x3*x4 - k12*t1*x2*y2 + k12*t2*x1*y2 + k12*t1*x2*y4 + k12*t1*x4*y2 - k12*t2*x1*y4 - k12*t2*x3*y2 + k12*t3*x2*y2 - k12*t4*x1*y2 - k12*t1*x4*y4 + k12*t2*x3*y4 - k12*t3*x2*y4 - k12*t3*x4*y2 + k12*t4*x1*y4 + k12*t4*x3*y2 + k12*t3*x4*y4 - k12*t4*x3*y4 - k21*t1*x2*y2 + k21*t2*x2*y1 + k21*t1*x2*y4 + k21*t1*x4*y2 - k21*t2*x2*y3 - k21*t2*x4*y1 + k21*t3*x2*y2 - k21*t4*x2*y1 - k21*t1*x4*y4 + k21*t2*x4*y3 - k21*t3*x2*y4 - k21*t3*x4*y2 + k21*t4*x2*y3 + k21*t4*x4*y1 + k21*t3*x4*y4 - k21*t4*x4*y3 - k11*t2*y1*y2 - 2*k11*t1*y2*y4 + k11*t2*y1*y4 + k11*t2*y2*y3 + k11*t4*y1*y2 - k11*t2*y3*y4 + 2*k11*t3*y2*y4 - k11*t4*y1*y4 - k11*t4*y2*y3 + k11*t4*y3*y4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));
            qint[n2->num-1] += -(k22*t2*x1*x1 + k22*t2*x3*x3 - k22*t4*x1*x1 - k22*t4*x3*x3 + k11*t2*y1*y1 + k11*t2*y3*y3 - k11*t4*y1*y1 - k11*t4*y3*y3 - k22*t1*x1*x2 + k22*t1*x1*x4 + k22*t1*x2*x3 - 2*k22*t2*x1*x3 + k22*t3*x1*x2 - k22*t1*x3*x4 - k22*t3*x1*x4 - k22*t3*x2*x3 + 2*k22*t4*x1*x3 + k22*t3*x3*x4 + k12*t1*x2*y1 - k12*t2*x1*y1 - k12*t1*x2*y3 - k12*t1*x4*y1 + k12*t2*x1*y3 + k12*t2*x3*y1 - k12*t3*x2*y1 + k12*t4*x1*y1 + k12*t1*x4*y3 - k12*t2*x3*y3 + k12*t3*x2*y3 + k12*t3*x4*y1 - k12*t4*x1*y3 - k12*t4*x3*y1 - k12*t3*x4*y3 + k12*t4*x3*y3 + k21*t1*x1*y2 - k21*t2*x1*y1 - k21*t1*x1*y4 - k21*t1*x3*y2 + k21*t2*x1*y3 + k21*t2*x3*y1 - k21*t3*x1*y2 + k21*t4*x1*y1 + k21*t1*x3*y4 - k21*t2*x3*y3 + k21*t3*x1*y4 + k21*t3*x3*y2 - k21*t4*x1*y3 - k21*t4*x3*y1 - k21*t3*x3*y4 + k21*t4*x3*y3 - k11*t1*y1*y2 + k11*t1*y1*y4 + k11*t1*y2*y3 - 2*k11*t2*y1*y3 + k11*t3*y1*y2 - k11*t1*y3*y4 - k11*t3*y1*y4 - k11*t3*y2*y3 + 2*k11*t4*y1*y3 + k11*t3*y3*y4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));
            qint[n3->num-1] += (k22*t1*x2*x2 + k22*t1*x4*x4 - k22*t3*x2*x2 - k22*t3*x4*x4 + k11*t1*y2*y2 + k11*t1*y4*y4 - k11*t3*y2*y2 - k11*t3*y4*y4 - k22*t2*x1*x2 - 2*k22*t1*x2*x4 + k22*t2*x1*x4 + k22*t2*x2*x3 + k22*t4*x1*x2 - k22*t2*x3*x4 + 2*k22*t3*x2*x4 - k22*t4*x1*x4 - k22*t4*x2*x3 + k22*t4*x3*x4 - k12*t1*x2*y2 + k12*t2*x1*y2 + k12*t1*x2*y4 + k12*t1*x4*y2 - k12*t2*x1*y4 - k12*t2*x3*y2 + k12*t3*x2*y2 - k12*t4*x1*y2 - k12*t1*x4*y4 + k12*t2*x3*y4 - k12*t3*x2*y4 - k12*t3*x4*y2 + k12*t4*x1*y4 + k12*t4*x3*y2 + k12*t3*x4*y4 - k12*t4*x3*y4 - k21*t1*x2*y2 + k21*t2*x2*y1 + k21*t1*x2*y4 + k21*t1*x4*y2 - k21*t2*x2*y3 - k21*t2*x4*y1 + k21*t3*x2*y2 - k21*t4*x2*y1 - k21*t1*x4*y4 + k21*t2*x4*y3 - k21*t3*x2*y4 - k21*t3*x4*y2 + k21*t4*x2*y3 + k21*t4*x4*y1 + k21*t3*x4*y4 - k21*t4*x4*y3 - k11*t2*y1*y2 - 2*k11*t1*y2*y4 + k11*t2*y1*y4 + k11*t2*y2*y3 + k11*t4*y1*y2 - k11*t2*y3*y4 + 2*k11*t3*y2*y4 - k11*t4*y1*y4 - k11*t4*y2*y3 + k11*t4*y3*y4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));
            qint[n4->num-1] += (k22*t2*x1*x1 + k22*t2*x3*x3 - k22*t4*x1*x1 - k22*t4*x3*x3 + k11*t2*y1*y1 + k11*t2*y3*y3 - k11*t4*y1*y1 - k11*t4*y3*y3 - k22*t1*x1*x2 + k22*t1*x1*x4 + k22*t1*x2*x3 - 2*k22*t2*x1*x3 + k22*t3*x1*x2 - k22*t1*x3*x4 - k22*t3*x1*x4 - k22*t3*x2*x3 + 2*k22*t4*x1*x3 + k22*t3*x3*x4 + k12*t1*x2*y1 - k12*t2*x1*y1 - k12*t1*x2*y3 - k12*t1*x4*y1 + k12*t2*x1*y3 + k12*t2*x3*y1 - k12*t3*x2*y1 + k12*t4*x1*y1 + k12*t1*x4*y3 - k12*t2*x3*y3 + k12*t3*x2*y3 + k12*t3*x4*y1 - k12*t4*x1*y3 - k12*t4*x3*y1 - k12*t3*x4*y3 + k12*t4*x3*y3 + k21*t1*x1*y2 - k21*t2*x1*y1 - k21*t1*x1*y4 - k21*t1*x3*y2 + k21*t2*x1*y3 + k21*t2*x3*y1 - k21*t3*x1*y2 + k21*t4*x1*y1 + k21*t1*x3*y4 - k21*t2*x3*y3 + k21*t3*x1*y4 + k21*t3*x3*y2 - k21*t4*x1*y3 - k21*t4*x3*y1 - k21*t3*x3*y4 + k21*t4*x3*y3 - k11*t1*y1*y2 + k11*t1*y1*y4 + k11*t1*y2*y3 - 2*k11*t2*y1*y3 + k11*t3*y1*y2 - k11*t1*y3*y4 - k11*t3*y1*y4 - k11*t3*y2*y3 + 2*k11*t4*y1*y3 + k11*t3*y3*y4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));

        }//end if element = quadrangle

    }//end loop over the elements


}


/*------------------------------------------------------TANGENT STIFFNESS MATRIX ROUTINE---------------------------------------------------*/
void Tangent_Stiffness_Matrix(std::vector<double> &theta_k, std::map<int, Parameter*> &region, std::vector<Element*> &elements,
                                gmm::row_matrix< gmm::wsvector<double> > &KT, map<Node*, Node*> &NodesCorresp,
                                std::vector<Node*> &nodes, FemFlag thermalOrElectrical)
{
    gmm::row_matrix< gmm::wsvector<double> > Tmp(nodes.size(), nodes.size());//Temporaray matrix in the building of KT
    gmm::dense_matrix<double> alpha(2,2); // matrice alpha
    gmm::dense_matrix<double> beta(2,2); // matrice beta pour la conductivité

    for(unsigned int i = 0; i < elements.size(); i++)//loop over the elements
    {
        if(elements[i]->type == 2)//If triangle
        {
            gmm::dense_matrix<double> Ke_1(3,3); //First term in KT
            gmm::dense_matrix<double> Ke_2(3,3); //Second term in KT

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
            }
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
            }

            std::vector<double> theta_nodes(3);//Temperature at nodes
            theta_nodes[0] = theta_k[n1->num-1];
            theta_nodes[1] = theta_k[n2->num-1];
            theta_nodes[2] = theta_k[n3->num-1];
            double thetam = (theta_nodes[0]+theta_nodes[1]+theta_nodes[2])/3;//Average temperature on the element
            double detJ = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);


            //Assembly of the conductivity tensor on the element
            gmm::dense_matrix<double> kappa(2,2);

            for (unsigned int k = 0; k<2;k++)
            {
                for(unsigned int j = 0 ; j<2; j++)
                {
                    kappa(k,j) = alpha(k,j)*thetam+beta(k,j);
                }
            }

            //Inverse of the Jacobian Matrix
            gmm::dense_matrix<double> inv_J(2,2);
            inv_J(0,0) = 1.0/detJ*(y3-y1);
            inv_J(0,1) = 1.0/detJ*(y1-y2);
            inv_J(1,0) = 1.0/detJ*(x1-x3);
            inv_J(1,1) = 1.0/detJ*(x2-x1);

            //Matrix containing the gradients of the shape functions in isoparametric coordinates
            gmm::dense_matrix<double> grad_phi(2,3);
            grad_phi(0,0) = -1;
            grad_phi(1,0) = -1;
            grad_phi(0,1) = 1;
            grad_phi(1,1) = 0;
            grad_phi(0,2) = 0;
            grad_phi(1,2) = 1;




            std::vector<double> tmp_vec(2,0);//Temporary vector in the computation of Matrix KT
            std::vector<double> tmp_vec_2(2,0);//Temporary vector in the computation of Matrix KT


            //First term of the elementary stiffness matrix : Ke_1

            for(unsigned int j = 0; j<3 ; j++)
            {
                tmp_vec[0] += theta_nodes[j]*grad_phi(0,j);
                tmp_vec[1] += theta_nodes[j]*grad_phi(1,j);
            }

            for(unsigned int j = 0; j<2; j++)
            {
                for(unsigned int k = 0; k<2 ; k++)
                {
                    tmp_vec_2[j] += inv_J(j,k)*tmp_vec[k];
                }
            }


            for(unsigned int j = 0; j<3; j++)//one contribution per shape function
            {
                tmp_vec[0]=0;
                tmp_vec[1]=0;
                for(unsigned int k = 0; k<2 ; k++)
                {
                    for(unsigned int l = 0 ; l<2; l++)
                    {
                        tmp_vec[k] += inv_J(k,l)*grad_phi(l,j);
                    }
                }
                for (unsigned int k =0;k<2;k++)
                {
                    for(unsigned int l = 0 ; l<2; l++)
                    {
                        Ke_1(j,0) += detJ/6*alpha(k,l)*tmp_vec[k]*tmp_vec_2[l];
                        Ke_1(j,1) += detJ/6*alpha(k,l)*tmp_vec[k]*tmp_vec_2[l];
                        Ke_1(j,2) += detJ/6*alpha(k,l)*tmp_vec[k]*tmp_vec_2[l];
                    }
                }
            }



            //Second term of the elementary stiffness matrix : Ke_2
             gmm::dense_matrix<double> tmp_mat(2,3);

             for(unsigned int j = 0; j<3; j++)//one contribution per shape function
            {

                for(unsigned int k = 0; k<2 ; k++)
                {
                    for(unsigned int l = 0 ; l<2; l++)
                    {
                        tmp_mat(k,j) += inv_J(k,l)*grad_phi(l,j);
                    }
                }
            }

            //gmm::mult(inv_J,grad_phi,tmp_mat)  ;
             for(unsigned int j = 0; j<3; j++)//one contribution per shape function
            {
                for(unsigned int k = 0; k<3 ; k++)
                {
                    for(unsigned int l = 0; l<2 ; l++)
                    {
                        for(unsigned int m = 0; m<2 ; m++)
                        {
                            Ke_2(j,k) += detJ/2*kappa(l,m)*tmp_mat(m,j)*tmp_mat(l,k);
                        }
                    }
                }
            }


            /*Utilisation du mapping NodesCorresp.
            Si une valeur doit être ajoutée à la ligne d'un noeud de droite, celle-ci est directement ajoutée à la ligne correspondant au noeud de gauche en vis-a-vis*/
            int num1 = NodesCorresp[n1]->num-1;
            int num2 = NodesCorresp[n2]->num-1;
            int num3 = NodesCorresp[n3]->num-1;

            Tmp(num1, n1->num-1) += (Ke_1(0,0) + Ke_2(0,0));
            Tmp(num1, n2->num-1) += (Ke_1(0,1) + Ke_2(0,1));
            Tmp(num1, n3->num-1) += (Ke_1(0,2) + Ke_2(0,2));
            Tmp(num2, n1->num-1) += (Ke_1(1,0) + Ke_2(1,0));
            Tmp(num2, n2->num-1) += (Ke_1(1,1) + Ke_2(1,1));
            Tmp(num2, n3->num-1) += (Ke_1(1,2) + Ke_2(1,2));
            Tmp(num3, n1->num-1) += (Ke_1(2,0) + Ke_2(2,0));
            Tmp(num3, n2->num-1) += (Ke_1(2,1) + Ke_2(2,1));
            Tmp(num3, n3->num-1) += (Ke_1(2,2) + Ke_2(2,2));

            gmm::clear(Ke_1);
            gmm::clear(Ke_2);
            gmm::clear(tmp_mat);
        }
        else if(elements[i]->type == 3)//If quad
        {
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
            }
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
            }

            std::vector<double> theta_nodes(4);//Temperature at nodes
            theta_nodes[0] = theta_k[n1->num-1];
            theta_nodes[1] = theta_k[n2->num-1];
            theta_nodes[2] = theta_k[n3->num-1];
            theta_nodes[3] = theta_k[n4->num-1];

            //Average temperature on the element
            double thetam = (theta_nodes[1]*((x1*y2)/4 - (x2*y1)/4 - (x1*y3)/4 + (x3*y1)/4 + (x2*y3)/4 - (x3*y2)/4)
                             + theta_nodes[0]*((x1*y2)/4 - (x2*y1)/4 - (x1*y4)/4 + (x4*y1)/4 + (x2*y4)/4 - (x4*y2)/4)
                             + theta_nodes[2]*((x1*y3)/4 - (x3*y1)/4 - (x1*y4)/4 + (x4*y1)/4 + (x3*y4)/4 - (x4*y3)/4)
                             + theta_nodes[3]*((x2*y3)/4 - (x3*y2)/4 - (x2*y4)/4 + (x4*y2)/4 + (x3*y4)/4 - (x4*y3)/4))
                             /((x1*y2)/2 - (x2*y1)/2 - (x1*y4)/2 + (x2*y3)/2 - (x3*y2)/2 + (x4*y1)/2 + (x3*y4)/2 - (x4*y3)/2);


            //Assembly of the conductivity tensor on the element
            gmm::dense_matrix<double> kappa(2,2);

            for (unsigned int k = 0; k<2;k++)
            {
                for(unsigned int j = 0 ; j<2; j++)
                {
                    kappa(k,j) = alpha(k,j)*thetam+beta(k,j);
                }
            }

            gmm::dense_matrix<double> Ke_1(4,4); //First term in KT
            gmm::dense_matrix<double> Ke_2(4,4); //Second term in KT

            double a11 = alpha(0,0);
            double a12 = alpha(0,1);
            double a21 = alpha(1,0);
            double a22 = alpha(1,1);

            double k11 = kappa(0,0);
            double k12 = kappa(0,1);
            double k21 = kappa(1,0);
            double k22 = kappa(1,1);

            double t1 = theta_nodes[0];
            double t2 = theta_nodes[1];
            double t3 = theta_nodes[2];
            double t4 = theta_nodes[3];

            Ke_1(0,0) = (a22*t1*x2*x2 + a22*t1*x4*x4 - a22*t3*x2*x2 - a22*t3*x4*x4 + a11*t1*y2*y2 + a11*t1*y4*y4 - a11*t3*y2*y2 - a11*t3*y4*y4 - a22*t2*x1*x2 - 2*a22*t1*x2*x4 + a22*t2*x1*x4 + a22*t2*x2*x3 + a22*t4*x1*x2 - a22*t2*x3*x4 + 2*a22*t3*x2*x4 - a22*t4*x1*x4 - a22*t4*x2*x3 + a22*t4*x3*x4 - a12*t1*x2*y2 + a12*t2*x1*y2 + a12*t1*x2*y4 + a12*t1*x4*y2 - a12*t2*x1*y4 - a12*t2*x3*y2 + a12*t3*x2*y2 - a12*t4*x1*y2 - a12*t1*x4*y4 + a12*t2*x3*y4 - a12*t3*x2*y4 - a12*t3*x4*y2 + a12*t4*x1*y4 + a12*t4*x3*y2 + a12*t3*x4*y4 - a12*t4*x3*y4 - a21*t1*x2*y2 + a21*t2*x2*y1 + a21*t1*x2*y4 + a21*t1*x4*y2 - a21*t2*x2*y3 - a21*t2*x4*y1 + a21*t3*x2*y2 - a21*t4*x2*y1 - a21*t1*x4*y4 + a21*t2*x4*y3 - a21*t3*x2*y4 - a21*t3*x4*y2 + a21*t4*x2*y3 + a21*t4*x4*y1 + a21*t3*x4*y4 - a21*t4*x4*y3 - a11*t2*y1*y2 - 2*a11*t1*y2*y4 + a11*t2*y1*y4 + a11*t2*y2*y3 + a11*t4*y1*y2 - a11*t2*y3*y4 + 2*a11*t3*y2*y4 - a11*t4*y1*y4 - a11*t4*y2*y3 + a11*t4*y3*y4)/(8*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));
            Ke_1(0,1) = (a22*t2*x1*x1 + a22*t2*x3*x3 - a22*t4*x1*x1 - a22*t4*x3*x3 + a11*t2*y1*y1 + a11*t2*y3*y3 - a11*t4*y1*y1 - a11*t4*y3*y3 - a22*t1*x1*x2 + a22*t1*x1*x4 + a22*t1*x2*x3 - 2*a22*t2*x1*x3 + a22*t3*x1*x2 - a22*t1*x3*x4 - a22*t3*x1*x4 - a22*t3*x2*x3 + 2*a22*t4*x1*x3 + a22*t3*x3*x4 + a12*t1*x2*y1 - a12*t2*x1*y1 - a12*t1*x2*y3 - a12*t1*x4*y1 + a12*t2*x1*y3 + a12*t2*x3*y1 - a12*t3*x2*y1 + a12*t4*x1*y1 + a12*t1*x4*y3 - a12*t2*x3*y3 + a12*t3*x2*y3 + a12*t3*x4*y1 - a12*t4*x1*y3 - a12*t4*x3*y1 - a12*t3*x4*y3 + a12*t4*x3*y3 + a21*t1*x1*y2 - a21*t2*x1*y1 - a21*t1*x1*y4 - a21*t1*x3*y2 + a21*t2*x1*y3 + a21*t2*x3*y1 - a21*t3*x1*y2 + a21*t4*x1*y1 + a21*t1*x3*y4 - a21*t2*x3*y3 + a21*t3*x1*y4 + a21*t3*x3*y2 - a21*t4*x1*y3 - a21*t4*x3*y1 - a21*t3*x3*y4 + a21*t4*x3*y3 - a11*t1*y1*y2 + a11*t1*y1*y4 + a11*t1*y2*y3 - 2*a11*t2*y1*y3 + a11*t3*y1*y2 - a11*t1*y3*y4 - a11*t3*y1*y4 - a11*t3*y2*y3 + 2*a11*t4*y1*y3 + a11*t3*y3*y4)/(8*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));
            Ke_1(0,2) = -(a22*t1*x2*x2 + a22*t1*x4*x4 - a22*t3*x2*x2 - a22*t3*x4*x4 + a11*t1*y2*y2 + a11*t1*y4*y4 - a11*t3*y2*y2 - a11*t3*y4*y4 - a22*t2*x1*x2 - 2*a22*t1*x2*x4 + a22*t2*x1*x4 + a22*t2*x2*x3 + a22*t4*x1*x2 - a22*t2*x3*x4 + 2*a22*t3*x2*x4 - a22*t4*x1*x4 - a22*t4*x2*x3 + a22*t4*x3*x4 - a12*t1*x2*y2 + a12*t2*x1*y2 + a12*t1*x2*y4 + a12*t1*x4*y2 - a12*t2*x1*y4 - a12*t2*x3*y2 + a12*t3*x2*y2 - a12*t4*x1*y2 - a12*t1*x4*y4 + a12*t2*x3*y4 - a12*t3*x2*y4 - a12*t3*x4*y2 + a12*t4*x1*y4 + a12*t4*x3*y2 + a12*t3*x4*y4 - a12*t4*x3*y4 - a21*t1*x2*y2 + a21*t2*x2*y1 + a21*t1*x2*y4 + a21*t1*x4*y2 - a21*t2*x2*y3 - a21*t2*x4*y1 + a21*t3*x2*y2 - a21*t4*x2*y1 - a21*t1*x4*y4 + a21*t2*x4*y3 - a21*t3*x2*y4 - a21*t3*x4*y2 + a21*t4*x2*y3 + a21*t4*x4*y1 + a21*t3*x4*y4 - a21*t4*x4*y3 - a11*t2*y1*y2 - 2*a11*t1*y2*y4 + a11*t2*y1*y4 + a11*t2*y2*y3 + a11*t4*y1*y2 - a11*t2*y3*y4 + 2*a11*t3*y2*y4 - a11*t4*y1*y4 - a11*t4*y2*y3 + a11*t4*y3*y4)/(8*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));
            Ke_1(0,3) = -(a22*t2*x1*x1 + a22*t2*x3*x3 - a22*t4*x1*x1 - a22*t4*x3*x3 + a11*t2*y1*y1 + a11*t2*y3*y3 - a11*t4*y1*y1 - a11*t4*y3*y3 - a22*t1*x1*x2 + a22*t1*x1*x4 + a22*t1*x2*x3 - 2*a22*t2*x1*x3 + a22*t3*x1*x2 - a22*t1*x3*x4 - a22*t3*x1*x4 - a22*t3*x2*x3 + 2*a22*t4*x1*x3 + a22*t3*x3*x4 + a12*t1*x2*y1 - a12*t2*x1*y1 - a12*t1*x2*y3 - a12*t1*x4*y1 + a12*t2*x1*y3 + a12*t2*x3*y1 - a12*t3*x2*y1 + a12*t4*x1*y1 + a12*t1*x4*y3 - a12*t2*x3*y3 + a12*t3*x2*y3 + a12*t3*x4*y1 - a12*t4*x1*y3 - a12*t4*x3*y1 - a12*t3*x4*y3 + a12*t4*x3*y3 + a21*t1*x1*y2 - a21*t2*x1*y1 - a21*t1*x1*y4 - a21*t1*x3*y2 + a21*t2*x1*y3 + a21*t2*x3*y1 - a21*t3*x1*y2 + a21*t4*x1*y1 + a21*t1*x3*y4 - a21*t2*x3*y3 + a21*t3*x1*y4 + a21*t3*x3*y2 - a21*t4*x1*y3 - a21*t4*x3*y1 - a21*t3*x3*y4 + a21*t4*x3*y3 - a11*t1*y1*y2 + a11*t1*y1*y4 + a11*t1*y2*y3 - 2*a11*t2*y1*y3 + a11*t3*y1*y2 - a11*t1*y3*y4 - a11*t3*y1*y4 - a11*t3*y2*y3 + 2*a11*t4*y1*y3 + a11*t3*y3*y4)/(8*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));

            Ke_1(1,0) = (a22*t1*x2*x2 + a22*t1*x4*x4 - a22*t3*x2*x2 - a22*t3*x4*x4 + a11*t1*y2*y2 + a11*t1*y4*y4 - a11*t3*y2*y2 - a11*t3*y4*y4 - a22*t2*x1*x2 - 2*a22*t1*x2*x4 + a22*t2*x1*x4 + a22*t2*x2*x3 + a22*t4*x1*x2 - a22*t2*x3*x4 + 2*a22*t3*x2*x4 - a22*t4*x1*x4 - a22*t4*x2*x3 + a22*t4*x3*x4 - a12*t1*x2*y2 + a12*t2*x1*y2 + a12*t1*x2*y4 + a12*t1*x4*y2 - a12*t2*x1*y4 - a12*t2*x3*y2 + a12*t3*x2*y2 - a12*t4*x1*y2 - a12*t1*x4*y4 + a12*t2*x3*y4 - a12*t3*x2*y4 - a12*t3*x4*y2 + a12*t4*x1*y4 + a12*t4*x3*y2 + a12*t3*x4*y4 - a12*t4*x3*y4 - a21*t1*x2*y2 + a21*t2*x2*y1 + a21*t1*x2*y4 + a21*t1*x4*y2 - a21*t2*x2*y3 - a21*t2*x4*y1 + a21*t3*x2*y2 - a21*t4*x2*y1 - a21*t1*x4*y4 + a21*t2*x4*y3 - a21*t3*x2*y4 - a21*t3*x4*y2 + a21*t4*x2*y3 + a21*t4*x4*y1 + a21*t3*x4*y4 - a21*t4*x4*y3 - a11*t2*y1*y2 - 2*a11*t1*y2*y4 + a11*t2*y1*y4 + a11*t2*y2*y3 + a11*t4*y1*y2 - a11*t2*y3*y4 + 2*a11*t3*y2*y4 - a11*t4*y1*y4 - a11*t4*y2*y3 + a11*t4*y3*y4)/(8*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));
            Ke_1(1,1) = (a22*t2*x1*x1 + a22*t2*x3*x3 - a22*t4*x1*x1 - a22*t4*x3*x3 + a11*t2*y1*y1 + a11*t2*y3*y3 - a11*t4*y1*y1 - a11*t4*y3*y3 - a22*t1*x1*x2 + a22*t1*x1*x4 + a22*t1*x2*x3 - 2*a22*t2*x1*x3 + a22*t3*x1*x2 - a22*t1*x3*x4 - a22*t3*x1*x4 - a22*t3*x2*x3 + 2*a22*t4*x1*x3 + a22*t3*x3*x4 + a12*t1*x2*y1 - a12*t2*x1*y1 - a12*t1*x2*y3 - a12*t1*x4*y1 + a12*t2*x1*y3 + a12*t2*x3*y1 - a12*t3*x2*y1 + a12*t4*x1*y1 + a12*t1*x4*y3 - a12*t2*x3*y3 + a12*t3*x2*y3 + a12*t3*x4*y1 - a12*t4*x1*y3 - a12*t4*x3*y1 - a12*t3*x4*y3 + a12*t4*x3*y3 + a21*t1*x1*y2 - a21*t2*x1*y1 - a21*t1*x1*y4 - a21*t1*x3*y2 + a21*t2*x1*y3 + a21*t2*x3*y1 - a21*t3*x1*y2 + a21*t4*x1*y1 + a21*t1*x3*y4 - a21*t2*x3*y3 + a21*t3*x1*y4 + a21*t3*x3*y2 - a21*t4*x1*y3 - a21*t4*x3*y1 - a21*t3*x3*y4 + a21*t4*x3*y3 - a11*t1*y1*y2 + a11*t1*y1*y4 + a11*t1*y2*y3 - 2*a11*t2*y1*y3 + a11*t3*y1*y2 - a11*t1*y3*y4 - a11*t3*y1*y4 - a11*t3*y2*y3 + 2*a11*t4*y1*y3 + a11*t3*y3*y4)/(8*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));
            Ke_1(1,2) = -(a22*t1*x2*x2 + a22*t1*x4*x4 - a22*t3*x2*x2 - a22*t3*x4*x4 + a11*t1*y2*y2 + a11*t1*y4*y4 - a11*t3*y2*y2 - a11*t3*y4*y4 - a22*t2*x1*x2 - 2*a22*t1*x2*x4 + a22*t2*x1*x4 + a22*t2*x2*x3 + a22*t4*x1*x2 - a22*t2*x3*x4 + 2*a22*t3*x2*x4 - a22*t4*x1*x4 - a22*t4*x2*x3 + a22*t4*x3*x4 - a12*t1*x2*y2 + a12*t2*x1*y2 + a12*t1*x2*y4 + a12*t1*x4*y2 - a12*t2*x1*y4 - a12*t2*x3*y2 + a12*t3*x2*y2 - a12*t4*x1*y2 - a12*t1*x4*y4 + a12*t2*x3*y4 - a12*t3*x2*y4 - a12*t3*x4*y2 + a12*t4*x1*y4 + a12*t4*x3*y2 + a12*t3*x4*y4 - a12*t4*x3*y4 - a21*t1*x2*y2 + a21*t2*x2*y1 + a21*t1*x2*y4 + a21*t1*x4*y2 - a21*t2*x2*y3 - a21*t2*x4*y1 + a21*t3*x2*y2 - a21*t4*x2*y1 - a21*t1*x4*y4 + a21*t2*x4*y3 - a21*t3*x2*y4 - a21*t3*x4*y2 + a21*t4*x2*y3 + a21*t4*x4*y1 + a21*t3*x4*y4 - a21*t4*x4*y3 - a11*t2*y1*y2 - 2*a11*t1*y2*y4 + a11*t2*y1*y4 + a11*t2*y2*y3 + a11*t4*y1*y2 - a11*t2*y3*y4 + 2*a11*t3*y2*y4 - a11*t4*y1*y4 - a11*t4*y2*y3 + a11*t4*y3*y4)/(8*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));
            Ke_1(1,3) = -(a22*t2*x1*x1 + a22*t2*x3*x3 - a22*t4*x1*x1 - a22*t4*x3*x3 + a11*t2*y1*y1 + a11*t2*y3*y3 - a11*t4*y1*y1 - a11*t4*y3*y3 - a22*t1*x1*x2 + a22*t1*x1*x4 + a22*t1*x2*x3 - 2*a22*t2*x1*x3 + a22*t3*x1*x2 - a22*t1*x3*x4 - a22*t3*x1*x4 - a22*t3*x2*x3 + 2*a22*t4*x1*x3 + a22*t3*x3*x4 + a12*t1*x2*y1 - a12*t2*x1*y1 - a12*t1*x2*y3 - a12*t1*x4*y1 + a12*t2*x1*y3 + a12*t2*x3*y1 - a12*t3*x2*y1 + a12*t4*x1*y1 + a12*t1*x4*y3 - a12*t2*x3*y3 + a12*t3*x2*y3 + a12*t3*x4*y1 - a12*t4*x1*y3 - a12*t4*x3*y1 - a12*t3*x4*y3 + a12*t4*x3*y3 + a21*t1*x1*y2 - a21*t2*x1*y1 - a21*t1*x1*y4 - a21*t1*x3*y2 + a21*t2*x1*y3 + a21*t2*x3*y1 - a21*t3*x1*y2 + a21*t4*x1*y1 + a21*t1*x3*y4 - a21*t2*x3*y3 + a21*t3*x1*y4 + a21*t3*x3*y2 - a21*t4*x1*y3 - a21*t4*x3*y1 - a21*t3*x3*y4 + a21*t4*x3*y3 - a11*t1*y1*y2 + a11*t1*y1*y4 + a11*t1*y2*y3 - 2*a11*t2*y1*y3 + a11*t3*y1*y2 - a11*t1*y3*y4 - a11*t3*y1*y4 - a11*t3*y2*y3 + 2*a11*t4*y1*y3 + a11*t3*y3*y4)/(8*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));

            Ke_1(2,0) = (a22*t1*x2*x2 + a22*t1*x4*x4 - a22*t3*x2*x2 - a22*t3*x4*x4 + a11*t1*y2*y2 + a11*t1*y4*y4 - a11*t3*y2*y2 - a11*t3*y4*y4 - a22*t2*x1*x2 - 2*a22*t1*x2*x4 + a22*t2*x1*x4 + a22*t2*x2*x3 + a22*t4*x1*x2 - a22*t2*x3*x4 + 2*a22*t3*x2*x4 - a22*t4*x1*x4 - a22*t4*x2*x3 + a22*t4*x3*x4 - a12*t1*x2*y2 + a12*t2*x1*y2 + a12*t1*x2*y4 + a12*t1*x4*y2 - a12*t2*x1*y4 - a12*t2*x3*y2 + a12*t3*x2*y2 - a12*t4*x1*y2 - a12*t1*x4*y4 + a12*t2*x3*y4 - a12*t3*x2*y4 - a12*t3*x4*y2 + a12*t4*x1*y4 + a12*t4*x3*y2 + a12*t3*x4*y4 - a12*t4*x3*y4 - a21*t1*x2*y2 + a21*t2*x2*y1 + a21*t1*x2*y4 + a21*t1*x4*y2 - a21*t2*x2*y3 - a21*t2*x4*y1 + a21*t3*x2*y2 - a21*t4*x2*y1 - a21*t1*x4*y4 + a21*t2*x4*y3 - a21*t3*x2*y4 - a21*t3*x4*y2 + a21*t4*x2*y3 + a21*t4*x4*y1 + a21*t3*x4*y4 - a21*t4*x4*y3 - a11*t2*y1*y2 - 2*a11*t1*y2*y4 + a11*t2*y1*y4 + a11*t2*y2*y3 + a11*t4*y1*y2 - a11*t2*y3*y4 + 2*a11*t3*y2*y4 - a11*t4*y1*y4 - a11*t4*y2*y3 + a11*t4*y3*y4)/(8*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));
            Ke_1(2,1) = (a22*t2*x1*x1 + a22*t2*x3*x3 - a22*t4*x1*x1 - a22*t4*x3*x3 + a11*t2*y1*y1 + a11*t2*y3*y3 - a11*t4*y1*y1 - a11*t4*y3*y3 - a22*t1*x1*x2 + a22*t1*x1*x4 + a22*t1*x2*x3 - 2*a22*t2*x1*x3 + a22*t3*x1*x2 - a22*t1*x3*x4 - a22*t3*x1*x4 - a22*t3*x2*x3 + 2*a22*t4*x1*x3 + a22*t3*x3*x4 + a12*t1*x2*y1 - a12*t2*x1*y1 - a12*t1*x2*y3 - a12*t1*x4*y1 + a12*t2*x1*y3 + a12*t2*x3*y1 - a12*t3*x2*y1 + a12*t4*x1*y1 + a12*t1*x4*y3 - a12*t2*x3*y3 + a12*t3*x2*y3 + a12*t3*x4*y1 - a12*t4*x1*y3 - a12*t4*x3*y1 - a12*t3*x4*y3 + a12*t4*x3*y3 + a21*t1*x1*y2 - a21*t2*x1*y1 - a21*t1*x1*y4 - a21*t1*x3*y2 + a21*t2*x1*y3 + a21*t2*x3*y1 - a21*t3*x1*y2 + a21*t4*x1*y1 + a21*t1*x3*y4 - a21*t2*x3*y3 + a21*t3*x1*y4 + a21*t3*x3*y2 - a21*t4*x1*y3 - a21*t4*x3*y1 - a21*t3*x3*y4 + a21*t4*x3*y3 - a11*t1*y1*y2 + a11*t1*y1*y4 + a11*t1*y2*y3 - 2*a11*t2*y1*y3 + a11*t3*y1*y2 - a11*t1*y3*y4 - a11*t3*y1*y4 - a11*t3*y2*y3 + 2*a11*t4*y1*y3 + a11*t3*y3*y4)/(8*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));
            Ke_1(2,2) = -(a22*t1*x2*x2 + a22*t1*x4*x4 - a22*t3*x2*x2 - a22*t3*x4*x4 + a11*t1*y2*y2 + a11*t1*y4*y4 - a11*t3*y2*y2 - a11*t3*y4*y4 - a22*t2*x1*x2 - 2*a22*t1*x2*x4 + a22*t2*x1*x4 + a22*t2*x2*x3 + a22*t4*x1*x2 - a22*t2*x3*x4 + 2*a22*t3*x2*x4 - a22*t4*x1*x4 - a22*t4*x2*x3 + a22*t4*x3*x4 - a12*t1*x2*y2 + a12*t2*x1*y2 + a12*t1*x2*y4 + a12*t1*x4*y2 - a12*t2*x1*y4 - a12*t2*x3*y2 + a12*t3*x2*y2 - a12*t4*x1*y2 - a12*t1*x4*y4 + a12*t2*x3*y4 - a12*t3*x2*y4 - a12*t3*x4*y2 + a12*t4*x1*y4 + a12*t4*x3*y2 + a12*t3*x4*y4 - a12*t4*x3*y4 - a21*t1*x2*y2 + a21*t2*x2*y1 + a21*t1*x2*y4 + a21*t1*x4*y2 - a21*t2*x2*y3 - a21*t2*x4*y1 + a21*t3*x2*y2 - a21*t4*x2*y1 - a21*t1*x4*y4 + a21*t2*x4*y3 - a21*t3*x2*y4 - a21*t3*x4*y2 + a21*t4*x2*y3 + a21*t4*x4*y1 + a21*t3*x4*y4 - a21*t4*x4*y3 - a11*t2*y1*y2 - 2*a11*t1*y2*y4 + a11*t2*y1*y4 + a11*t2*y2*y3 + a11*t4*y1*y2 - a11*t2*y3*y4 + 2*a11*t3*y2*y4 - a11*t4*y1*y4 - a11*t4*y2*y3 + a11*t4*y3*y4)/(8*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));
            Ke_1(2,3) = -(a22*t2*x1*x1 + a22*t2*x3*x3 - a22*t4*x1*x1 - a22*t4*x3*x3 + a11*t2*y1*y1 + a11*t2*y3*y3 - a11*t4*y1*y1 - a11*t4*y3*y3 - a22*t1*x1*x2 + a22*t1*x1*x4 + a22*t1*x2*x3 - 2*a22*t2*x1*x3 + a22*t3*x1*x2 - a22*t1*x3*x4 - a22*t3*x1*x4 - a22*t3*x2*x3 + 2*a22*t4*x1*x3 + a22*t3*x3*x4 + a12*t1*x2*y1 - a12*t2*x1*y1 - a12*t1*x2*y3 - a12*t1*x4*y1 + a12*t2*x1*y3 + a12*t2*x3*y1 - a12*t3*x2*y1 + a12*t4*x1*y1 + a12*t1*x4*y3 - a12*t2*x3*y3 + a12*t3*x2*y3 + a12*t3*x4*y1 - a12*t4*x1*y3 - a12*t4*x3*y1 - a12*t3*x4*y3 + a12*t4*x3*y3 + a21*t1*x1*y2 - a21*t2*x1*y1 - a21*t1*x1*y4 - a21*t1*x3*y2 + a21*t2*x1*y3 + a21*t2*x3*y1 - a21*t3*x1*y2 + a21*t4*x1*y1 + a21*t1*x3*y4 - a21*t2*x3*y3 + a21*t3*x1*y4 + a21*t3*x3*y2 - a21*t4*x1*y3 - a21*t4*x3*y1 - a21*t3*x3*y4 + a21*t4*x3*y3 - a11*t1*y1*y2 + a11*t1*y1*y4 + a11*t1*y2*y3 - 2*a11*t2*y1*y3 + a11*t3*y1*y2 - a11*t1*y3*y4 - a11*t3*y1*y4 - a11*t3*y2*y3 + 2*a11*t4*y1*y3 + a11*t3*y3*y4)/(8*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));

            Ke_1(3,0) = (a22*t1*x2*x2 + a22*t1*x4*x4 - a22*t3*x2*x2 - a22*t3*x4*x4 + a11*t1*y2*y2 + a11*t1*y4*y4 - a11*t3*y2*y2 - a11*t3*y4*y4 - a22*t2*x1*x2 - 2*a22*t1*x2*x4 + a22*t2*x1*x4 + a22*t2*x2*x3 + a22*t4*x1*x2 - a22*t2*x3*x4 + 2*a22*t3*x2*x4 - a22*t4*x1*x4 - a22*t4*x2*x3 + a22*t4*x3*x4 - a12*t1*x2*y2 + a12*t2*x1*y2 + a12*t1*x2*y4 + a12*t1*x4*y2 - a12*t2*x1*y4 - a12*t2*x3*y2 + a12*t3*x2*y2 - a12*t4*x1*y2 - a12*t1*x4*y4 + a12*t2*x3*y4 - a12*t3*x2*y4 - a12*t3*x4*y2 + a12*t4*x1*y4 + a12*t4*x3*y2 + a12*t3*x4*y4 - a12*t4*x3*y4 - a21*t1*x2*y2 + a21*t2*x2*y1 + a21*t1*x2*y4 + a21*t1*x4*y2 - a21*t2*x2*y3 - a21*t2*x4*y1 + a21*t3*x2*y2 - a21*t4*x2*y1 - a21*t1*x4*y4 + a21*t2*x4*y3 - a21*t3*x2*y4 - a21*t3*x4*y2 + a21*t4*x2*y3 + a21*t4*x4*y1 + a21*t3*x4*y4 - a21*t4*x4*y3 - a11*t2*y1*y2 - 2*a11*t1*y2*y4 + a11*t2*y1*y4 + a11*t2*y2*y3 + a11*t4*y1*y2 - a11*t2*y3*y4 + 2*a11*t3*y2*y4 - a11*t4*y1*y4 - a11*t4*y2*y3 + a11*t4*y3*y4)/(8*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));
            Ke_1(3,1) = (a22*t2*x1*x1 + a22*t2*x3*x3 - a22*t4*x1*x1 - a22*t4*x3*x3 + a11*t2*y1*y1 + a11*t2*y3*y3 - a11*t4*y1*y1 - a11*t4*y3*y3 - a22*t1*x1*x2 + a22*t1*x1*x4 + a22*t1*x2*x3 - 2*a22*t2*x1*x3 + a22*t3*x1*x2 - a22*t1*x3*x4 - a22*t3*x1*x4 - a22*t3*x2*x3 + 2*a22*t4*x1*x3 + a22*t3*x3*x4 + a12*t1*x2*y1 - a12*t2*x1*y1 - a12*t1*x2*y3 - a12*t1*x4*y1 + a12*t2*x1*y3 + a12*t2*x3*y1 - a12*t3*x2*y1 + a12*t4*x1*y1 + a12*t1*x4*y3 - a12*t2*x3*y3 + a12*t3*x2*y3 + a12*t3*x4*y1 - a12*t4*x1*y3 - a12*t4*x3*y1 - a12*t3*x4*y3 + a12*t4*x3*y3 + a21*t1*x1*y2 - a21*t2*x1*y1 - a21*t1*x1*y4 - a21*t1*x3*y2 + a21*t2*x1*y3 + a21*t2*x3*y1 - a21*t3*x1*y2 + a21*t4*x1*y1 + a21*t1*x3*y4 - a21*t2*x3*y3 + a21*t3*x1*y4 + a21*t3*x3*y2 - a21*t4*x1*y3 - a21*t4*x3*y1 - a21*t3*x3*y4 + a21*t4*x3*y3 - a11*t1*y1*y2 + a11*t1*y1*y4 + a11*t1*y2*y3 - 2*a11*t2*y1*y3 + a11*t3*y1*y2 - a11*t1*y3*y4 - a11*t3*y1*y4 - a11*t3*y2*y3 + 2*a11*t4*y1*y3 + a11*t3*y3*y4)/(8*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));
            Ke_1(3,2) = -(a22*t1*x2*x2 + a22*t1*x4*x4 - a22*t3*x2*x2 - a22*t3*x4*x4 + a11*t1*y2*y2 + a11*t1*y4*y4 - a11*t3*y2*y2 - a11*t3*y4*y4 - a22*t2*x1*x2 - 2*a22*t1*x2*x4 + a22*t2*x1*x4 + a22*t2*x2*x3 + a22*t4*x1*x2 - a22*t2*x3*x4 + 2*a22*t3*x2*x4 - a22*t4*x1*x4 - a22*t4*x2*x3 + a22*t4*x3*x4 - a12*t1*x2*y2 + a12*t2*x1*y2 + a12*t1*x2*y4 + a12*t1*x4*y2 - a12*t2*x1*y4 - a12*t2*x3*y2 + a12*t3*x2*y2 - a12*t4*x1*y2 - a12*t1*x4*y4 + a12*t2*x3*y4 - a12*t3*x2*y4 - a12*t3*x4*y2 + a12*t4*x1*y4 + a12*t4*x3*y2 + a12*t3*x4*y4 - a12*t4*x3*y4 - a21*t1*x2*y2 + a21*t2*x2*y1 + a21*t1*x2*y4 + a21*t1*x4*y2 - a21*t2*x2*y3 - a21*t2*x4*y1 + a21*t3*x2*y2 - a21*t4*x2*y1 - a21*t1*x4*y4 + a21*t2*x4*y3 - a21*t3*x2*y4 - a21*t3*x4*y2 + a21*t4*x2*y3 + a21*t4*x4*y1 + a21*t3*x4*y4 - a21*t4*x4*y3 - a11*t2*y1*y2 - 2*a11*t1*y2*y4 + a11*t2*y1*y4 + a11*t2*y2*y3 + a11*t4*y1*y2 - a11*t2*y3*y4 + 2*a11*t3*y2*y4 - a11*t4*y1*y4 - a11*t4*y2*y3 + a11*t4*y3*y4)/(8*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));
            Ke_1(3,3) = -(a22*t2*x1*x1 + a22*t2*x3*x3 - a22*t4*x1*x1 - a22*t4*x3*x3 + a11*t2*y1*y1 + a11*t2*y3*y3 - a11*t4*y1*y1 - a11*t4*y3*y3 - a22*t1*x1*x2 + a22*t1*x1*x4 + a22*t1*x2*x3 - 2*a22*t2*x1*x3 + a22*t3*x1*x2 - a22*t1*x3*x4 - a22*t3*x1*x4 - a22*t3*x2*x3 + 2*a22*t4*x1*x3 + a22*t3*x3*x4 + a12*t1*x2*y1 - a12*t2*x1*y1 - a12*t1*x2*y3 - a12*t1*x4*y1 + a12*t2*x1*y3 + a12*t2*x3*y1 - a12*t3*x2*y1 + a12*t4*x1*y1 + a12*t1*x4*y3 - a12*t2*x3*y3 + a12*t3*x2*y3 + a12*t3*x4*y1 - a12*t4*x1*y3 - a12*t4*x3*y1 - a12*t3*x4*y3 + a12*t4*x3*y3 + a21*t1*x1*y2 - a21*t2*x1*y1 - a21*t1*x1*y4 - a21*t1*x3*y2 + a21*t2*x1*y3 + a21*t2*x3*y1 - a21*t3*x1*y2 + a21*t4*x1*y1 + a21*t1*x3*y4 - a21*t2*x3*y3 + a21*t3*x1*y4 + a21*t3*x3*y2 - a21*t4*x1*y3 - a21*t4*x3*y1 - a21*t3*x3*y4 + a21*t4*x3*y3 - a11*t1*y1*y2 + a11*t1*y1*y4 + a11*t1*y2*y3 - 2*a11*t2*y1*y3 + a11*t3*y1*y2 - a11*t1*y3*y4 - a11*t3*y1*y4 - a11*t3*y2*y3 + 2*a11*t4*y1*y3 + a11*t3*y3*y4)/(8*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));


            Ke_2(0,0) = (k22*x2*x2 + k22*x4*x4 + k11*y2*y2 + k11*y4*y4 - 2*k22*x2*x4 - k12*x2*y2 + k12*x2*y4 + k12*x4*y2 - k12*x4*y4 - k21*x2*y2 + k21*x2*y4 + k21*x4*y2 - k21*x4*y4 - 2*k11*y2*y4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));
            Ke_2(0,1) = -(k22*x1*x2 - k22*x1*x4 - k22*x2*x3 + k22*x3*x4 - k12*x2*y1 + k12*x2*y3 + k12*x4*y1 - k12*x4*y3 - k21*x1*y2 + k21*x1*y4 + k21*x3*y2 - k21*x3*y4 + k11*y1*y2 - k11*y1*y4 - k11*y2*y3 + k11*y3*y4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));
            Ke_2(0,2) = -(k22*x2*x2 + k22*x4*x4 + k11*y2*y2 + k11*y4*y4 - 2*k22*x2*x4 - k12*x2*y2 + k12*x2*y4 + k12*x4*y2 - k12*x4*y4 - k21*x2*y2 + k21*x2*y4 + k21*x4*y2 - k21*x4*y4 - 2*k11*y2*y4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));
            Ke_2(0,3) = (k22*x1*x2 - k22*x1*x4 - k22*x2*x3 + k22*x3*x4 - k12*x2*y1 + k12*x2*y3 + k12*x4*y1 - k12*x4*y3 - k21*x1*y2 + k21*x1*y4 + k21*x3*y2 - k21*x3*y4 + k11*y1*y2 - k11*y1*y4 - k11*y2*y3 + k11*y3*y4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));

            Ke_2(1,0) = -(k22*x1*x2 - k22*x1*x4 - k22*x2*x3 + k22*x3*x4 - k12*x1*y2 + k12*x1*y4 + k12*x3*y2 - k12*x3*y4 - k21*x2*y1 + k21*x2*y3 + k21*x4*y1 - k21*x4*y3 + k11*y1*y2 - k11*y1*y4 - k11*y2*y3 + k11*y3*y4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));
            Ke_2(1,1) = (k22*x1*x1 + k22*x3*x3 + k11*y1*y1 + k11*y3*y3 - 2*k22*x1*x3 - k12*x1*y1 + k12*x1*y3 + k12*x3*y1 - k12*x3*y3 - k21*x1*y1 + k21*x1*y3 + k21*x3*y1 - k21*x3*y3 - 2*k11*y1*y3)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));
            Ke_2(1,2) = (k22*x1*x2 - k22*x1*x4 - k22*x2*x3 + k22*x3*x4 - k12*x1*y2 + k12*x1*y4 + k12*x3*y2 - k12*x3*y4 - k21*x2*y1 + k21*x2*y3 + k21*x4*y1 - k21*x4*y3 + k11*y1*y2 - k11*y1*y4 - k11*y2*y3 + k11*y3*y4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));
            Ke_2(1,3) = -(k22*x1*x1 + k22*x3*x3 + k11*y1*y1 + k11*y3*y3 - 2*k22*x1*x3 - k12*x1*y1 + k12*x1*y3 + k12*x3*y1 - k12*x3*y3 - k21*x1*y1 + k21*x1*y3 + k21*x3*y1 - k21*x3*y3 - 2*k11*y1*y3)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));

            Ke_2(2,0) = -(k22*x2*x2 + k22*x4*x4 + k11*y2*y2 + k11*y4*y4 - 2*k22*x2*x4 - k12*x2*y2 + k12*x2*y4 + k12*x4*y2 - k12*x4*y4 - k21*x2*y2 + k21*x2*y4 + k21*x4*y2 - k21*x4*y4 - 2*k11*y2*y4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));
            Ke_2(2,1) = (k22*x1*x2 - k22*x1*x4 - k22*x2*x3 + k22*x3*x4 - k12*x2*y1 + k12*x2*y3 + k12*x4*y1 - k12*x4*y3 - k21*x1*y2 + k21*x1*y4 + k21*x3*y2 - k21*x3*y4 + k11*y1*y2 - k11*y1*y4 - k11*y2*y3 + k11*y3*y4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));
            Ke_2(2,2) = (k22*x2*x2 + k22*x4*x4 + k11*y2*y2 + k11*y4*y4 - 2*k22*x2*x4 - k12*x2*y2 + k12*x2*y4 + k12*x4*y2 - k12*x4*y4 - k21*x2*y2 + k21*x2*y4 + k21*x4*y2 - k21*x4*y4 - 2*k11*y2*y4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));
            Ke_2(2,3) = -(k22*x1*x2 - k22*x1*x4 - k22*x2*x3 + k22*x3*x4 - k12*x2*y1 + k12*x2*y3 + k12*x4*y1 - k12*x4*y3 - k21*x1*y2 + k21*x1*y4 + k21*x3*y2 - k21*x3*y4 + k11*y1*y2 - k11*y1*y4 - k11*y2*y3 + k11*y3*y4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));

            Ke_2(3,0) = (k22*x1*x2 - k22*x1*x4 - k22*x2*x3 + k22*x3*x4 - k12*x1*y2 + k12*x1*y4 + k12*x3*y2 - k12*x3*y4 - k21*x2*y1 + k21*x2*y3 + k21*x4*y1 - k21*x4*y3 + k11*y1*y2 - k11*y1*y4 - k11*y2*y3 + k11*y3*y4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));
            Ke_2(3,1) = -(k22*x1*x1 + k22*x3*x3 + k11*y1*y1 + k11*y3*y3 - 2*k22*x1*x3 - k12*x1*y1 + k12*x1*y3 + k12*x3*y1 - k12*x3*y3 - k21*x1*y1 + k21*x1*y3 + k21*x3*y1 - k21*x3*y3 - 2*k11*y1*y3)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));
            Ke_2(3,2) = -(k22*x1*x2 - k22*x1*x4 - k22*x2*x3 + k22*x3*x4 - k12*x1*y2 + k12*x1*y4 + k12*x3*y2 - k12*x3*y4 - k21*x2*y1 + k21*x2*y3 + k21*x4*y1 - k21*x4*y3 + k11*y1*y2 - k11*y1*y4 - k11*y2*y3 + k11*y3*y4)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));
            Ke_2(3,3) = (k22*x1*x1 + k22*x3*x3 + k11*y1*y1 + k11*y3*y3 - 2*k22*x1*x3 - k12*x1*y1 + k12*x1*y3 + k12*x3*y1 - k12*x3*y3 - k21*x1*y1 + k21*x1*y3 + k21*x3*y1 - k21*x3*y3 - 2*k11*y1*y3)/(2*(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3));

            /*Utilisation du mapping NodesCorresp.
            Si une valeur doit être ajoutée à la ligne d'un noeud de droite, celle-ci est directement ajoutée à la ligne correspondant au noeud de gauche en vis-a-vis*/
            int num1 = NodesCorresp[n1]->num-1;
            int num2 = NodesCorresp[n2]->num-1;
            int num3 = NodesCorresp[n3]->num-1;
            int num4 = NodesCorresp[n4]->num-1;

            Tmp(num1, n1->num-1) += (Ke_1(0,0) + Ke_2(0,0));
            Tmp(num1, n2->num-1) += (Ke_1(0,1) + Ke_2(0,1));
            Tmp(num1, n3->num-1) += (Ke_1(0,2) + Ke_2(0,2));
            Tmp(num1, n4->num-1) += (Ke_1(0,3) + Ke_2(0,3));

            Tmp(num2, n1->num-1) += (Ke_1(1,0) + Ke_2(1,0));
            Tmp(num2, n2->num-1) += (Ke_1(1,1) + Ke_2(1,1));
            Tmp(num2, n3->num-1) += (Ke_1(1,2) + Ke_2(1,2));
            Tmp(num2, n4->num-1) += (Ke_1(1,3) + Ke_2(1,3));

            Tmp(num3, n1->num-1) += (Ke_1(2,0) + Ke_2(2,0));
            Tmp(num3, n2->num-1) += (Ke_1(2,1) + Ke_2(2,1));
            Tmp(num3, n3->num-1) += (Ke_1(2,2) + Ke_2(2,2));
            Tmp(num3, n4->num-1) += (Ke_1(2,3) + Ke_2(2,3));

            Tmp(num4, n1->num-1) += (Ke_1(3,0) + Ke_2(3,0));
            Tmp(num4, n2->num-1) += (Ke_1(3,1) + Ke_2(3,1));
            Tmp(num4, n3->num-1) += (Ke_1(3,2) + Ke_2(3,2));
            Tmp(num4, n4->num-1) += (Ke_1(3,3) + Ke_2(3,3));

        }//end if element = quadrangle
    }
    //cout << "In Tangent Stiffness Matrix, before copy" << endl;
    gmm::copy(Tmp,KT);
}//end of Tangent Stiffness Matrix routine


/*----------------------END CRITERION ROUTINE---------------------------*/
bool End_Criterion(std::vector<double> &RHS, double normRHS0 , double eps, Type type)
{
    //Threshold value
    double criterion = gmm::vect_norm2(RHS);
    if(normRHS0 > eps)
    {
         criterion = criterion/normRHS0;
    }
    //cout << criterion << endl;
    //cout << normRHS0 << " " << gmm::vect_norm2(RHS) << endl;
    if(type == DIRICHLET || type == VONNEUMANN || type == PERIODIC) cout << "Newton Raphson relative residue = " << gmm::vect_norm2(RHS)/normRHS0 <<endl;
    if(criterion >eps)
        return false;

    return true;
}

