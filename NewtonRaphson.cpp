#include <iostream>
#include <cstdio>
#include <vector>
#include <map>
#include <string>
#include "gmm/gmm.h"
#include "gmshio.h"
#include "physicalio.h"
#include "NewtonRaphson.h"
#include "fem.h"

using namespace std;

/*---------------------------------------------NEWTON RAPHSON ROUTINE--------------------------------------------------------------*/

void NewtonRaphson(std::vector<Node*> &nodes, std::vector<Element*> &elements, std::vector<Physical*> &physicals, std::vector<Parameter*> &parameters,
					std::map<Node*, std::vector<double> > &solution, std::vector<double> &theta_k,gmm::row_matrix< gmm::wsvector<double> > Tmp,
					std::vector<double> &qext,FemFlag method,std::map<int, double> linesRegion, std::map<Node*, Node*> &NodesCorresp,
					std::vector<double> &delta_theta_k,std::map<int, std::vector<double> > &surfaceRegion,std::vector<double> &RHS)
{

	std::vector<double> qint(nodes.size()); //Internal flux vector
	gmm::row_matrix< gmm::wsvector<double> > KT_tmp(nodes.size(), nodes.size());//Tangent matrix used for build-in
    gmm::csr_matrix<double> KT;//Tangent matrix used in system solving

    //Internal Flux : qint
	Internal_flux(theta_k, surfaceRegion, elements, qint);

    //TANGENT STIFFNESS MATRIX :  KT 
    Tangent_Stiffness_Matrix(theta_k,surfaceRegion,elements,
                             KT_tmp,NodesCorresp,nodes);
							 
    //Building of the Right Hand Side                      
    for (unsigned int i=0 ; i<nodes.size(); i++)
        RHS[i] = qint[i]+qext[i];

    //Including the Dirichlet boundary conditions
    if(method == DIRICHLETFLAG)
    {
        //Dirichlet conditions in KT_tmp
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

        //Including the Dirichlet boundary conditions in RHS
        for(unsigned int i = 0; i < elements.size(); i++)
        {
            if(elements[i]->type == 1)//If line
            {
                if(linesRegion.count(elements[i]->region) == 1)//If linesRegion contains elements[i]->region
                {
                    for(unsigned int j = 0; j < elements[i]->nodes.size(); j++)
                    {
                        RHS[elements[i]->nodes[j]->num-1] = 0;
                    }
                }
            }
        }
    }



    //Solving the system
    gmm::copy(KT_tmp, KT);
    gmm::lu_solve(KT, delta_theta_k, RHS);

    for (unsigned int i=0; i<nodes.size(); i++)
	{
		theta_k[i]+=delta_theta_k[i];
	}
       
	
}//end of Newton Raphson


/*----------------------------------------------------------INTERNAL FLUX ROUTINE----------------------------------------------------------------------*/
void Internal_flux(std::vector<double> &theta_k, std::map<int,std::vector<double> > &surfaceRegion, std::vector<Element*> &elements, std::vector<double> &qint)
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

            alpha(0,0) = surfaceRegion[elements[i]->region][1];
            alpha(0,1) = surfaceRegion[elements[i]->region][2];
            alpha(1,0) = surfaceRegion[elements[i]->region][3];
            alpha(1,1) = surfaceRegion[elements[i]->region][4];

            beta(0,0) = surfaceRegion[elements[i]->region][5];
            beta(0,1) = surfaceRegion[elements[i]->region][6];
            beta(1,0) = surfaceRegion[elements[i]->region][7];
            beta(1,1) = surfaceRegion[elements[i]->region][8];

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


    }//end loop over the elements


}


/*------------------------------------------------------TANGENT STIFFNESS MATRIX ROUTINE---------------------------------------------------*/
void Tangent_Stiffness_Matrix(std::vector<double> &theta_k,std::map<int,std::vector<double> > &surfaceRegion, std::vector<Element*> &elements, 
                                gmm::row_matrix< gmm::wsvector<double> > &KT,map<Node*, Node*> &NodesCorresp,std::vector<Node*> &nodes)
{

    gmm::row_matrix< gmm::wsvector<double> > Tmp(nodes.size(), nodes.size());//Temporaray matrix in the building of KT
    gmm::dense_matrix<double> alpha(2,2); // matrice alpha
    gmm::dense_matrix<double> beta(2,2); // matrice beta pour la conductivité
    gmm::dense_matrix<double> Ke_1(3,3); //First term in KT
    gmm::dense_matrix<double> Ke_2(3,3); //Second term in KT

    for(unsigned int i = 0; i < elements.size(); i++)//loop over the elements
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


            //alpha tensor
            alpha(0,0) = surfaceRegion[elements[i]->region][1];
            alpha(0,1) = surfaceRegion[elements[i]->region][2];
            alpha(1,0) = surfaceRegion[elements[i]->region][3];
            alpha(1,1) = surfaceRegion[elements[i]->region][4];

            //beta tensor
            beta(0,0) = surfaceRegion[elements[i]->region][5];
            beta(0,1) = surfaceRegion[elements[i]->region][6];
            beta(1,0) = surfaceRegion[elements[i]->region][7];
            beta(1,1) = surfaceRegion[elements[i]->region][8];

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
    }
    //cout << "In Tangent Stiffness Matrix, before copy" << endl;
    gmm::copy(Tmp,KT);
}//end of Tangent Stiffness Matrix routine


/*----------------------END CRITERION ROUTINE---------------------------*/
bool End_Criterion(std::vector<double> &RHS,double normRHS0)
{
    //Threshold value
    double eps = 1e-4;
    
    cout << "Relative residue = "<<gmm::vect_norm2(RHS)/normRHS0<<endl;
    if(gmm::vect_norm2(RHS)/normRHS0 >eps)
        return false; 

    return true;
}

