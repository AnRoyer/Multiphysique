#include <iostream>
#include <cstdio>
#include <vector>
#include <map>
#include <string>
#include <stdlib.h>
#include <mpi.h>
#include "gmm/gmm.h"
#include "fem.h"
#include "FE2.h"
#ifdef GMM_USES_MUMPS
#include "gmm/gmm_MUMPS_interface.h"
#endif

// On va faire du MPI !

using namespace std;

void FE2(std::vector<Node*> &nodes_micro, std::vector<Element*> &elements_micro, std::vector<Physical*> &physicals_micro,
         std::vector<Parameter*> &parameters_micro, std::map<Node*, std::vector<double> > &solutionTemperature_micro,
         std::map<Node*, std::vector<double> > &solutionFlux_micro, Periodique &conditions_micro, std::vector<Node*> &nodes_macro,
         std::vector<Element*> &elements_macro, std::vector<Physical*> &physicals_macro,std::vector<Parameter*> &parameters_macro,
         std::map<Node*, std::vector<double> > &solutionTemperature_macro, double eps, int argc, char ** argv)
{

// MPI initialization
MPI_Init(&argc, &argv);
MPI_Status status;
int nbproc, myrank ;
MPI_Comm_rank( MPI_COMM_WORLD, &myrank);
MPI_Comm_size( MPI_COMM_WORLD, &nbproc);

if (myrank == 0)
{
    //criterion for the overall FEM2 method.
    double criterionFEM2 = 10000; //criterion used for the convergence of the overall FEM2 method.
    double criterionFEM2_min = 0.1; //When convergence is reached.
    double criterionFEM0 =0;
}
// Declarations by Corentin :

if (myrank != 0)
{
	double T_mean;
	
	double u1, u2, u3; // Nodal temperatures of an element
	
	gmm::dense_matrix<double> J(2, 2); // Jacobian matrix
	gmm::dense_matrix<double> inverse_J(2, 2); // Its inverse
	double det_J; // Its determinant
	
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
	
	std::vector<double> q_Me(2);   // Flux moyen sur l'élément
	
	std::vector<double> gradT(2);// Gradient moyen
	gmm::dense_matrix<double> kappa_e(2,2);// Conductivité élémentaire.
	
	gmm::dense_matrix<double> element_stiffness(3, 3); // Matrix of the elementary stiffness K_e
	gmm::dense_matrix<double> total_stiffness(nodes_macro.size(), nodes_macro.size()); // Matrix of the stiffness for the whole domain K
	
	gmm::dense_matrix<double> a0(2, 2); // kappa_e_e times inverse of J in Eq. (9) of ddl 4 group B
	std::vector<double> a1(2); // a1 times gradPhi_1
	std::vector<double> a2(2); // a1 times gradPhi_2
	std::vector<double> a3(2); // a1 times gradPhi_3
	std::vector<double> b1(2); // inverse of J times gradPhi_i
	std::vector<double> b2(2); // inverse of J times gradPhi_2
	std::vector<double> b3(2); // inverse of J times gradPhi_3
	
	std::map<int, Parameter*> region_micro;//Stock le lien entre le numéro du physical de msh (stocker dans "physicals") et la valeur du parametre de "parametres" pour les régions de dimension 1 (ligne)

}


if (myrank == 0)
{

	/*Chargement des paramètres de la struc Parameter contenant les paramètres d'entrées.
	A chaque éléments de physicals, on lui associe l'élément de "parameters" correspondant. La correspondance est mappé dans linesRegion et surfaceRegion en fonction du type de paramètre
	*/
	
	for(unsigned int i = 0; i < physicals_micro.size(); i++)
	{
	    for(unsigned int j = 0; j < parameters_micro.size(); j++)
	    {
	        if(parameters_micro[j]->name == physicals_micro[i]->name)
	        {
	            if(parameters_micro[j]->dim != physicals_micro[i]->dim)//Verification si erreurs entre les deux fichiers
	            {
	                cout << "Error: file.phy and file.msh do not correspond" << endl;
	            }
	            else
	            {
	                region_micro[physicals_micro[i]->num] = parameters_micro[j];
	            }
	        }
	    }
	}
	
	double vol_micro=0;
	std::vector<double> c(nodes_micro.size());
	
	f_function(c, nodes_micro, elements_micro, region_micro, THERMALFLAG, 1);

	for(unsigned int j=0; j<nodes_micro.size(); j++)
	{
	    vol_micro += c[j];
	}
	
	cout << "vol_micro " << vol_micro << endl;
	
	//Computation of the f_i vector function of the heat generation in the macro domain. region_macro has to be initialized as in the fem function.
	map<int, Parameter*> region_macro;
	for(unsigned int i = 0; i < physicals_macro.size(); i++)
	{
	    for(unsigned int j = 0; j < parameters_macro.size(); j++)
	    {
	        if(parameters_macro[j]->name == physicals_macro[i]->name)
	        {
	            if(parameters_macro[j]->dim != physicals_macro[i]->dim)//Verification si erreurs entre les deux fichiers
	            {
	                cout << "Error: file.phy and file.msh do not correspond" << endl;
	            }
	            else
	            {
	                region_macro[physicals_macro[i]->num] = parameters_macro[j];
	            }
	        }
	    }
	}
	
	std::vector<double> f_i(nodes_macro.size());
	f_function(f_i, nodes_macro, elements_macro, region_macro, THERMALFLAG, 0); //utilisation of f_function from fem.cpp file.

}

//-----------------DEBUT DE LA BOUCLE A PARALLELISER-------------------------------

int nbElem = elements_macro.size();
std::vector<double> sol_u_tmp;
int source;	// Rank of the process sending a submatrix Ke
int numToSend; // Number of the element to send to the slave. The slave will be responsible for this macroscopic element.
int numElem; // Contains the number of the finite element associated to the Ke matrix a slave sends.
std::vector<Node*> nodesSlave(3); // Nodes : pointeur ou pas??? cfr. gmshio.h sinon 3 nodes à la main
std::vector<int> numNodesSlave(3); // Number of the nodes sent to the slave by the master
std::vector<int> numNodesMaster(3); // Number of the nodes sent by the slave to the master
int elementNumber;
std::vector<double> temperaturesSlave(3); // Temperatures of the nodes sent to the slave
std::vector<double> stiffnessMaster(9); // Will contain the 9 elements of the submatrice Ke (sent by process source)

int i_while = 0;

while(criterionFEM2 > criterionFEM2_min)
{
	if (myrank == 0) // Travail du maître
	{
		
		// INITIALIZATION : Now an element number is sent to each slave, so that all the slaves are busy at the beginning.
		int p; // Process number
		for (p = 1; p <= nbproc-1; p++)
		{
				numToSend = p - 1; // Process 0 will be busy with finite element 0,...
	    		
	    		// Nodes of the elements
	    		nodesSlave[0] = elements_macro[numToSend]->nodes[0];
	        	nodesSlave[1] = elements_macro[numToSend]->nodes[1];
        		nodesSlave[2] = elements_macro[numToSend]->nodes[2];
        		
        		// Numbers of these nodes
		        numNodesSlave[0] = nodesSlave[0]->num;
		        numNodesSlave[1] = nodesSlave[1]->num;
		        numNodesSlave[2] = nodesSlave[2]->num;
		        
		        // Temperatures at these nodes
			sol_u_tmp = solutionTemperature_macro[elements_macro[numToSend] -> nodes[0]];
			temperaturesSlave[0] = sol_u[0];
			sol_u_tmp = solutionTemperature_macro[elements_macro[numToSend] -> nodes[1]];
			temperaturesSlave[1] = sol_u[0];
			sol_u_tmp = solutionTemperature_macro[elements_macro[numToSend] -> nodes[2]];
			temperaturesSlave[2] = sol_u[0];
		        
			MPI_Send (numToSend, 1, MPI_INT, p, 38, MPI_COMM_WORLD);
			MPI_Send (temperaturesSlave, 3, MPI_DOUBLE, p, 42, MPI_COMM_WORLD);
		}
		// End of initialization
		
		/* Now the master loops over the remainder of the finite elements. When the slave process number source sends back a Ke matrix,
		master process asks it to work now with finite element i and put the submatrix in K.
		*/
		
		for (int i = nbproc - 1; i < nbElem; i++)
		{
			// Reception d'une sous matrice de la part d'un process a priori inconnu :
			MPI_Recv (stiffnessMaster, 9, MPI_DOUBLE, MPI_ANY_SOURCE, 39, MPI_COMM_WORLD, & status);
			source = status.MPI_SOURCE;
			// Réception des noeuds correspondants
			MPI_Recv (numNodesMaster, 3, MPI_INT, source, 40, MPI_COMM_WORLD, & status);
			
			// On demande maintenant au process source de se charger de l'élément i
			numToSend = i; // Same procedure than for the initializatino but with element number i
	    		
	    		// Nodes of the elements
	    		nodesSlave[0] = elements_macro[numToSend]->nodes[0];
	        	nodesSlave[1] = elements_macro[numToSend]->nodes[1];
        		nodesSlave[2] = elements_macro[numToSend]->nodes[2];
        		
        		// Numbers of these nodes
		        numNodesSlave[0] = nodesSlave[0]->num;
		        numNodesSlave[1] = nodesSlave[1]->num;
		        numNodesSlave[2] = nodesSlave[2]->num;
		        
		        // Temperatures at these nodes
			sol_u_tmp = solutionTemperature_macro[elements_macro[numToSend] -> nodes[0]];
			temperaturesSlave[0] = sol_u[0];
			sol_u_tmp = solutionTemperature_macro[elements_macro[numToSend] -> nodes[1]];
			temperaturesSlave[1] = sol_u[0];
			sol_u_tmp = solutionTemperature_macro[elements_macro[numToSend] -> nodes[2]];
			temperaturesSlave[2] = sol_u[0];
		        
			MPI_Send (numToSend, 1, MPI_INT, p, 38, MPI_COMM_WORLD);
			MPI_Send (temperaturesSlave, 3, MPI_DOUBLE, p, 42, MPI_COMM_WORLD);
			
			// Now compute the stiffness matrix with the one sent by the slave
		        total_stiffness(numNodesMaster[0], numNodesMaster[0]) = total_stiffness(numNodesMaster[0], numNodesMaster[0]) + stiffnessMaster [0];
		        total_stiffness(numNodesMaster[0], numNodesMaster[1]) = total_stiffness(numNodesMaster[0], numNodesMaster[1]) + stiffnessMaster [1];
		        total_stiffness(numNodesMaster[0], numNodesMaster[2]) = total_stiffness(numNodesMaster[0], numNodesMaster[2]) + stiffnessMaster [2];
		        total_stiffness(numNodesMaster[1], numNodesMaster[0]) = total_stiffness(numNodesMaster[1], numNodesMaster[0]) + stiffnessMaster [3];
		        total_stiffness(numNodesMaster[1], numNodesMaster[1]) = total_stiffness(numNodesMaster[1], numNodesMaster[1]) + stiffnessMaster [4];
		        total_stiffness(numNodesMaster[1], numNodesMaster[2]) = total_stiffness(numNodesMaster[1], numNodesMaster[2]) + stiffnessMaster [5];
		        total_stiffness(numNodesMaster[2], numNodesMaster[0]) = total_stiffness(numNodesMaster[2], numNodesMaster[0]) + stiffnessMaster [6];
		        total_stiffness(numNodesMaster[2], numNodesMaster[1]) = total_stiffness(numNodesMaster[2], numNodesMaster[1]) + stiffnessMaster [7];
		        total_stiffness(numNodesMaster[2], numNodesMaster[2]) = total_stiffness(numNodesMaster[2], numNodesMaster[2]) + stiffnessMaster [8];
	
		}
		
		/* Now the last element matrices have been "ordered", on need to wait that each of the slaves gives backs its Ke matrix,
		   then add it to the total stifness matrix and send the finite element number -1 to each slave to tell them their work is finished.
		*/
		int nodeFlag = -1;
		std::vector<double> temperatureFlag(3) = {-1.0, -1.0, -1.0};
		
		for (p = 1; p <= nbproc - 1; p++)
		{
			MPI_Recv (stiffnessMaster, 9, MPI_DOUBLE, p, 39, MPI_COMM_WORLD, & status);
			MPI_Recv (numNodesMaster, 3, MPI_INT, p, 40, MPI_COMM_WORLD, & status);
			// Dire aux process que c'est terminé
			MPI_Send (nodeFlag, 1, MPI_INT, p, 38, MPI_COMM_WORLD);
			MPI_Send (temperatureFlag, 3, MPI_DOUBLE, p, 42, MPI_COMM_WORLD);
			
		        total_stiffness(numNodesMaster[0], numNodesMaster[0]) = total_stiffness(numNodesMaster[0], numNodesMaster[0]) + stiffnessMaster [0];
		        total_stiffness(numNodesMaster[0], numNodesMaster[1]) = total_stiffness(numNodesMaster[0], numNodesMaster[1]) + stiffnessMaster [1];
		        total_stiffness(numNodesMaster[0], numNodesMaster[2]) = total_stiffness(numNodesMaster[0], numNodesMaster[2]) + stiffnessMaster [2];
		        total_stiffness(numNodesMaster[1], numNodesMaster[0]) = total_stiffness(numNodesMaster[1], numNodesMaster[0]) + stiffnessMaster [3];
		        total_stiffness(numNodesMaster[1], numNodesMaster[1]) = total_stiffness(numNodesMaster[1], numNodesMaster[1]) + stiffnessMaster [4];
		        total_stiffness(numNodesMaster[1], numNodesMaster[2]) = total_stiffness(numNodesMaster[1], numNodesMaster[2]) + stiffnessMaster [5];
		        total_stiffness(numNodesMaster[2], numNodesMaster[0]) = total_stiffness(numNodesMaster[2], numNodesMaster[0]) + stiffnessMaster [6];
		        total_stiffness(numNodesMaster[2], numNodesMaster[1]) = total_stiffness(numNodesMaster[2], numNodesMaster[1]) + stiffnessMaster [7];
		        total_stiffness(numNodesMaster[2], numNodesMaster[2]) = total_stiffness(numNodesMaster[2], numNodesMaster[2]) + stiffnessMaster [8];
			}// End of the loop over the processes.
	
	}// End of the master.
	
	std::vector<int> numNodesSlave(3); // Number of the nodes sent to the slave
	std::vector<Node*> nodesSlave(3);
	
    for(unsigned int i = 0; i < elements_macro.size(); i++)
    {
        if(elements_macro[i]->type == 2)//If triangle
        {
        	// vérifier qu'on a pas un -1 qqpart
        if(myrank != 0){
        	MPI_Recv (numElem, 1, MPI_INT, 0, 38, MPI_COMM_WORLD, & status);
			MPI_Recv (temperaturesSlave, 3, MPI_DOUBLE, 0, 42, MPI_COMM_WORLD, & status);

			if(numElem == -1)
			{
				break;
			}
			
			// On récupère les noeuds relatifs à l'élément envoyé
			Node *n1 = elements_macro[numElem]->nodes[0];
            Node *n2 = elements_macro[numElem]->nodes[1];
            Node *n3 = elements_macro[numElem]->nodes[2];
			
			// Leurs coordonnées
            double x1 = n1->x;
            double y1 = n1->y;
            double x2 = n2->x;
            double y2 = n2->y;
            double x3 = n3->x;
            double y3 = n3->y;

			// Leurs numéros
            int num_n1, num_n2, num_n3;
            num_n1 = n1->num;
            num_n2 = n2->num;
            num_n3 = n3->num;

            // Jacobian
            J (0, 0) = x2 - x1;
            J (0, 1) = y2 - y1;
            J (1, 0) = x3 - x1;
            J (1, 1) = y3 - y1;

            // Determinant of the Jacobian 
            det_J = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
            
            //Inverse of the Jacobian Matrix
            inverse_J(0,0) = 1.0/det_J*(y3-y1);
            inverse_J(0,1) = 1.0/det_J*(y1-y2);
            inverse_J(1,0) = 1.0/det_J*(x1-x3);
            inverse_J(1,1) = 1.0/det_J*(x2-x1);

            gmm::mult(inverse_J, gradPhi1_red, gradPhi1);
            gmm::mult(inverse_J, gradPhi2_red, gradPhi2);
            gmm::mult(inverse_J, gradPhi3_red, gradPhi3);

            u1 = temperaturesSlave[0];
            u2 = temperaturesSlave[1];
            u3 = temperaturesSlave[2];
            
            T_mean = (1.0/3.0)*(u1 + u2 + u3);
            
            conditions_micro.meanTemperature = T_mean;
            conditions_micro.xGradient = u1*gradPhi1[0] +  u2*gradPhi2[0] + u3*gradPhi3[0];
            conditions_micro.yGradient = u1*gradPhi1[1] +  u2*gradPhi2[1] + u3*gradPhi3[1];

            gradT[0] = u1*gradPhi1[0] +  u2*gradPhi2[0] + u3*gradPhi3[0];
            gradT[1] = u1*gradPhi1[1] +  u2*gradPhi2[1] + u3*gradPhi3[1];

            fem(nodes_micro, elements_micro, physicals_micro, parameters_micro, solutionTemperature_micro, solutionFlux_micro,
                THERMALFLAG, PERIODICFLAG, conditions_micro, eps);

            q_Me[0] = 0.0;
            q_Me[1] = 0.0;
            Average_flux(solutionTemperature_micro, region_micro, elements_micro, q_Me, vol_micro);

            conductivityTensor(q_Me, gradT, kappa_e);
        
        	// Computing elementary stiffness matrix
            gmm::mult(kappa_e, inverse_J, a0);

            gmm::mult(a0, gradPhi1_red, a1);
            gmm::mult(a0, gradPhi2_red, a2);
            gmm::mult(a0, gradPhi3_red, a3);

            gmm::mult(inverse_J, gradPhi1_red, b1);
            gmm::mult(inverse_J, gradPhi2_red, b2);
            gmm::mult(inverse_J, gradPhi3_red, b3);

            double element_stiffness_temp;
            element_stiffness_temp = gmm::vect_sp(a1, b1);
            element_stiffness(0, 0) = element_stiffness_temp;
            element_stiffness_temp = gmm::vect_sp(a2, b1);
            element_stiffness(0, 1) = element_stiffness_temp;
            element_stiffness_temp = gmm::vect_sp(a3, b1);
            element_stiffness(0, 2) = element_stiffness_temp;
            element_stiffness_temp = gmm::vect_sp(a1, b2);
            element_stiffness(1, 0) = element_stiffness_temp;
            element_stiffness_temp = gmm::vect_sp(a2, b2);
            element_stiffness(1, 1) = element_stiffness_temp;
            element_stiffness_temp = gmm::vect_sp(a3, b2);
            element_stiffness(1, 2) = element_stiffness_temp;
            element_stiffness_temp = gmm::vect_sp(a1, b3);
            element_stiffness(2, 0) = element_stiffness_temp;
            element_stiffness_temp = gmm::vect_sp(a2, b3);
            element_stiffness(2, 1) = element_stiffness_temp;
            element_stiffness_temp = gmm::vect_sp(a3, b3);
            element_stiffness(2, 2) = element_stiffness_temp;
            
            gmm::scale(element_stiffness, 0.5*det_J);
            
        	stiffnessMaster[0] = element_stiffness(0, 0);
        	stiffnessMaster[1] = element_stiffness(0, 1);
        	stiffnessMaster[2] = element_stiffness(0, 2);
        	stiffnessMaster[3] = element_stiffness(1, 0);
        	stiffnessMaster[4] = element_stiffness(1, 1);
        	stiffnessMaster[5] = element_stiffness(1, 2);
        	stiffnessMaster[6] = element_stiffness(2, 0);
        	stiffnessMaster[7] = element_stiffness(2, 1);
            stiffnessMaster[8] = element_stiffness(2, 2);

            MPI_Send (stiffnessMaster, 9, MPI_DOUBLE, 0, 39, MPI_COMM_WORLD, & status);
			MPI_Send (numNodesMaster, 3, MPI_INT, 0, 40, MPI_COMM_WORLD, & status);
        }// end if myrank

            // Assembling the K_ij
            // total_stiffness(num_n1-1, num_n1-1) = total_stiffness(num_n1-1, num_n1-1) + element_stiffness(0, 0);
            // total_stiffness(num_n1-1, num_n2-1) = total_stiffness(num_n1-1, num_n2-1) + element_stiffness(0, 1);
            // total_stiffness(num_n1-1, num_n3-1) = total_stiffness(num_n1-1, num_n3-1) + element_stiffness(0, 2);
            // total_stiffness(num_n2-1, num_n1-1) = total_stiffness(num_n2-1, num_n1-1) + element_stiffness(1, 0);
            // total_stiffness(num_n2-1, num_n2-1) = total_stiffness(num_n2-1, num_n2-1) + element_stiffness(1, 1);
            // total_stiffness(num_n2-1, num_n3-1) = total_stiffness(num_n2-1, num_n3-1) + element_stiffness(1, 2);
            // total_stiffness(num_n3-1, num_n1-1) = total_stiffness(num_n3-1, num_n1-1) + element_stiffness(2, 0);
            // total_stiffness(num_n3-1, num_n2-1) = total_stiffness(num_n3-1, num_n2-1) + element_stiffness(2, 1);
            // total_stiffness(num_n3-1, num_n3-1) = total_stiffness(num_n3-1, num_n3-1) + element_stiffness(2, 2);
        }//end if triangle
    }//end for i

    for(unsigned int i = 0; i < elements_macro.size(); i++)
    {
        if(elements_macro[i]->type == 1)//If line
        {
            if(region_macro.count(elements_macro[i]->region) == 1 && region_macro[elements_macro[i]->region]->temperature != -1)
            {
                for(unsigned int j = 0; j < elements_macro[i]->nodes.size(); j++)
                {
                    for(unsigned int k = 0; k < nodes_macro.size(); k++)
                    {
                        if(k == elements_macro[i]->nodes[j]->num-1)
                        {
                            total_stiffness(elements_macro[i]->nodes[j]->num-1, k) = 1;
                        }
                        else
                        {
                            total_stiffness(elements_macro[i]->nodes[j]->num-1, k) = 0;
                        }
                    }
                }
            }
        }
    }

    for(unsigned int i = 0; i < elements_macro.size(); i++)
    {
        if(elements_macro[i]->type == 1 && region_macro[elements_macro[i]->region]->temperature != -1)//If line
        {
            if(region_macro.count(elements_macro[i]->region) == 1)//If linesRegion contains elements[i]->region
            {
                for(unsigned int j = 0; j < elements_macro[i]->nodes.size(); j++)
                {
                    f_i[elements_macro[i]->nodes[j]->num-1] = region_macro[elements_macro[i]->region]->temperature;
                }
            }
        }
    }

    std::vector<double> u_guess_vec(1);
    std::vector<double> u_guess(nodes_macro.size());

    for(unsigned int l=0; l<nodes_macro.size(); l++)
    {
        u_guess_vec = solutionTemperature_macro[nodes_macro[l]];
        u_guess[l] = u_guess_vec[0];
        //cout << u_guess[l] << " ";
    }
    //cout << endl;

    std::vector<double> error(nodes_macro.size());
    std::vector<double> temporary(nodes_macro.size());

    gmm::mult(total_stiffness, u_guess, temporary);

    /*for(unsigned int m=0; m<nodes_macro.size(); m++)
    {
        for(unsigned int n=0; n<nodes_macro.size(); n++)
        {
            double out = total_stiffness(m,n);
            cout << out << "		";
        }
        cout << endl;
    }
    cout << endl;

    for(unsigned int l=0; l<nodes_macro.size(); l++)
    {
        cout << f_i[l] << " ";
    }
    cout << endl;*/

    gmm::add(temporary, gmm::scaled(f_i, -1.0), error);
    criterionFEM2 = gmm::vect_norm2(error);
    if(i_while ==0)
        criterionFEM0 = gmm::vect_norm2(error);

    criterionFEM2 = criterionFEM2/criterionFEM0;
    cout << "criterionFEM2 " << criterionFEM2 << endl; //acorriger
    if(criterionFEM2 > criterionFEM2_min)
    {
        #ifdef GMM_USES_MUMPS
        std::cout << "solving linear system with MUMPS\n";
            gmm::csr_matrix<double> total_stiffness_csr; //à changer
            gmm::copy(total_stiffness,total_stiffness_csr);
            gmm::MUMPS_solve(total_stiffness_csr, u_guess, f_i);
        #else
        std::cout << "solving linear system with gmm::lu_solve\n";
            gmm::lu_solve(total_stiffness, u_guess, f_i);
        #endif
        //gmm::lu_solve(total_stiffness, u_guess, f_i);
        for(unsigned int l=0; l<nodes_macro.size(); l++)
        {
            u_guess_vec[0] = u_guess[l];
            solutionTemperature_macro[nodes_macro[l]] = u_guess_vec;
            //cout << u_guess[l] << " ";
        }
        //cout << endl;
    }
    i_while ++;

}//end while.

}//end function.


void conductivityTensor(std::vector<double> &q, std::vector<double> &gradT, gmm::dense_matrix<double> &kappa)
{

	// Prendre l'opposé du gradient

	std::vector<double> gradT2;

	gradT2 = gradT;

	gmm::scale (gradT2, -1.0);

	// Computation of the norms :

	double n_gradT2 = gmm::vect_norm2(gradT2);

        if(n_gradT2 < 1e-10)
        {
    	    kappa (0,0) = 0.;
	    kappa (0, 1) = 0.;
	    kappa (1, 0) = 0.;
	    kappa (1, 1) = 0.;
            return;
        }

	double n_q = gmm::vect_norm2(q);

	double C = n_q/n_gradT2;

	// Orientation of the opposite of the gradient and the heat flux :

	double prod_scal = gmm::vect_sp(q, gradT2);

	//double alpha_gradT2 = atan (gradT2[1]/gradT2[0]);
	//double alpha_q = atan (q[1]/q[0]);



	// angle between these vector :

	double cosinus = prod_scal/(n_q*n_gradT2);

	if (cosinus >1)
		cosinus = 1;
	if (cosinus <-1)
		cosinus = -1;

        double prod_vect = gradT2[0]*q[1]-gradT2[1]*q[0];

	double sinus = prod_vect/(n_q*n_gradT2);

	if (sinus >1)
		sinus = 1;
	if (sinus <-1)
		sinus = -1;

	kappa (0,0) = C * cosinus;
	kappa (0, 1) = -C * sinus;
	kappa (1, 0) = C * sinus;
	kappa (1, 1) = C * cosinus;

            /*if(fabs(kappa (0,0)-50)>1e-3 || fabs(kappa (0,1))>1e-3 || fabs(kappa (1,0))>1e-3 || fabs(kappa (1,1)-50)>1e-3)
            {
                cout << "error" << endl;
                cout << "graT " << gradT << endl;
                cout << "graT2 " << gradT2 << endl;
                cout << "q " << q << endl;
                cout << "sinus " << sinus << endl;
                cout << "cosinus " << cosinus << endl;
            for(unsigned int m=0; m<2; m++)
            {
                for(unsigned int n=0; n<2; n++)
                {
                    double out = kappa(m,n);
                    cout << out << "		";
                }
                cout << endl;
            }
                exit(1);
            }*/
}


