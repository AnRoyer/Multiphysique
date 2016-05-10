#include <iostream>
#include <cstdio>
#include <vector>
#include <map>
#include <string>
#include <mpi.h>
#include "gmm/gmm.h"
#include "fem.h"
#include "FE2.h"
#include "physicalio.h"
#include "NewtonRaphson.h"
#ifdef GMM_USES_MUMPS
#include "gmm/gmm_MUMPS_interface.h"
#endif

// On va faire du MPI !

using namespace std;

void FE2(std::vector<Node*> &nodes_micro, std::vector<Element*> &elements_micro, std::vector<Physical*> &physicals_micro,
         std::vector<Parameter*> &parameters_micro, std::map<Node*, std::vector<double> > &solutionTemperature_micro,
         std::map<Node*, std::vector<double> > &solutionFlux_micro, Periodique &conditions_micro, std::vector<Node*> &nodes_macro,
         std::vector<Element*> &elements_macro, std::vector<Physical*> &physicals_macro,std::vector<Parameter*> &parameters_macro,
         std::map<Node*, std::vector<double> > &solutionTemperature_macro, std::map<Node*, std::vector<double> > &solutionFlux_macro, 
		 double eps, int argc, char ** argv, int &methodFE2, FemFlag thermalOrElectrical, Type type)
{

// MPI initialization
MPI_Init(&argc, &argv);
MPI_Status status;
int nbproc, myrank ;
MPI_Comm_rank( MPI_COMM_WORLD, &myrank);
MPI_Comm_size( MPI_COMM_WORLD, &nbproc);

//Displaying information
if(myrank == 0)
{
	cout << endl;
    cout << "\t############################################################" << endl;
    cout << "\t############################################################" << endl;
    cout << "\t##                                                        ##" << endl;
    cout << "\t##                                             222222     ##" << endl;
    cout << "\t##                                                 22     ##" << endl;
    cout << "\t##     FFFFFFFFFFFFFFFF    EEEEEEEEEEEEEEEE    222222     ##" << endl;
    cout << "\t##     FF                  EE                  22         ##" << endl;
    cout << "\t##     FF                  EE                  222222     ##" << endl;
    cout << "\t##     FF                  EE                             ##" << endl;
    cout << "\t##     FF                  EE                             ##" << endl;
    cout << "\t##     FFFFFFFFFFFF        EEEEEEEEEEEE                   ##" << endl;
    cout << "\t##     FF                  EE                             ##" << endl;
    cout << "\t##     FF                  EE                             ##" << endl;
    cout << "\t##     FF                  EE                             ##" << endl;
    cout << "\t##     FF                  EE                             ##" << endl;
    cout << "\t##     FF                  EEEEEEEEEEEEEEEE               ##" << endl;
    cout << "\t##                                                        ##" << endl;
    cout << "\t############################################################" << endl;
    cout << "\t############################################################" << endl << endl << endl;

    //Affichage des infos principales

	if(type == FE2withDIRICHLET) cout << "Calling the FE2 method with DIRICHLET conditions on the macroscopic domain." << endl;
	if(type == FE2withVONNEUMANN) cout << "Calling the FE2 method with VONNEUMANN conditions on the macroscopic domain." << endl;
	if(type == FE2withPERIODIC) cout << "Calling the FE2 method with PERIODIC conditions on the macroscopic domain." << endl;
    cout << endl;
    //cout << "Reads " << nodes_macro.size() << " nodes and " << elements_macro.size() << " elements for the macroscopic domain." << endl;
    //cout << "Reads " << nodes_micro.size() << " nodes and " << elements_micro.size() << " elements for the microscopic domain." << endl << endl;

    /*cout << "For the macroscopic domain:" << endl;

    for(unsigned int i = 0; i < parameters_macro.size(); i++)
    {
        if(parameters_macro[i]->dim == 1 )
        {
            if(parameters_macro[i]->temperature != -1 && parameters_macro[i]->voltage != -1)
            {
                cout << "Parameter " << parameters_macro[i]->name << " has temperature of " << parameters_macro[i]->temperature << " K and an electrical potential of " << parameters_macro[i]->voltage << " V." << endl;
            }
            else if(parameters_macro[i]->temperature != -1 && parameters_macro[i]->voltage == -1)
            {
                cout << "Parameter " << parameters_macro[i]->name << " has temperature of " << parameters_macro[i]->temperature << " K." << endl;
            }
            else if(parameters_macro[i]->temperature == -1 && parameters_macro[i]->voltage != -1)
            {
                cout << "Parameter " << parameters_macro[i]->name << " has an electrical potential of " << parameters_macro[i]->voltage << " V." << endl;
            }

            if(parameters_macro[i]->fluxTemperature != -1 && parameters_macro[i]->fluxVoltage != -1)
            {
                cout << "Parameter " << parameters_macro[i]->name << " has heat flux of " << parameters_macro[i]->fluxTemperature << " W/m and an electrical flux of " << parameters_macro[i]->fluxVoltage << "." << endl;
            }
            else if(parameters_macro[i]->fluxTemperature != -1 && parameters_macro[i]->fluxVoltage == -1)
            {
                cout << "Parameter " << parameters_macro[i]->name << " has heat flux of " << parameters_macro[i]->fluxTemperature << " W/m." << endl;
            }
            else if(parameters_macro[i]->fluxTemperature == -1 && parameters_macro[i]->fluxVoltage != -1)
            {
                cout << "Parameter " << parameters_macro[i]->name << " has an electrical potential of " << parameters_macro[i]->fluxVoltage << "." << endl;
            }
        }
        else if(parameters_macro[i]->dim == 2)
        {
            cout << "Parameter " << parameters_macro[i]->name << " has:" << endl;

            cout << "\t - A thermal conductivity [W/mK] of:" << endl;
            for(unsigned int j = 0; j < parameters_macro[i]->thermalConductivity.size(); j++)
            {
                cout << "\t\t" << parameters_macro[i]->thermalConductivity[j]->name << " = " << endl;
                cout << "\t\t|" << parameters_macro[i]->thermalConductivity[j]->conductivity[0][0] << "\t" << parameters_macro[i]->thermalConductivity[j]->conductivity[0][1] << "|" << endl;
                cout << "\t\t|" << parameters_macro[i]->thermalConductivity[j]->conductivity[1][0] << "\t" << parameters_macro[i]->thermalConductivity[j]->conductivity[1][1] << "|" << endl << endl;
            }

            cout << "\t - A heat production of " << parameters_macro[i]->thermalGeneration << " W/m^2" << endl;

            cout << "\t - An electrical conductivity of :" << endl;
            for(unsigned int j = 0; j < parameters_macro[i]->electricalConductivity.size(); j++)
            {
                cout << "\t\t" << parameters_macro[i]->electricalConductivity[j]->name << " = " << endl;
                cout << "\t\t|" << parameters_macro[i]->electricalConductivity[j]->conductivity[0][0] << "\t" << parameters_macro[i]->electricalConductivity[j]->conductivity[0][1] << "|" << endl;
                cout << "\t\t|" << parameters_macro[i]->electricalConductivity[j]->conductivity[1][0] << "\t" << parameters_macro[i]->electricalConductivity[j]->conductivity[1][1] << "|" << endl << endl;
            }
            cout << "\t - A change in charge density of " << parameters_macro[i]->electricalGeneration << " A/m^3." << endl;
        }
    }
	cout << endl << endl;
    cout << "For the microscopic domain:" << endl;

	for(unsigned int i = 0; i < parameters_micro.size(); i++)
	{
	    if(parameters_micro[i]->dim == 2)
	    {
	        cout << "Parameter " << parameters_micro[i]->name << " has:" << endl;

	        cout << "\t - A thermal conductivity [W/mK] of:" << endl;
	        for(unsigned int j = 0; j < parameters_micro[i]->thermalConductivity.size(); j++)
	        {
	            cout << "\t\t" << parameters_micro[i]->thermalConductivity[j]->name << " = " << endl;
	            cout << "\t\t|" << parameters_micro[i]->thermalConductivity[j]->conductivity[0][0] << "\t" << parameters_micro[i]->thermalConductivity[j]->conductivity[0][1] << "|" << endl;
	            cout << "\t\t|" << parameters_micro[i]->thermalConductivity[j]->conductivity[1][0] << "\t" << parameters_micro[i]->thermalConductivity[j]->conductivity[1][1] << "|" << endl << endl;
	        }

	        cout << "\t - An electrical conductivity of:" << endl;
	        for(unsigned int j = 0; j < parameters_micro[i]->electricalConductivity.size(); j++)
	        {
	            cout << "\t\t" << parameters_micro[i]->electricalConductivity[j]->name << " = " << endl;
	            cout << "\t\t|" << parameters_micro[i]->electricalConductivity[j]->conductivity[0][0] << "\t" << parameters_micro[i]->electricalConductivity[j]->conductivity[0][1] << "|" << endl;
	            cout << "\t\t|" << parameters_micro[i]->electricalConductivity[j]->conductivity[1][0] << "\t" << parameters_micro[i]->electricalConductivity[j]->conductivity[1][1] << "|" << endl << endl;
	        }
	    }
	}
	cout << endl << endl;*/
}// end of displaying

//criterion for the overall FEM2 method.
double criterionFE2 = 10000; 
double criterionFE2_min = 1e-5;
double criterionFE2_0 =0;
std::vector<double> error(nodes_macro.size());

double T_mean;
int method = methodFE2;

if(myrank == 0)
{
	if(method ==1) cout << "Using method 1 for assembling the stiffness matrix." << endl;
	else if(method ==2) cout << "Using method 2 for assembling the stiffness matrix." << endl;
	cout << endl;
}
	
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

std::vector <double> q_int_e(3); // Elementary q_int vector
std::vector <double> q_int(nodes_macro.size());

std::map<int, Parameter*> region_micro;//Stock le lien entre le numéro du physical de msh (stocker dans "physicals") et la valeur du parametre de "parametres" pour les régions de dimension 1 (ligne)


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

//cout << "vol_micro " << vol_micro << endl;

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
if(method == 2) gmm::scale(f_i,-1);

//-----------------DEBUT DE LA BOUCLE A PARALLELISER-------------------------------

int nbElem = elements_macro.size();
std::vector<double> sol_u_tmp;
int source;	// Rank of the process sending a submatrix Ke
int numToSend; // Number of the element to send to the slave. The slave will be responsible for this macroscopic element.
int numElem; // Contains the number of the finite element associated to the Ke matrix a slave sends.
std::vector<int> numNodesMaster(3); // Number of the nodes sent by the slave to the master
int elementNumber;
std::vector<double> temperaturesSlave(3); // Temperatures of the nodes sent to the slave
std::vector<double> stiffnessMaster(9); // Will contain the 9 elements of the submatrice Ke (sent by process source)


if (myrank == 0) // Travail du maître
{
    int i_while = 0;
	while(criterionFE2 > criterionFE2_min)
	{
		cout << "Running Newton Raphson on the RVEs..." << endl;
		gmm::clear(total_stiffness);
		std::vector<double> f_i(nodes_macro.size());
		f_function(f_i, nodes_macro, elements_macro, region_macro, THERMALFLAG, 0); //utilisation of f_function from fem.cpp file.
		if(method == 2) gmm::scale(f_i,-1);

		std::vector <double> q_int(nodes_macro.size());
		// INITIALIZATION : Now an element number is sent to each slave, so that all the slaves are busy at the beginning.
		int p; // Process number
        int k;
		for (k = 0; k< elements_macro.size(); k++)
        {
            if(elements_macro[k]->type == 2)
                break;
        } 
		for (p = 1; p <= nbproc-1; p++)
		{
	        numToSend = p-1+k; // Process 1 will be busy with finite element 0,...
    		// Nodes of the elements
		    // Temperatures at these nodes
			sol_u_tmp = solutionTemperature_macro[elements_macro[numToSend] -> nodes[0]];
			temperaturesSlave[0] = sol_u_tmp[0];
			sol_u_tmp = solutionTemperature_macro[elements_macro[numToSend] -> nodes[1]];
			temperaturesSlave[1] = sol_u_tmp[0];
			sol_u_tmp = solutionTemperature_macro[elements_macro[numToSend] -> nodes[2]];
			temperaturesSlave[2] = sol_u_tmp[0];

			MPI_Send (& numToSend, 1, MPI_INT, p, 38, MPI_COMM_WORLD);
			MPI_Send (&temperaturesSlave[0], 3, MPI_DOUBLE, p, 42, MPI_COMM_WORLD);
		}
		// End of initialization
		
		/* Now the master loops over the remainder of the finite elements. When the slave process number source sends back a Ke matrix,
		master process asks it to work now with finite element i and put the submatrix in K.
		*/
		
		for (int i = nbproc - 1+k; i < nbElem + nbproc-1; i++)
		{
			//cout << "i " << i << endl;

			// Reception d'une sous matrice de la part d'un process a priori inconnu :
			MPI_Recv (&stiffnessMaster[0], 9, MPI_DOUBLE, MPI_ANY_SOURCE, 39, MPI_COMM_WORLD, & status);
			source = status.MPI_SOURCE;

			// Réception des noeuds correspondants
			MPI_Recv (&numNodesMaster[0], 3, MPI_INT, source, 40, MPI_COMM_WORLD, & status);

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
			
			MPI_Recv (&q_int_e[0], 3, MPI_DOUBLE, source, 41, MPI_COMM_WORLD, & status);

			// Assemling q_int : 
			q_int[numNodesMaster[0]] = q_int[numNodesMaster[0]] + q_int_e[0];
			q_int[numNodesMaster[1]] = q_int[numNodesMaster[1]] + q_int_e[1];
			q_int[numNodesMaster[2]] = q_int[numNodesMaster[2]] + q_int_e[2];

			if(i<nbElem)
			{
				// On demande maintenant au process source de se charger de l'élément i
				numToSend = i; // Same procedure than for the initializatino but with element number i
				    
				// Temperatures at these nodes
				sol_u_tmp = solutionTemperature_macro[elements_macro[numToSend] -> nodes[0]];
				temperaturesSlave[0] = sol_u_tmp[0];
				sol_u_tmp = solutionTemperature_macro[elements_macro[numToSend] -> nodes[1]];
				temperaturesSlave[1] = sol_u_tmp[0];
				sol_u_tmp = solutionTemperature_macro[elements_macro[numToSend] -> nodes[2]];
				temperaturesSlave[2] = sol_u_tmp[0];
				    
				MPI_Send (&numToSend, 1, MPI_INT, source, 38, MPI_COMM_WORLD);
				MPI_Send (&temperaturesSlave[0], 3, MPI_DOUBLE, source, 42, MPI_COMM_WORLD);
			}
	
		}

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

		//cout << "qint" << q_int << endl;
		//cout << "f_i" << f_i << endl;
		error.clear();
		for(unsigned int i=0;i<nodes_macro.size();i++)
		{
			error.push_back(q_int[i] - f_i[i]);
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
						error[elements_macro[i]->nodes[j]->num-1] = 0.;
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
		}

		int flag =1;
		if(flag)
		{
			std::vector<double> flag_temp(nodes_macro.size());
			gmm::mult(total_stiffness,u_guess,error);
            gmm::add(flag_temp,f_i,flag_temp);
			gmm::scale(flag_temp,-1);
			gmm::add(error,flag_temp,error);
		}

		if(i_while ==0)
		{
		    criterionFE2_0 = gmm::vect_norm2(error);
			if(criterionFE2_0 < 1e-6) criterionFE2_0 =1; //a etudier
		}
		criterionFE2 = gmm::vect_norm2(error);
		criterionFE2 = criterionFE2/criterionFE2_0;
		cout << "FE2 relative residue = " << criterionFE2 << endl; //acorriger

		if(criterionFE2 > criterionFE2_min)
		{
		    #ifdef GMM_USES_MUMPS
		        std::cout << "solving linear system with MUMPS\n";
		        gmm::csr_matrix<double> total_stiffness_csr; //à changer
		        gmm::copy(total_stiffness,total_stiffness_csr);
		        gmm::MUMPS_solve(total_stiffness_csr, u_guess, f_i);
				//cout << u_guess << endl;
		    #else
		        std::cout << "solving linear system with gmm::lu_solve\n";
		        gmm::lu_solve(total_stiffness, u_guess, f_i);
		    #endif

		    for(unsigned int l=0; l<nodes_macro.size(); l++)
		    {
		        u_guess_vec[0] = u_guess[l];
		        solutionTemperature_macro[nodes_macro[l]] = u_guess_vec;
		    }
		}
		//computation of the flux in the macroscopic domain.
		std::vector<double> q_m_x(nodes_macro.size());
		std::vector<double> q_m_y(nodes_macro.size());
		std::vector<double> qint(nodes_macro.size());

		Internal_flux(u_guess, region_macro, elements_macro, qint, q_m_x, q_m_y);//Chargement de la solution pour la solution du flux

		for(unsigned int i = 0; i < nodes_macro.size(); i++)
		{
		    std::vector<double> val2(3);
		    val2[0] = q_m_x[i];
		    val2[1] = q_m_y[i];
		    val2[2] = 0;
		    solutionFlux_macro[nodes_macro[i]] = val2;
		}

		//writing the solution.
		writeMSH((char*)"solutionTemperature.pos", solutionTemperature_macro);
		writeMSH((char*)"solutionFlux.pos", solutionFlux_macro);
		
        i_while ++;

		if(criterionFE2 > 0.1 && i_while == 5)
		{
			cout << endl;
			cout << "Anisotropy on the microscopic domain is too important for method 1." << endl;
			i_while = 0;
			method = 2;
			cout << "Press enter to use method 2.";
			cout << endl;
  			std::cin.ignore();
			int elementFlag = -2;
    		for(unsigned int i=0; i<nbproc;i++)
    		{
				MPI_Send (&elementFlag, 1, MPI_INT, i, 38, MPI_COMM_WORLD);
   			}
		}
    }//end while

	//flag used to stop the subprocesse from waiting when the FE2 method is over.
	int elementFlag = -1;
    for(unsigned int i=0; i<nbproc;i++)
    {
		MPI_Send (&elementFlag, 1, MPI_INT, i, 38, MPI_COMM_WORLD);
    }
	cout << endl;
	if(type == FE2withDIRICHLET) cout << "The problem has been solved in two scales with DIRICHLET conditions on the macroscopic domain." << endl;
	if(type == FE2withVONNEUMANN) cout <<"The problem has been solved in two scales with VONNEUMANN conditions on the macroscopic domain." << endl;
	if(type == FE2withPERIODIC) cout <<  "The problem has been solved in two scales with PERIODIC conditions on the macroscopic domain."  << endl;
	cout << endl;	

}//end master.

if(myrank !=0)
{
    while(1)
    {
    	MPI_Recv (&numElem, 1, MPI_INT, 0, 38, MPI_COMM_WORLD, & status);
        if(numElem == -1) break;
		else if(numElem == -2)
		{
			method = 2;
			MPI_Recv (&numElem, 1, MPI_INT, 0, 38, MPI_COMM_WORLD, & status);
		}
		MPI_Recv (&temperaturesSlave[0], 3, MPI_DOUBLE, 0, 42, MPI_COMM_WORLD, & status);
	
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
        numNodesMaster[0] = n1->num-1;
        numNodesMaster[1] = n2->num-1;
        numNodesMaster[2] = n3->num-1;

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
		    
		fem(nodes_micro, elements_micro, physicals_micro, parameters_micro, solutionTemperature_micro, solutionFlux_micro, THERMALFLAG, PERIODICFLAG, conditions_micro, eps, type);

        q_Me[0] = 0.0;
        q_Me[1] = 0.0;
        Average_flux(solutionTemperature_micro, region_micro, elements_micro, q_Me, vol_micro, thermalOrElectrical);

		// Computation of the q_int_e vector (q_int of the element) :
		q_int_e[0] = gmm::vect_sp(q_Me, gradPhi1);
		q_int_e[1] = gmm::vect_sp(q_Me, gradPhi2);
		q_int_e[2] = gmm::vect_sp(q_Me, gradPhi3);
        gmm::scale(q_int_e,det_J/2);

		if(method == 1)
		{
		    conductivityTensor(q_Me, gradT, kappa_e);
		
			// Computing elementary stiffness matrix
		    gmm::mult(kappa_e, inverse_J, a0);

		    gmm::mult(a0, gradPhi1_red, a1);
		    gmm::mult(a0, gradPhi2_red, a2);
		    gmm::mult(a0, gradPhi3_red, a3);

		    gmm::mult(inverse_J, gradPhi1_red, b1);
		    gmm::mult(inverse_J, gradPhi2_red, b2);
		    gmm::mult(inverse_J, gradPhi3_red, b3);

			//Computation of the elementary stiffness matrix
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
		}// end if method ==1

		else if(method == 2)
		{
			double delta_u = 1e-5;

			std::vector <double> q_int1_plus(3); // Elementary q_int vectors
			std::vector <double> q_int1_minus(3);
			std::vector <double> q_int2_plus(3);
			std::vector <double> q_int2_minus(3);
			std::vector <double> q_int3_plus(3);
			std::vector <double> q_int3_minus(3);

		    u1 = temperaturesSlave[0]+delta_u;
		    u2 = temperaturesSlave[1];
		    u3 = temperaturesSlave[2];
		    
		    T_mean = (1.0/3.0)*(u1 + u2 + u3);
		    
		    conditions_micro.meanTemperature = T_mean;
		    conditions_micro.xGradient = u1*gradPhi1[0] + u2*gradPhi2[0] + u3*gradPhi3[0];
		    conditions_micro.yGradient = u1*gradPhi1[1] + u2*gradPhi2[1] + u3*gradPhi3[1];

			fem(nodes_micro, elements_micro, physicals_micro, parameters_micro, solutionTemperature_micro, solutionFlux_micro, THERMALFLAG, PERIODICFLAG, conditions_micro, eps, type);

		    Average_flux(solutionTemperature_micro, region_micro, elements_micro, q_Me, vol_micro, thermalOrElectrical);

			// Computation of the q_int_e vector (q_int of the element) :
			q_int1_plus[0] = gmm::vect_sp(q_Me, gradPhi1);
			q_int1_plus[1] = gmm::vect_sp(q_Me, gradPhi2);
			q_int1_plus[2] = gmm::vect_sp(q_Me, gradPhi3);
		    gmm::scale(q_int1_plus,det_J/2);

		    u1 = temperaturesSlave[0]-delta_u;
		    u2 = temperaturesSlave[1];
		    u3 = temperaturesSlave[2];
		    
			//cout <<  "t1 " << u1 << endl;
			//cout <<  "t2 " << u2 << endl;
			//cout <<  "t3 " << u3 << endl;
		    //T_mean = (1.0/3.0)*(u1 + u2 + u3);
			//cout << "Tmean " << T_mean << endl;
		    
		    conditions_micro.meanTemperature = T_mean;
		    conditions_micro.xGradient = u1*gradPhi1[0] + u2*gradPhi2[0] + u3*gradPhi3[0];
		    conditions_micro.yGradient = u1*gradPhi1[1] + u2*gradPhi2[1] + u3*gradPhi3[1];
			//cout << "xgrad " << conditions_micro.xGradient << endl;
			//cout << "ygrad " << conditions_micro.yGradient << endl;

			fem(nodes_micro, elements_micro, physicals_micro, parameters_micro, solutionTemperature_micro, solutionFlux_micro, THERMALFLAG, PERIODICFLAG, conditions_micro, eps, type);

		    Average_flux(solutionTemperature_micro, region_micro, elements_micro, q_Me, vol_micro, thermalOrElectrical);

			// Computation of the q_int_e vector (q_int of the element) :
			q_int1_minus[0] = gmm::vect_sp(q_Me, gradPhi1);
			q_int1_minus[1] = gmm::vect_sp(q_Me, gradPhi2);
			q_int1_minus[2] = gmm::vect_sp(q_Me, gradPhi3);
		    gmm::scale(q_int1_minus,det_J/2);

		    u1 = temperaturesSlave[0];
		    u2 = temperaturesSlave[1]+delta_u;
		    u3 = temperaturesSlave[2];
		    
		    T_mean = (1.0/3.0)*(u1 + u2 + u3);
		    
		    conditions_micro.meanTemperature = T_mean;
		    conditions_micro.xGradient = u1*gradPhi1[0] + u2*gradPhi2[0] + u3*gradPhi3[0];
		    conditions_micro.yGradient = u1*gradPhi1[1] + u2*gradPhi2[1] + u3*gradPhi3[1];

			fem(nodes_micro, elements_micro, physicals_micro, parameters_micro, solutionTemperature_micro, solutionFlux_micro, THERMALFLAG, PERIODICFLAG, conditions_micro, eps, type);

		    Average_flux(solutionTemperature_micro, region_micro, elements_micro, q_Me, vol_micro, thermalOrElectrical);

			// Computation of the q_int_e vector (q_int of the element) :
			q_int2_plus[0] = gmm::vect_sp(q_Me, gradPhi1);
			q_int2_plus[1] = gmm::vect_sp(q_Me, gradPhi2);
			q_int2_plus[2] = gmm::vect_sp(q_Me, gradPhi3);
		    gmm::scale(q_int2_plus,det_J/2);

		    u1 = temperaturesSlave[0];
		    u2 = temperaturesSlave[1]-delta_u;
		    u3 = temperaturesSlave[2];
		    
		    T_mean = (1.0/3.0)*(u1 + u2 + u3);
		    
		    conditions_micro.meanTemperature = T_mean;
		    conditions_micro.xGradient = u1*gradPhi1[0] + u2*gradPhi2[0] + u3*gradPhi3[0];
		    conditions_micro.yGradient = u1*gradPhi1[1] + u2*gradPhi2[1] + u3*gradPhi3[1];

			fem(nodes_micro, elements_micro, physicals_micro, parameters_micro, solutionTemperature_micro, solutionFlux_micro, THERMALFLAG, PERIODICFLAG, conditions_micro, eps, type);

		    Average_flux(solutionTemperature_micro, region_micro, elements_micro, q_Me, vol_micro, thermalOrElectrical);

			// Computation of the q_int_e vector (q_int of the element) :
			q_int2_minus[0] = gmm::vect_sp(q_Me, gradPhi1);
			q_int2_minus[1] = gmm::vect_sp(q_Me, gradPhi2);
			q_int2_minus[2] = gmm::vect_sp(q_Me, gradPhi3);
		    gmm::scale(q_int2_minus,det_J/2);

		    u1 = temperaturesSlave[0];
		    u2 = temperaturesSlave[1];
		    u3 = temperaturesSlave[2]+delta_u;
		    
		    T_mean = (1.0/3.0)*(u1 + u2 + u3);
		    
		    conditions_micro.meanTemperature = T_mean;
		    conditions_micro.xGradient = u1*gradPhi1[0] + u2*gradPhi2[0] + u3*gradPhi3[0];
		    conditions_micro.yGradient = u1*gradPhi1[1] + u2*gradPhi2[1] + u3*gradPhi3[1];

			fem(nodes_micro, elements_micro, physicals_micro, parameters_micro, solutionTemperature_micro, solutionFlux_micro, THERMALFLAG, PERIODICFLAG, conditions_micro, eps, type);

		    Average_flux(solutionTemperature_micro, region_micro, elements_micro, q_Me, vol_micro, thermalOrElectrical);

			// Computation of the q_int_e vector (q_int of the element) :
			q_int3_plus[0] = gmm::vect_sp(q_Me, gradPhi1);
			q_int3_plus[1] = gmm::vect_sp(q_Me, gradPhi2);
			q_int3_plus[2] = gmm::vect_sp(q_Me, gradPhi3);
		    gmm::scale(q_int3_plus,det_J/2);

		    u1 = temperaturesSlave[0];
		    u2 = temperaturesSlave[1];
		    u3 = temperaturesSlave[2]-delta_u;
		    
		    T_mean = (1.0/3.0)*(u1 + u2 + u3);
		    
		    conditions_micro.meanTemperature = T_mean;
		    conditions_micro.xGradient = u1*gradPhi1[0] + u2*gradPhi2[0] + u3*gradPhi3[0];
		    conditions_micro.yGradient = u1*gradPhi1[1] + u2*gradPhi2[1] + u3*gradPhi3[1];

			fem(nodes_micro, elements_micro, physicals_micro, parameters_micro, solutionTemperature_micro, solutionFlux_micro, THERMALFLAG, PERIODICFLAG, conditions_micro, eps, type);

		    Average_flux(solutionTemperature_micro, region_micro, elements_micro, q_Me, vol_micro, thermalOrElectrical);

			// Computation of the q_int_e vector (q_int of the element) :
			q_int3_minus[0] = gmm::vect_sp(q_Me, gradPhi1);
			q_int3_minus[1] = gmm::vect_sp(q_Me, gradPhi2);
			q_int3_minus[2] = gmm::vect_sp(q_Me, gradPhi3);
		    gmm::scale(q_int3_minus,det_J/2);

			//Computation of the elementary stiffness matrix
		    element_stiffness(0, 0) = q_int1_plus[0]-q_int1_minus[0];
		    element_stiffness(0, 1) = q_int2_plus[0]-q_int2_minus[0];
		    element_stiffness(0, 2) = q_int3_plus[0]-q_int3_minus[0];
		    element_stiffness(1, 0) = q_int1_plus[1]-q_int1_minus[1];
		    element_stiffness(1, 1) = q_int2_plus[1]-q_int2_minus[1];
		    element_stiffness(1, 2) = q_int3_plus[1]-q_int3_minus[1];
		    element_stiffness(2, 0) = q_int1_plus[2]-q_int1_minus[2];
		    element_stiffness(2, 1) = q_int2_plus[2]-q_int2_minus[2];
		    element_stiffness(2, 2) = q_int3_plus[2]-q_int3_minus[2];

			gmm::scale(element_stiffness, 1./(2*delta_u));
		}//end if method ==2
        
    	stiffnessMaster[0] = element_stiffness(0, 0);
    	stiffnessMaster[1] = element_stiffness(0, 1);
    	stiffnessMaster[2] = element_stiffness(0, 2);
    	stiffnessMaster[3] = element_stiffness(1, 0);
    	stiffnessMaster[4] = element_stiffness(1, 1);
    	stiffnessMaster[5] = element_stiffness(1, 2);
    	stiffnessMaster[6] = element_stiffness(2, 0);
    	stiffnessMaster[7] = element_stiffness(2, 1);
        stiffnessMaster[8] = element_stiffness(2, 2);

        MPI_Send (&stiffnessMaster[0], 9, MPI_DOUBLE, 0, 39, MPI_COMM_WORLD);
		MPI_Send (&numNodesMaster[0], 3, MPI_INT, 0, 40, MPI_COMM_WORLD);
		MPI_Send (&q_int_e[0], 3, MPI_DOUBLE, 0, 41, MPI_COMM_WORLD);

    }//end while subprocesses.
}// end if myrank != 0.

MPI_Finalize();

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
}

