#include <iostream>
#include <cstdio>
#include <stdlib.h>
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


using namespace std;

void FE2(std::vector<Node*> &nodes_micro, std::vector<Element*> &elements_micro, std::vector<Physical*> &physicals_micro,
         std::vector<Parameter*> &parameters_micro, std::map<Node*, std::vector<double> > &solutionTemperature_micro,
         std::map<Node*, std::vector<double> > &solutionFlux_micro, Periodique &conditions_micro, std::vector<Node*> &nodes_macro,
         std::vector<Element*> &elements_macro, std::vector<Physical*> &physicals_macro,std::vector<Parameter*> &parameters_macro,
         std::map<Node*, std::vector<double> > &solutionTemperature_macro, std::map<Node*, std::vector<double> > &solutionFlux_macro,
	     Periodique &conditions_macro, double eps, std::vector<int> &methodFE2, Type type_thermic, Type type_electric,
		 int argc, char **argv, int &nature)
{

// MPI initialization
MPI_Init(&argc, &argv);
MPI_Status status;
int nbproc, myrank ;
MPI_Comm_rank( MPI_COMM_WORLD, &myrank);
MPI_Comm_size( MPI_COMM_WORLD, &nbproc);

FemFlag thermalOrElectrical = DEFAULTFLAG;
Type type = DEFAULTTYPE;
//Allows to run acoupling scheme.
int JouleEffect = 0;
int i_FE2 = 0;

if(nbproc ==1)
{
	cout << endl << "Error: FE2 method must be run with at least 2 processes !" << endl << endl;
	MPI_Finalize();
	return;
}

if (nature == 1)
{
	thermalOrElectrical = THERMALFLAG;
	type = type_thermic;
	i_FE2 = 1;
}
else if (nature == 2)
{
	thermalOrElectrical = ELECTRICFLAG;
	type = type_electric;
	i_FE2 = 1;
}
else if (nature == 3)
{
	thermalOrElectrical = ELECTRICFLAG;
	type = type_electric;
	i_FE2 = 2;
}

/*cout << "thermalOrElectrical" << thermalOrElectrical << endl;
cout << "type " << type << endl;
cout << "nature " << nature << endl;
cout << "i_FE2 " << i_FE2 << endl;*/

//Flag defining which FE2 method will be used.
int method;
if(thermalOrElectrical == THERMALFLAG) method = methodFE2[0];
else if(thermalOrElectrical == ELECTRICFLAG) method = methodFE2[1];

//Error and criterion for the FE2 method.
std::vector<double> error(nodes_macro.size());
double criterionFE2 = 1;
double criterionFE2_0 =0;
double criterionFE2_old = 0;
double criterionFE2_min = 1e-6; //Is initiated to 1e-6 but will automatically grow if method 2 is not converging.

double u1, u2, u3; // Nodal temperatures of an element.
double T_mean; // Mean temperature of an element.
std::vector<double> gradT(2);// Mean gradient over an element.
std::vector<double> q_Me(2);   // Mean flux over an element.
double s_Me;
gmm::dense_matrix<double> kappa_e(2,2);// Conductivity of an element.
gmm::dense_matrix<double> element_stiffness(3, 3); // Stiffness matrix of an element.
std::vector<double> q_int_e(3); // Elementary q_int vector
int elementFlag = 0;

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

gmm::dense_matrix<double> a0(2, 2); // kappa_e_e times inverse of J in Eq. (9) of ddl 4 group B
std::vector<double> a1(2); // a1 times gradPhi_1
std::vector<double> a2(2); // a1 times gradPhi_2
std::vector<double> a3(2); // a1 times gradPhi_3
std::vector<double> b1(2); // inverse of J times gradPhi_i
std::vector<double> b2(2); // inverse of J times gradPhi_2
std::vector<double> b3(2); // inverse of J times gradPhi_3

gmm::dense_matrix<double> total_stiffness(nodes_macro.size(), nodes_macro.size()); // Matrix of the stiffness for the whole domain K
std::vector <double> q_int(nodes_macro.size());

std::map<int, Parameter*> region_micro;//Stock le lien entre le numéro du physical de msh (stocker dans "physicals") et la valeur du parametre de "parametres" pour les régions de dimension 1 (ligne)

//Loading parameters.
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


//Computation of the domains volumes.
double vol_micro=0;
std::vector<double> c(nodes_micro.size());

f_function(c, nodes_micro, elements_micro, region_micro, thermalOrElectrical, 1, 0); //FemFlag is obsolete !

for(unsigned int j=0; j<nodes_micro.size(); j++)
{
    vol_micro += c[j];
}
double vol_macro=0;
std::vector<double> c2(nodes_macro.size());

f_function(c2, nodes_macro, elements_macro, region_macro, thermalOrElectrical, 1, 0); //FemFlag is obsolete !

for(unsigned int j=0; j<nodes_macro.size(); j++)
{
    vol_macro += c2[j];
}

//Computation of the f_i vector function of the heat generation in the macro domain.
std::vector<double> f_i(nodes_macro.size());
//f_function(f_i, nodes_macro, elements_macro, region_macro, thermalOrElectrical, 0, JouleEffect); //utilisation of f_function from fem.cpp file.

//-----------------Variables from the parallel loop-------------------------------

int nbElem = elements_macro.size();
std::vector<double> sol_u_tmp;
int source;	// Rank of the process sending a submatrix Ke
int numToSend; // Number of the element to send to the slave. The slave will be responsible for this macroscopic element.
int numElem; // Contains the number of the finite element associated to the Ke matrix a slave sends.
std::vector<int> numNodesMaster(3); // Number of the nodes sent by the slave to the master
int elementNumber;
std::vector<double> temperaturesSlave(3); // Temperatures of the nodes sent to the slave
std::vector<double> stiffnessMaster(9); // Will contain the 9 elements of the submatrice Ke (sent by process source)
//if(method == 2) gmm::scale(f_i,-1);//due to sign convention differing in the two methods.
int k;
int OkFE2 = 0;
if(myrank == 0)// process 0 will take care of all the displaying.
{
	while(i_FE2>0) //Allows coupling.
	{
		if(nature == 3) //if coupling.
		{
			if(i_FE2 == 2 && type == FE2withDIRICHLET) OkFE2 = 1;
			else
			{
				if(i_FE2 == 1 && type == FE2withDIRICHLET) OkFE2 = 1;
				else OkFE2 = 0;
			}
		}
		else if(nature == 1 || nature == 2) OkFE2 =1; //if not coupling.

		if(OkFE2 == 1)
		{
			//Displaying information
			if(thermalOrElectrical == ELECTRICFLAG)
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
				cout << "\t##     FFFFFFFFFFFF        EEEEEEEEEEEE             ,/    ##" << endl;
				cout << "\t##     FF                  EE                     ,'/     ##" << endl;
				cout << "\t##     FF                  EE                   ,' /      ##" << endl;
				cout << "\t##     FF                  EE                 ,'  /_____, ##" << endl;
				cout << "\t##     FF                  EE               .'____    ,'  ##" << endl;
				cout << "\t##     FF                  EEEEEEEEEEEEEEEE      /  ,'    ##" << endl;
				cout << "\t##                                              / ,'      ##" << endl;
				cout << "\t###############################################/,' #########" << endl;
				cout << "\t##############################################/'############" << endl << endl << endl;
			}
			else if(thermalOrElectrical == THERMALFLAG)
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
				cout << "\t##     FF                  EE                     (       ##" << endl;
				cout << "\t##     FFFFFFFFFFFF        EEEEEEEEEEEE           ))      ##" << endl;
				cout << "\t##     FF                  EE                     {_}     ##" << endl;
				cout << "\t##     FF                  EE                    .-;-.    ##" << endl;
				cout << "\t##     FF                  EE                   |'-=-'|   ##" << endl;
				cout << "\t##     FF                  EE                   |     |   ##" << endl;
				cout << "\t##     FF                  EEEEEEEEEEEEEEEE     |     |   ##" << endl;
				cout << "\t##                                              |     |   ##" << endl;
				cout << "\t################################################|     |#####" << endl;
				cout << "\t################################################'.___.'#####" << endl << endl << endl;
			}

			//Affichage des infos principales

			//if(thermalOrElectrical == THERMALFLAG) cout << "SOLVING A THERMAL PROBLEM." << endl;
			//if(thermalOrElectrical == ELECTRICFLAG) cout << "SOLVING AN ELECTRICAL PROBLEM." << endl;
			if(type == FE2withDIRICHLET) cout << "Calling the FE2 method with DIRICHLET conditions." << endl;
			if(type == FE2withVONNEUMANN) cout << "Calling the FE2 method with VONNEUMANN conditions." << endl;
			if(type == FE2withPERIODIC) cout << "Calling the FE2 method with PERIODIC conditions." << endl;
			cout << endl;

			if(method ==1) cout << "Using method 1 for assembling the stiffness matrix." << endl;
			else if(method ==2) cout << "Using method 2 for assembling the stiffness matrix." << endl;
			cout << endl;
			cout << "-----------------------------" << endl;
			// end of displaying (master)

			//Initial guess.
			if(type == FE2withDIRICHLET)
			{
				fem(nodes_macro, elements_macro, physicals_macro, parameters_macro, solutionTemperature_macro, solutionFlux_macro, thermalOrElectrical, DIRICHLETFLAG, conditions_macro, eps, type, JouleEffect);
				/*for(int k=0;k<nodes_macro.size();k++)
				{
					cout << solutionTemperature_macro[nodes_macro[k]] << " " ;
				}
				cout << endl;*/
			}
			else if(type == FE2withVONNEUMANN)
			{
				fem(nodes_macro, elements_macro, physicals_macro, parameters_macro, solutionTemperature_macro, solutionFlux_macro, thermalOrElectrical, VONNEUMANNFLAG, conditions_macro, eps, type, JouleEffect);
			}
			/*for(int k=0;k<nodes_macro.size();k++)
			{
				cout << solutionTemperature_macro[nodes_macro[k]] << endl;
			}
			cout << endl;*/

			/*elementFlag = -4;
			for(unsigned int i=1; i<nbproc;i++)
			{
				MPI_Send (&elementFlag, 1, MPI_INT, i, 31, MPI_COMM_WORLD);
				//This will tell the subprocesses to make an initial guess too.
			}*/

			//int flagFlux = 0;
			int i_while = 0;
			while(criterionFE2 > criterionFE2_min)
			{
				cout << "Running Newton Raphson on the RVEs for macroscopic iteration " << i_while+1 << " ..." << endl;

				gmm::clear(total_stiffness);
				gmm::clear(f_i);
				f_function(f_i, nodes_macro, elements_macro, region_macro, thermalOrElectrical, 0, JouleEffect); //utilisation of f_function from fem.cpp file.
				if(method == 2) gmm::scale(f_i,-1);
				//if(i_while == 1) cout << "f_i" << f_i << endl;

				std::vector <double> q_int(nodes_macro.size());
				// INITIALIZATION : Now an element number is sent to each slave, so that all the slaves are busy at the beginning.
				int p; // Process number
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
					MPI_Send (& numToSend, 1, MPI_INT, p, 31, MPI_COMM_WORLD);
					MPI_Send (&temperaturesSlave[0], 3, MPI_DOUBLE, p, 32, MPI_COMM_WORLD);
				}
				// End of initialization

				/* Now the master loops over the remainder of the finite elements. When the slave process number source sends back a Ke matrix,
				master process asks it to work now with finite element i and put the submatrix in K.
				*/

				for (int i = nbproc - 1+k; i < nbElem + nbproc-1; i++)
				{
					//cout << "i " << i << endl;
					if(i>=(nbElem-k-1)/5+k && i<(nbElem-k-1)/5 +1+k) cout << " 20% Done..." << endl;
					else if(i>=2*(nbElem-k-1)/5+k && i<2*(nbElem-k-1)/5 +1+k) cout << " 40% Done..." << endl;
					else if(i>=3*(nbElem-k-1)/5+k && i<3*(nbElem-k-1)/5 +1+k) cout << " 60% Done..." << endl;
					else if(i>=4*(nbElem-k-1)/5+k && i<4*(nbElem-k-1)/5 +1+k) cout << " 80% Done..." << endl;
					else if(i>=(nbElem-k-1)+k && i<(nbElem-k-1)+1+k) cout << " 100% Done. Macroscopic iteration " << i_while +1 << " is finished." << endl;

					// Reception d'une sous matrice de la part d'un process a priori inconnu :
					MPI_Recv (&stiffnessMaster[0], 9, MPI_DOUBLE, MPI_ANY_SOURCE, 33, MPI_COMM_WORLD, & status);
					source = status.MPI_SOURCE;

					// Réception des noeuds correspondants
					MPI_Recv (&numNodesMaster[0], 3, MPI_INT, source, 34, MPI_COMM_WORLD, & status);
					MPI_Recv (&numElem, 1, MPI_INT, source, 38, MPI_COMM_WORLD, & status);

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

					//Assemling q_int :
					MPI_Recv (&q_int_e[0], 3, MPI_DOUBLE, source, 35, MPI_COMM_WORLD, & status);

					q_int[numNodesMaster[0]] = q_int[numNodesMaster[0]] + q_int_e[0];
					q_int[numNodesMaster[1]] = q_int[numNodesMaster[1]] + q_int_e[1];
					q_int[numNodesMaster[2]] = q_int[numNodesMaster[2]] + q_int_e[2];

					//Computation of the macroscopic flux using q_Me
					MPI_Recv (&q_Me[0], 2, MPI_DOUBLE, source, 36, MPI_COMM_WORLD, & status);
					elements_macro[numElem]->q_Mex = q_Me[0];
					elements_macro[numElem]->q_Mey = q_Me[1];
					std::vector<double> val2(3);
					val2[0] = q_Me[0];
					val2[1] = q_Me[1];
					val2[2] = 0;
					solutionFlux_macro[nodes_macro[numNodesMaster[0]]] = val2;
					solutionFlux_macro[nodes_macro[numNodesMaster[1]]] = val2;
					solutionFlux_macro[nodes_macro[numNodesMaster[2]]] = val2;

					MPI_Recv (&s_Me, 1, MPI_DOUBLE, source, 37, MPI_COMM_WORLD, & status);
					if(thermalOrElectrical == ELECTRICFLAG && nature == 3) elements_macro[i-nbproc+1]->s_Me = s_Me;

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

						MPI_Send (&numToSend, 1, MPI_INT, source, 31, MPI_COMM_WORLD);
						MPI_Send (&temperaturesSlave[0], 3, MPI_DOUBLE, source, 32, MPI_COMM_WORLD);
					}

				}//end for.

				for(unsigned int i = 0; i < elements_macro.size(); i++)
				{
					if(elements_macro[i]->type == 1)//If line
					{
						if(region_macro.count(elements_macro[i]->region) == 1 && (region_macro[elements_macro[i]->region]->temperature != -1 && thermalOrElectrical == THERMALFLAG || region_macro[elements_macro[i]->region]->voltage != -1 && thermalOrElectrical == ELECTRICFLAG))
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

				//Part used to compute the energy balance in the macroscopic domain.
				//cout << "q_int " << q_int << endl;
				FILE *fp2 = fopen("q_int.dat", "w");
				for(int l = 0; l < q_int.size(); l++)
				fprintf(fp2, "%.15f \n", q_int[l]);
				fclose(fp2);

				fp2 = fopen("x.dat", "w");
				for(int l = 0; l < nodes_macro.size(); l++)
				fprintf(fp2, "%.15f \n", nodes_macro[l]->x);
				fclose(fp2);

				fp2 = fopen("y.dat", "w");
				for(int l = 0; l < nodes_macro.size(); l++)
				fprintf(fp2, "%.15f \n", nodes_macro[l]->y);
				fclose(fp2);

				//cout << "total_stiffness " << total_stiffness << endl;
				error.clear();
				for(unsigned int i=0;i<nodes_macro.size();i++)
				{
					error.push_back(q_int[i] - f_i[i]);
				}


				if(thermalOrElectrical == THERMALFLAG)
				{
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
				}

				else if(thermalOrElectrical == ELECTRICFLAG)
				{
					for(unsigned int i = 0; i < elements_macro.size(); i++)
					{
						if(elements_macro[i]->type == 1 && region_macro[elements_macro[i]->region]->voltage != -1)//If line
						{
							if(region_macro.count(elements_macro[i]->region) == 1)//If linesRegion contains elements[i]->region
							{
								for(unsigned int j = 0; j < elements_macro[i]->nodes.size(); j++)
								{
								    f_i[elements_macro[i]->nodes[j]->num-1] = region_macro[elements_macro[i]->region]->voltage;
									error[elements_macro[i]->nodes[j]->num-1] = 0.;
								}
							}

						}
					}
				}

				//cout << "f_i cond" << f_i << endl;
				std::vector<double> u_guess_vec(1);
				std::vector<double> u_guess(nodes_macro.size());

				for(unsigned int l=0; l<nodes_macro.size(); l++)
				{
					u_guess_vec = solutionTemperature_macro[nodes_macro[l]];
					u_guess[l] = u_guess_vec[0];
				}

				if(i_while ==0)
				{
					criterionFE2 = 1000;
					criterionFE2_0 = gmm::vect_norm2(error);
				}

				if(i_while !=0) criterionFE2_old = criterionFE2;
				if(i_while ==0) criterionFE2_old = 2000;
				criterionFE2 = gmm::vect_norm2(error);

				if(criterionFE2_0 > criterionFE2_min)
				{
					criterionFE2 = criterionFE2/criterionFE2_0;
				}

				cout << "FE2 relative residue = " << criterionFE2 << endl << endl;

				if(criterionFE2 > criterionFE2_min)
				{
					//cout << "ok solving !" << endl;
					#ifdef GMM_USES_MUMPS
						//std::cout << "solving linear system with MUMPS\n";
						gmm::csr_matrix<double> total_stiffness_csr; //à changer
						gmm::copy(total_stiffness,total_stiffness_csr);
						gmm::MUMPS_solve(total_stiffness_csr, u_guess, f_i);
					#else
						std::cout << "solving linear system with gmm::lu_solve !\n";
						gmm::lu_solve(total_stiffness, u_guess, f_i);
					#endif

					for(unsigned int l=0; l<nodes_macro.size(); l++)
					{
						u_guess_vec[0] = u_guess[l];
						solutionTemperature_macro[nodes_macro[l]] = u_guess_vec;
					}
				}

				/*if(i_while == 1)
					{
					for(int k=0;k<nodes_macro.size();k++)
					{
						cout << solutionTemperature_macro[nodes_macro[k]] << endl;
					}
					cout << endl;
				}*/

				i_while ++;

				//This section takes care of the special cases when the criterion is not converging
				//-----------------------------------------------------------------------------------------

				if(method == 2 && criterionFE2>criterionFE2_old && criterionFE2 > criterionFE2_min)
				{
					if(i_while == 0) cout << "!!! Error: probably due to non-linearity in the microscopic domain !!!" << endl;
					else cout << "The method stopped converging but couldn't reach the minimal criterion (" << criterionFE2_min << ")." << endl;
					while(criterionFE2 > criterionFE2_min) criterionFE2_min = 10*criterionFE2_min;
					//cout << "!!! The minimal criterion has been rised to " << criterionFE2_min << "!!!" << endl;
					cout << "Press enter to continue." << endl;
		  			std::cin.ignore();
				}

				if(method ==1 && ((criterionFE2 > 0.01 && i_while == 3) || (criterionFE2 > 0.001 && i_while == 10) ))//To avoid looping when not converging.
				{
					cout << endl;
					cout << "Method 1 is a very bad method :-(..." << endl;
					i_while = 0;
					cout << "Press enter to try method 2.";
					cout << endl;
		  			cin.ignore();
					method = 2;//Automatically switches to method 2.
					elementFlag = -2;
					for(unsigned int i=1; i<nbproc;i++)
					{
						MPI_Send (&elementFlag, 1, MPI_INT, i, 31, MPI_COMM_WORLD);
						//This will tell the subprocesses they have to switch to method 2.
		   			}
					cout << "-----------------------------" << endl;
					//Initial guess.
					if(type == FE2withDIRICHLET)
					{
						fem(nodes_macro, elements_macro, physicals_macro, parameters_macro, solutionTemperature_macro, solutionFlux_macro, thermalOrElectrical, DIRICHLETFLAG, conditions_macro, eps, type, JouleEffect);
						/*for(int k=0;k<nodes_macro.size();k++)
						{
							cout << solutionTemperature_macro[nodes_macro[k]] << " " ;
						}
						cout << endl;*/
					}
					else if(type == FE2withVONNEUMANN)
					{
						fem(nodes_macro, elements_macro, physicals_macro, parameters_macro, solutionTemperature_macro, solutionFlux_macro, thermalOrElectrical, VONNEUMANNFLAG, conditions_macro, eps, type, JouleEffect);
					}
				}//end if

				//-----------------------------------------------------------------------------------------

			}//end while

		}//end if OkFE2.

		else if(OkFE2 == 0)//If no fe2D needed.
		{

			if(thermalOrElectrical == ELECTRICFLAG)
			{
				cout << endl;
				cout << "\t############################################################" << endl;
				cout << "\t############################################################" << endl;
				cout << "\t##                                                        ##" << endl;
				cout << "\t##                                                 11     ##" << endl;
				cout << "\t##                                               1111     ##" << endl;
				cout << "\t##     FFFFFFFFFFFFFFFF    EEEEEEEEEEEEEEEE     11 11     ##" << endl;
				cout << "\t##     FF                  EE                      11     ##" << endl;
				cout << "\t##     FF                  EE                      11     ##" << endl;
				cout << "\t##     FF                  EE                             ##" << endl;
				cout << "\t##     FF                  EE                             ##" << endl;
				cout << "\t##     FFFFFFFFFFFF        EEEEEEEEEEEE             ,/    ##" << endl;
				cout << "\t##     FF                  EE                     ,'/     ##" << endl;
				cout << "\t##     FF                  EE                   ,' /      ##" << endl;
				cout << "\t##     FF                  EE                 ,'  /_____, ##" << endl;
				cout << "\t##     FF                  EE               .'____    ,'  ##" << endl;
				cout << "\t##     FF                  EEEEEEEEEEEEEEEE      /  ,'    ##" << endl;
				cout << "\t##                                              / ,'      ##" << endl;
				cout << "\t###############################################/,' #########" << endl;
				cout << "\t##############################################/'############" << endl << endl << endl;
			}
			else if(thermalOrElectrical == THERMALFLAG)
			{
				cout << endl;
				cout << "\t############################################################" << endl;
				cout << "\t############################################################" << endl;
				cout << "\t##                                                        ##" << endl;
				cout << "\t##                                                 11     ##" << endl;
				cout << "\t##                                               1111     ##" << endl;
				cout << "\t##     FFFFFFFFFFFFFFFF    EEEEEEEEEEEEEEEE     11 11     ##" << endl;
				cout << "\t##     FF                  EE                      11     ##" << endl;
				cout << "\t##     FF                  EE                      11     ##" << endl;
				cout << "\t##     FF                  EE                             ##" << endl;
				cout << "\t##     FF                  EE                     (       ##" << endl;
				cout << "\t##     FFFFFFFFFFFF        EEEEEEEEEEEE           ))      ##" << endl;
				cout << "\t##     FF                  EE                     {_}     ##" << endl;
				cout << "\t##     FF                  EE                    .-;-.    ##" << endl;
				cout << "\t##     FF                  EE                   |'-=-'|   ##" << endl;
				cout << "\t##     FF                  EE                   |     |   ##" << endl;
				cout << "\t##     FF                  EEEEEEEEEEEEEEEE     |     |   ##" << endl;
				cout << "\t##                                              |     |   ##" << endl;
				cout << "\t################################################|     |#####" << endl;
				cout << "\t################################################'.___.'#####" << endl << endl << endl;
			}

			cout << endl;
			//if(thermalOrElectrical == THERMALFLAG) cout << "SOLVING A THERMIC PROBLEM." << endl << endl;
			//else if(thermalOrElectrical == ELECTRICFLAG) cout << "SOLVING AN ELECTRIC PROBLEM." << endl << endl;
			//cout << "Read " << nodes.size() << " nodes and " << elements.size() << " elements." << endl << endl;
			if(type == DIRICHLET) cout << "Calling the FE1 method with DIRICHLET conditions." << endl;
			else if(type == VONNEUMANN) cout << "Calling the FE1 method with VONNEUMANN conditions." << endl;
			else if(type == PERIODIC) cout << "Calling the FE1 method with PERIODIC conditions." << endl;
			cout << endl;
			cout << "-----------------------------" << endl;

			if(type == PERIODIC)
			fem(nodes_macro, elements_macro, physicals_macro, parameters_macro, solutionTemperature_macro, solutionFlux_macro, thermalOrElectrical, PERIODICFLAG, conditions_macro, eps, type, JouleEffect);

			else if(type == DIRICHLET)
			fem(nodes_macro, elements_macro, physicals_macro, parameters_macro, solutionTemperature_macro, solutionFlux_macro, thermalOrElectrical, DIRICHLETFLAG, conditions_macro, eps, type, JouleEffect);

			else if(type == VONNEUMANN)
			fem(nodes_macro, elements_macro, physicals_macro, parameters_macro, solutionTemperature_macro, solutionFlux_macro, thermalOrElectrical, VONNEUMANNFLAG, conditions_macro, eps, type, JouleEffect);

			/*for(int e=0;e<elements_macro.size();e++)
			{
				cout << elements_macro[e]->s_Me;
			}*/
			Average_Joule_FE1(solutionTemperature_macro,elements_macro,region_macro);
			/*for(int e=0;e<elements_macro.size();e++)
			{
				cout << elements_macro[e]->s_Me;
			}
			cout<< endl<< endl;*/

		} //end if OkFe2 ==0.

		//writing the solution.
		if(thermalOrElectrical == THERMALFLAG)
		{
			writeMSH((char*)"solutionTemperature.pos", solutionTemperature_macro);
			writeMSH((char*)"solutionFlux.pos", solutionFlux_macro);

			//Write in .dat file
			FILE *fp = fopen("dataMatlabT.dat", "w");
			std::map<Node*, std::vector<double> >::iterator itT = solutionTemperature_macro.begin();
			for(itT = solutionTemperature_macro.begin(); itT != solutionTemperature_macro.end(); itT++)
			fprintf(fp, "%.15f \t %.15f \t %.15f \n", itT->first->x, itT->first->y, itT->second[0]);

			fclose(fp);
			//cout << "q_int " << q_int << endl;

		}
		else if(thermalOrElectrical == ELECTRICFLAG)
		{
			writeMSH((char*)"solutionPotential.pos", solutionTemperature_macro);
			writeMSH((char*)"solutionCurrent.pos", solutionFlux_macro);

			//Write in .dat file
			FILE *fp = fopen("dataMatlabV.dat", "w");
			std::map<Node*, std::vector<double> >::iterator itT = solutionTemperature_macro.begin();
			for(itT = solutionTemperature_macro.begin(); itT != solutionTemperature_macro.end(); itT++)
			fprintf(fp, "%.15f \t %.15f \t %.15f \n", itT->first->x, itT->first->y, itT->second[0]);

			fclose(fp);
		}

		if(i_FE2 == 2) //if first part of coupling scheme: writes the Joule heat generation.
		{
			JouleEffect = 1;
			thermalOrElectrical = THERMALFLAG;
			type = type_thermic;
			method = methodFE2[0];
			criterionFE2 = 1;
			criterionFE2_0 =0;
			criterionFE2_old = 0;
			cout << endl;
			cout << "ELECTRICAL COMPUTATION FINISHED, STARTING THERMAL COMPUTATION..." << endl;
			if (type == PERIODIC || type == DIRICHLET || type == VONNEUMANN)
			{
				elementFlag = -1; //flag used to stop the subprocesse from waiting when the FE2 method is over.
				for(unsigned int i=1; i<nbproc;i++)
				{
					MPI_Send (&elementFlag, 1, MPI_INT, i, 31, MPI_COMM_WORLD);
					//cout << "sent !" << endl;
				}
			}
			else if(type == FE2withPERIODIC || type == FE2withDIRICHLET || type == FE2withVONNEUMANN)
			{
				elementFlag = -3; //flag used to tell the subrocesses to switch to thermic.
				for(unsigned int i=1; i<nbproc;i++)
				{
					MPI_Send (&elementFlag, 1, MPI_INT, i, 31, MPI_COMM_WORLD);
				}
			}
		}

		if(i_FE2 == 1) //if second part of coupling or normal implementation, finishes.
		{
			if (type == FE2withPERIODIC || type == FE2withDIRICHLET || type == FE2withVONNEUMANN)
			{
				elementFlag = -1; //flag used to stop the subprocesse from waiting when the FE2 method is over.
				for(unsigned int i=1; i<nbproc;i++)
				{
					MPI_Send (&elementFlag, 1, MPI_INT, i, 31, MPI_COMM_WORLD);
					//cout << "sent !" << endl;
				}
			}

			if(nature == 1 || nature == 2) //If not in a coupling scheme.
			{
				cout << "-----------------------------" << endl;
				if(type == FE2withDIRICHLET && thermalOrElectrical == THERMALFLAG)
				cout << "The THERMIC problem has been solved in TWO SCALES with DIRICHLET conditions." << endl;
				if(type == FE2withVONNEUMANN && thermalOrElectrical == THERMALFLAG)
				cout <<"The THERMIC problem has been solved in TWO SCALES with VON NEUMANN conditions." << endl;
				if(type == FE2withPERIODIC && thermalOrElectrical == THERMALFLAG)
				cout <<  "The THERMIC problem has been solved in TWO SCALES with PERIODIC conditions."  << endl;
				if(type == FE2withDIRICHLET && thermalOrElectrical == ELECTRICFLAG)
				cout << "The ELECTRIC problem has been solved in TWO SCALES with DIRICHLET conditions." << endl;
				if(type == FE2withVONNEUMANN && thermalOrElectrical == ELECTRICFLAG)
				cout <<"The ELECTRIC problem has been solved in TWO SCALES with VON NEUMANN conditions." << endl;
				if(type == FE2withPERIODIC && thermalOrElectrical == ELECTRICFLAG)
				cout <<  "The ELECTRIC problem has been solved in TWO SCALES with PERIODIC conditions."  << endl;

				cout << endl;

				elementFlag = -1; //flag used to stop the subprocesse from waiting when the FE2 method is over.
				for(unsigned int i=1; i<nbproc;i++)
				{
					MPI_Send (&elementFlag, 1, MPI_INT, i, 31, MPI_COMM_WORLD);
					//cout << "sent !" << endl;
				}
			}

			else if(nature == 3) //If in a coupling scheme.
			{
				cout << "-----------------------------" << endl;
				cout <<  "The COUPLED problem has been solved in TWO SCALES."  << endl;
				cout << endl;
			}
		}//end if i_FE2 == 1.

		i_FE2 --;

	}//end while i_FE2>0.
}//end master

if(myrank !=0)
{
    while(1)
    {
    	MPI_Recv (&numElem, 1, MPI_INT, 0, 31, MPI_COMM_WORLD, & status);
		//cout << "numelem" << numElem << endl;
        if(numElem == -1) break;
		else if(numElem == -2)
		{
			method = 2;
			MPI_Recv (&numElem, 1, MPI_INT, 0, 31, MPI_COMM_WORLD, & status);
		}
		else if(numElem == -3)
		{
			thermalOrElectrical = THERMALFLAG;
			type = type_thermic;
			method = methodFE2[0];
			MPI_Recv (&numElem, 1, MPI_INT, 0, 31, MPI_COMM_WORLD, & status);
		}


		MPI_Recv (&temperaturesSlave[0], 3, MPI_DOUBLE, 0, 32, MPI_COMM_WORLD, & status);
		//cout << temperaturesSlave << endl;
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

		fem(nodes_micro, elements_micro, physicals_micro, parameters_micro, solutionTemperature_micro, solutionFlux_micro, thermalOrElectrical, PERIODICFLAG, conditions_micro, eps, type, 0);

		//Computation of average flux over RVE. Will be used in method one and for computing the flux in the macroscopic domain.
		q_Me[0] = 0.0;
		q_Me[1] = 0.0;
		Average_flux(solutionTemperature_micro, region_micro, elements_micro, q_Me, vol_micro, thermalOrElectrical);

		// Computation of the q_int_e vector (q_int of the element) :
		q_int_e[0] = gmm::vect_sp(q_Me, gradPhi1);
		q_int_e[1] = gmm::vect_sp(q_Me, gradPhi2);
		q_int_e[2] = gmm::vect_sp(q_Me, gradPhi3);
	    gmm::scale(q_int_e,det_J/2);

		//If coupling, computation of s_Me
		if(nature ==3) Average_Joule(solutionTemperature_micro, elements_micro, region_micro, s_Me, vol_micro);

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
			double delta_u = 1e-1;

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

			fem(nodes_micro, elements_micro, physicals_micro, parameters_micro, solutionTemperature_micro, solutionFlux_micro, thermalOrElectrical, PERIODICFLAG, conditions_micro, eps, type, 0);

		    Average_flux(solutionTemperature_micro, region_micro, elements_micro, q_Me, vol_micro, thermalOrElectrical);

			// Computation of the q_int_e vector (q_int of the element) :
			q_int1_plus[0] = gmm::vect_sp(q_Me, gradPhi1);
			q_int1_plus[1] = gmm::vect_sp(q_Me, gradPhi2);
			q_int1_plus[2] = gmm::vect_sp(q_Me, gradPhi3);
		    gmm::scale(q_int1_plus,det_J/2);

		    u1 = temperaturesSlave[0]-delta_u;
		    u2 = temperaturesSlave[1];
		    u3 = temperaturesSlave[2];

		    conditions_micro.meanTemperature = T_mean;
		    conditions_micro.xGradient = u1*gradPhi1[0] + u2*gradPhi2[0] + u3*gradPhi3[0];
		    conditions_micro.yGradient = u1*gradPhi1[1] + u2*gradPhi2[1] + u3*gradPhi3[1];

			fem(nodes_micro, elements_micro, physicals_micro, parameters_micro, solutionTemperature_micro, solutionFlux_micro, thermalOrElectrical, PERIODICFLAG, conditions_micro, eps, type, 0);

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

			fem(nodes_micro, elements_micro, physicals_micro, parameters_micro, solutionTemperature_micro, solutionFlux_micro, thermalOrElectrical, PERIODICFLAG, conditions_micro, eps, type, 0);

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

			fem(nodes_micro, elements_micro, physicals_micro, parameters_micro, solutionTemperature_micro, solutionFlux_micro, thermalOrElectrical, PERIODICFLAG, conditions_micro, eps, type, 0);

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

			fem(nodes_micro, elements_micro, physicals_micro, parameters_micro, solutionTemperature_micro, solutionFlux_micro, thermalOrElectrical, PERIODICFLAG, conditions_micro, eps, type, 0);

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

			fem(nodes_micro, elements_micro, physicals_micro, parameters_micro, solutionTemperature_micro, solutionFlux_micro, thermalOrElectrical, PERIODICFLAG, conditions_micro, eps, type, 0);

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

		/*cout << numElem << " " << endl;
		cout << stiffnessMaster << endl;
		cout << numNodesMaster << endl;
		cout << q_int_e << endl;
		cout << q_Me << endl;*/

        MPI_Send (&stiffnessMaster[0], 9, MPI_DOUBLE, 0, 33, MPI_COMM_WORLD);
		MPI_Send (&numNodesMaster[0], 3, MPI_INT, 0, 34, MPI_COMM_WORLD);
		MPI_Send (&q_int_e[0], 3, MPI_DOUBLE, 0, 35, MPI_COMM_WORLD);
		MPI_Send (&q_Me[0], 2, MPI_DOUBLE, 0, 36, MPI_COMM_WORLD);
		MPI_Send (&s_Me, 1, MPI_DOUBLE, 0, 37, MPI_COMM_WORLD);
		MPI_Send (&numElem, 1, MPI_INT, 0, 38, MPI_COMM_WORLD);

    }//end while subprocesses.
}// end if myrank != 0.

//cout << "myrank " << myrank <<endl;
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

