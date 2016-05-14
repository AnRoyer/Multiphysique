#ifndef FE2_H
#define FE2_H
#include <mpi.h>
#include "fem.h"
#include "gmm/gmm.h"
#include "gmshio.h"
#include "physicalio.h"

void FE2(std::vector<Node*> &nodes_micro, std::vector<Element*> &elements_micro, std::vector<Physical*> &physicals_micro,
         std::vector<Parameter*> &parameters_micro, std::map<Node*, std::vector<double> > &solutionTemperature_micro,
         std::map<Node*, std::vector<double> > &solutionFlux_micro, Periodique &conditions_micro, std::vector<Node*> &nodes_macro,
         std::vector<Element*> &elements_macro,std::vector<Physical*> &physicals_macro,std::vector<Parameter*> &parameters_macro,
         std::map<Node*, std::vector<double> > &solutionTemperature_macro, std::map<Node*, std::vector<double> > &solutionFlux_macro,
	     Periodique &conditions_macro, double eps, std::vector<int> &methodFE2, Type type_thermic, Type type_electric,
		 int argc, char **argv, int &natureFlag);

void conductivityTensor(std::vector<double> &q, std::vector<double> &gradT, gmm::dense_matrix<double> &kappa);

#endif // FE2_H
