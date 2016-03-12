#ifndef NEWTON_RAPHSON_H
#define NEWTON_RAPHSON_H


#include "gmm/gmm.h"
#include "gmshio.h"
#include "physicalio.h"

void NewtonRaphson(std::vector<Node*> &nodes, std::vector<Element*> &elements, std::vector<Physical*> &physicals, std::vector<Parameter*> &parameters,
					std::map<Node*, std::vector<double> > &solution, std::vector<double> &theta_k,gmm::row_matrix< gmm::wsvector<double> > Tmp,
					std::vector<double> &qext,FemFlag method,map<int, double> linesRegion, map<Node*, Node*> &NodesCorresp,
					std::vector<double> &delta_theta_k);

void Internal_flux(std::vector<double> &theta_k, std::map<int,std::vector<double> > &surfaceRegion, std::vector<Element*> &elements, std::vector<double> &qint);

void Tangent_Stiffness_Matrix(std::vector<double> &theta_k,std::map<int,std::vector<double> > &surfaceRegion, std::vector<Element*> &elements, gmm::csr_matrix<double> &KT,map<Node*, Node*> &NodesCorresp);

bool End_Criterion(std::vector<double> &theta_k, std::vector<double> &delta_theta_k);
#endif // NEWTON_RAPHSON_H
