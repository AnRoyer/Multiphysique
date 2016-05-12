#ifndef NEWTON_RAPHSON_H
#define NEWTON_RAPHSON_H

#include "fem.h"
#include "gmm/gmm.h"
#include "gmshio.h"
#include "physicalio.h"

void NewtonRaphson(std::vector<Node*> &nodes, std::vector<Element*> &elements, std::vector<Physical*> &physicals, std::vector<double> &theta_k,
					std::vector<double> &qext,FemFlag method, std::map<Node*, Node*> &NodesCorresp,
					std::vector<double> &delta_theta_k,std::map<int, Parameter*> &region,std::vector<double> &RHS, NodeCorner &corner,
					NodeBorder &border, Periodique &conditions, FemFlag thermalOrElectrical);

void Internal_flux(std::vector<double> &theta_k, std::map<int, Parameter*> &region, std::vector<Element*> &elements,
                   std::vector<double> &qint, std::vector<double> &q_m_x, std::vector<double> &q_m_y, FemFlag thermalOrElectrical);

void Tangent_Stiffness_Matrix(std::vector<double> &theta_k,std::map<int, Parameter*> &region,
                              std::vector<Element*> &elements, gmm::row_matrix< gmm::wsvector<double> >&KT,
                              std::map<Node*, Node*> &NodesCorresp,std::vector<Node*> &nodes, FemFlag thermalOrElectrical);

bool End_Criterion(std::vector<double> &RHS, double normRHS0, double eps, Type type);
#endif // NEWTON_RAPHSON_H
