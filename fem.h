#ifndef FEM_H_INCLUDED
#define FEM_H_INCLUDED

#include "gmshio.h"
#include "physicalio.h"

enum FemFlag
{
    THERMALFLAG,
    ELECTRICFLAG,
    DIRICHLETFLAG,
    PERIODICFLAG,
    VONNEUMANNFLAG,
	COUPLINGFLAG,
	DEFAULTFLAG
};

struct NodeCorner
{
	NodeCorner()
	{
		C1 = NULL;
		C2 = NULL;
		C3 = NULL;
		C4 = NULL;
    }
    Node *C1;
    Node *C2;
    Node *C3;
    Node *C4;
};

struct NodeBorder
{
    std::vector<Node*> LeftNodes;
    std::vector<Node*> RightNodes;
    std::vector<Node*> TopNodes;
    std::vector<Node*> BottomNodes;
};

void fem(std::vector<Node*> &nodes, std::vector<Element*> &elements, std::vector<Physical*> &physicals, std::vector<Parameter*> &parameters, std::map<Node*, std::vector<double> > &solutionTemperature,
        std::map<Node*, std::vector<double> > &solutionFlux, FemFlag thermalOrElectrical, FemFlag method, Periodique &conditions, double eps, Type type, int JouleEffect);
void f_function(std::vector<double> &f, std::vector<Node*> &nodes, std::vector<Element*> &elements, std::map<int,Parameter*> &region, FemFlag thermalOrElectrical, int constantProperty, int JouleEffect);

void Average_flux(std::map<Node*, std::vector<double> > &solutionTemperature, std::map<int,Parameter*> &region, std::vector<Element*> &elements, std::vector<double> &q_Me, double vol, FemFlag thermalOrElectrical);

void Average_Joule(std::map<Node*, std::vector<double> > &solutionTemperature, std::vector<Element*> &elements, std::map<int, Parameter*> &region, double &s_Me, double vol);

void Average_Joule_FE1(std::map<Node*, std::vector<double> > &solutionTemperature, std::vector<Element*> &elements, std::map<int, Parameter*> &region);

#endif // FEM_H_INCLUDED
