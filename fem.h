#ifndef FEM_H_INCLUDED
#define FEM_H_INCLUDED

#include "gmshio.h"
#include "physicalio.h"

enum FemFlag
{
    THERMALFLAG,
    ELECTRICFLAG,
    DIRICHLETFLAG,
    PERIODICFLAG
};

struct NodeCorner
{
    Node *C1 = NULL;
    Node *C2 = NULL;
    Node *C3 = NULL;
    Node *C4 = NULL;
};

struct NodeBorder
{
    std::vector<Node*> LeftNodes;
    std::vector<Node*> RightNodes;
    std::vector<Node*> TopNodes;
    std::vector<Node*> BottomNodes;
};

void fem(std::vector<Node*> &nodes, std::vector<Element*> &elements, std::vector<Physical*> &physicals, std::vector<Parameter*> &parameters, std::map<Node*, std::vector<double> > &solutionTemperature, std::map<Node*, std::vector<double> > &solutionFlux, FemFlag type, FemFlag method, Periodique &conditions);
void f_function(std::vector<double> &f, std::vector<Node*> &nodes, std::vector<Element*> &elements, std::map<int,Parameter*> &region, FemFlag type, int constantProperty);

#endif // FEM_H_INCLUDED
