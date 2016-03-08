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

void fem(std::vector<Node*> &nodes, std::vector<Element*> &elements, std::vector<Physical*> &physicals, std::vector<Parameter*> &parameters, std::map<Node*, std::vector<double> > &solution, FemFlag type, FemFlag method, Periodique &conditions);
void f_function(std::vector<double> &f, std::vector<Node*> &nodes, std::vector<Element*> &elements, std::map<int,std::vector<double> > &surfaceRegion, int constantProperty);

#endif // FEM_H_INCLUDED
