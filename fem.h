#ifndef FEM_H_INCLUDED
#define FEM_H_INCLUDED

#include "gmshio.h"
#include "physicalio.h"

void f_function(std::vector<double> &f, std::vector<Node*> &nodes, std::vector<Element*> &elements, std::map<int,std::vector<double> > &surfaceRegion, int constantProperty);
void fem(std::vector<Node*> &nodes, std::vector<Element*> &elements, std::vector<Physical*> &physicals, std::vector<Parameter*> &parameters, std::map<Node*, std::vector<double> > &solution, int type, int conditions, std::vector<double> cond_periodiques);

#endif // FEM_H_INCLUDED
