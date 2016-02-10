#ifndef PHYSICALIO_H_INCLUDED
#define PHYSICALIO_H_INCLUDED

//Structure contenant les infos des paramètres du fichier .phy
struct Parameter
{
    std::string name;
    int dim;
    std::vector<double> value;//Ex : Temperature boundary or thermal diffusivity; -1 for null
};

void readPHY(const char *fileName, std::vector<Parameter*> &parameters, std::vector<double> &cond_periodiques);

#endif // PHYSICALIO_H_INCLUDED
