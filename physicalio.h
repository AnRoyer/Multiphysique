#ifndef PHYSICALIO_H_INCLUDED
#define PHYSICALIO_H_INCLUDED

//Structure contenant les infos des paramètres du fichier .phy
struct Parameter
{
    std::string name;
    int dim;
    //line parameters
    double temperature;
    double voltage;
    //surface parameters
    std::vector<double> thermalConductivity;
    std::vector<double> electricalConductivity;
    double thermalGeneration;
    double electricalGeneration;
};

struct Periodique
{
    double meanTemperature;
    double xGradient;
    double yGradient;
    bool exist;
};

struct XMLparam
{
    std::string name;
    std::string value;
};

enum Type
{
    DIRICHLET,
    PERIODIC
};

enum Nature
{
    THERMAL,
    ELECTRICAL
};

enum Dim
{
    LINE,
    SURFACE,
    GLOBAL
};

void readPHY(const char *fileName, std::vector<Parameter*> &parameters, Periodique &conditions);
XMLparam readParam(std::ifstream& fp);
double readValue(std::ifstream& fp);

#endif // PHYSICALIO_H_INCLUDED
