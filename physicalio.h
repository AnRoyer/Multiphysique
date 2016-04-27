#ifndef PHYSICALIO_H_INCLUDED
#define PHYSICALIO_H_INCLUDED

enum Type
{
    DIRICHLET,
    PERIODIC,
    VONNEUMANN,
    FE2withDIRICHLET,
    FE2withVONNEUMANN,
    FE2withPERIODIC,
    DEFAULTTYPE
};

enum Nature
{
    THERMAL,
    ELECTRICAL,
    DEFAULTNATURE
};

enum Dim
{
    LINE,
    SURFACE,
    GLOBAL,
    DEFAULTDIM
};

//Structure contenant les infos des paramètres du fichier .phy
struct Conductivity
{
        std::string name;
        double conductivity[2][2];
};

struct Parameter
{
	Parameter()
	{
		dim = -1;
		temperature = -1;
		voltage = -1;
		thermalGeneration = -1;
		electricalGeneration = -1;
		fluxTemperature = -1;
		fluxVoltage = -1;
	}

    std::string name;
    int dim;
    //line parameters
    double temperature;
    double voltage;

    double fluxTemperature;
    double fluxVoltage;// **************** A VERIFIER
    //surface parameters
    std::vector<Conductivity*> thermalConductivity;
    std::vector<Conductivity*> electricalConductivity;
    double thermalGeneration;
    double electricalGeneration;
};

struct Periodique
{
	Periodique()
	{
		meanTemperature = -1;
		xGradient = -1;
		yGradient = -1;
		exist = false;
	}
    double meanTemperature;
    double xGradient;
    double yGradient;
    bool exist;
};

struct Micro
{
    std::string fileMsh;
    std::string filePhy;
};

struct XMLparam
{
    std::string name;
    std::string value;
};

void readPHY(const char *fileName, std::vector<Parameter*> &parameters, Periodique &conditions, Micro &micro, Type &typeUsed, double &eps);
XMLparam readParam(std::ifstream& fp);
double readValue(std::ifstream& fp);

#endif // PHYSICALIO_H_INCLUDED
