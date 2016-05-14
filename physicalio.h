#ifndef PHYSICALIO_H_INCLUDED
#define PHYSICALIO_H_INCLUDED

#define ELECTRICALDATA 0x0001
#define THERMALDATA 0x0002

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
        meanVoltage = -1;
		xGradient = -1;
		yGradient = -1;
		exist = false;
	}
    double meanTemperature;
    double meanVoltage;
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

void readPHY(const char *fileName, std::vector<Parameter*> &parameters, Periodique &conditions, Micro &micro,
             Type &typeUsed_thermic, Type &typeUsed_electric, int &nature, double &eps, std::vector<int> &methodFE2, int &natureFlag);
XMLparam readParam(std::ifstream& fp);
double readValue(std::ifstream& fp);

#endif // PHYSICALIO_H_INCLUDED
