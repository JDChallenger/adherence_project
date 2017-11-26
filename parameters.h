// Header file

//Guard
#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <vector>

using namespace std;

//struct for PK parameters?
struct params
{
	//Body Weight (kg)
double BW;

//Artemether constants
double CLAM;//Write as function of BW?
double VCAM;
double VMAM;
double kaAM;
double k23AM;
double CLmetAM;

//Lumefantrine constants
double CLL;//Write as a function?
double VCL;
double VML;
double CLmetL;
double kaL;
double k23L;

};

//separate struct for PD parameters?
struct paramsPD
{
double kmax_AM;
double c50_AM;

//double kmax_DHA; fix DHA params to those for AM???
//double c50_DHA;

double kmax_L;
double c50_L;

double kmax_DLF;
double c50_DLF;

double Gkmax_AMyoung;
double Gkmax_AMold;

//fix DHA values to AM values

double Gkmax_Lyoung;
double Gkmax_Lold;
};

// User defined functions, for the Runge Kutta steps
double f1AM(double x, double y, double z, params* theta);//(x,y,z)=(gut,central,metabolite)
double f2AM(double x, double y, double z, params* theta);
double f3AM(double x, double y, double z, params* theta);

double f1L(double x, double y, double z, params* theta);
double f2L(double x, double y, double z, params* theta);
double f3L(double x, double y, double z, params* theta);

//Implementation of the analytical solution for the Artemether PK equation
vector<double> ART(double Time, double dose, params* theta);

//////////////////////////////////////////////////
/* Asexual parasite parameters here */
//////////////////////////////////////////////////

const double beta = 0.01; //Max. strength of general adaptive response
const double C = 1; //Parasitaemia threshold for general-adpative response to strength

const int kaC = 3;//Hill slope of innate immune response
const int kaV = 3;//Hill slope of EVS immune response
const int kaM= 1 ;//Hill slope of general adaptive immune response

const double Pvs = 10100;
const double lambda = 0.9996; //Controls duration of memory for EVS immune response

const double km = 0.021;
const double kc = 0.164;

const double gompcst1 = 0.0311;
const double gompcst2 = 0.0004;

//Constants for stochastic parasite growth rates
const double f = 0.64; // Autocorrelation strength
const double mu = 16.0; // Mean growth rate
const double sigma = 8.7;//8.0; // Standard Deviation of growth rate

//End of header guard

#endif