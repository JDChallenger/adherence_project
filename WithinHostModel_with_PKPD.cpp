#include <iostream>
#include <cmath>
#include <stdio.h>
#include <fstream>
#include <random>
#include <cstdlib>
#include "parameters.h"
#include <cstring>
#include <vector>
#include <chrono>

using namespace std;

/*Function & structure declarations in the header file parameters.h*/

int main () {

//Output file for results
const char output_File1[] = "WHM_Results_with_drugs.txt";
const char output_File2[] = "WHM_Results_with_PKPD.txt";


//Vectors to store random growth rates
vector<double> R1(450,0.0); //These will be uncorrelated random numbers
vector<double> R ; //Random numbers correlated in time

// construct a trivial random generator engine from a time-based seed:
unsigned seed = chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator (seed);

//Normal distribution, mean 0, standard deviation 1
normal_distribution<double> distribution (0,1);

//Uniform distribution, defined between 0 and 1
uniform_real_distribution<double> distribution2(0.0,1.0);

//Uniform distribution, used for the pyrogenic threshold
uniform_real_distribution<double> distribution3(-3.69,-0.0);

for(int j=0;j<450;j++){
double number = distribution(generator);
R1[j] = number;
}

/* Here we generate correlated random variables (used for the growth rates) from the uncorrelated random variables */
double sum=0;
double rv=0;
for(int j=0;j<450;j++){
	sum=0;
	for(int j1=1;j1<j+1;j1++){
		sum+=R1[j1]*pow(f,j-j1);
	}
rv = pow(f,j)*R1[0]+sqrt(1-pow(f,2))*sum;
if(mu + sigma * rv>1 && mu + sigma * rv <35){
  R.push_back (mu + sigma * rv);
}

}

//////////////////////////////////////////////////////////////////
/* Random numbers for the innate & general-adaptive immunities */
/////////////////////////////////////////////////////////////////

double Pms;
double Pcs;

double meanLN = 4.79*log(10);
double sigmaLN =1.2;

double rand1 = meanLN + sigmaLN * distribution(generator);
double rand2 = exp(rand1);
//cout<< rand2 <<" upper truncation: "<<pow(10,5.5)<<endl;

while(rand2 > pow(10,5.5)){
rand1 = meanLN + sigmaLN * distribution(generator);
rand2 = exp(rand1);
}

Pcs = kc * rand2;
cout<<"Pcs: " << Pcs <<endl;

//Draw from a Gompertz distribution
Pms = km * log(1-(1/gompcst2)*(log(1-distribution2(generator))))*(1/gompcst1);
cout<< "Pms: " << Pms <<endl;

//////////////////////////////////////////////////
/* Variables for drugs, and ODE time step */
//////////////////////////////////////////////////

//Choose a value of dt (often called h in RK4)
double dt = 0.025;
int TL = (1/dt) * 24 * 60;// For 60 days' worth
cout <<"For 60 days we require "<< TL <<" increments" <<endl;

//Also useful to have the number of intervals for one & two days (one generation)
int oneday = (1/dt) * 24;
int twoday = (1/dt) * 48;

//Timings for the six doses. Convert from hours into increments.
int D1 = 0 * (1 / dt);
int D2 = 8 * (1 / dt);
int D3 = 24 * (1 / dt);
int D4 = 36 * (1 / dt);
int D5 = 48 * (1 / dt);
int D6 = 60 * (1 / dt);

/*Variable for timing of 1st dose, DELTA, which will depend on fever timing + delay to 
access treatment. This initialisation means treatment won't happen over timespan
 of the model. However, value will be updated after pyrogenic threshold exceded */
int DELTA = TL + 5;

double DRUGL;
double DRUGAM;

//Arrays for Artmenter model
double GUTAM[TL];
double CENTRALAM[TL];
double METABOLITEAM[TL];
GUTAM[0] = 0.0;
CENTRALAM[0] = 0.0;
METABOLITEAM[0] = 0.0;

//Arrays for Lumefantrine variable
double GUTL[TL];
double CENTRALL[TL];
double METABOLITEL[TL];
GUTL[0] = 0.0;
CENTRALL[0] = 0.0;
METABOLITEL[0] = 0.0;

//Then convert to the right units
double CENTRALAMx[TL];
double METABOLITEAMx[TL];
CENTRALAMx[0] = 0.0;
METABOLITEAMx[0] = 0.0;

double CENTRALLx[TL];
double METABOLITELx[TL];
CENTRALLx[0] = 0.0;
METABOLITELx[0] = 0.0;

///////////////////////////////////////////////////////////////////////
/* Now choose body weight and, therefore, number of tablets per dose */
///////////////////////////////////////////////////////////////////////
params theta;
//Body Weight (kg)
double BW = 11.1 + 2.8 * distribution(generator);//SHOULD BE RANDOM, AND DOSE SHOULD BE ADJUSTED ACCORDINGLY
while(BW<5||BW>25){
	BW = 11.1 + 2.8 * distribution(generator);
}

int TABLETS; //How many tablets per dose?
if(5 < BW && BW < 15){
	TABLETS = 1;
}
if(15 < BW && BW < 25){
	TABLETS = 2;
}
if(25 < BW && BW < 35){
	TABLETS = 3;
}
if(BW > 35){
	TABLETS = 4;
}
cout<<"Body weight: "<< BW <<". How many tablets / dose ? " << TABLETS <<endl;

//Now introduce the inter-individual variation for the PK model 
double r1, r2, r3, r4, r5, r6;

r1 = 0.44 * distribution(generator); // CLAM
r2 = 1.19 * distribution(generator); // kaAM
r3 = 0.68 * distribution(generator); // k23AM
r4 = 0.38 * distribution(generator); // CLL
r5 = 0.6 * distribution(generator); // VCL
r6 = 0.38 * distribution(generator); // k23L

//Artemether constants
theta.CLAM = 24.7 * pow(BW , 0.75) * exp(r1);
theta.VCAM = 129;
theta.VMAM = theta.VCAM;
theta.kaAM = 0.27 * exp(r2);
theta.k23AM = 5.86 * exp(r3);
theta.CLmetAM = 419;

//Lumefantrine constants
theta.CLL = 0.84 * pow(BW , 0.52) * exp(r4);
theta.VCL = 59.9 * pow(BW , 0.35) * exp(r5);
theta.VML = theta.VCL;
theta.CLmetL = 4.8;
theta.kaL = 0.54;
theta.k23L = 3.7 * pow(10,-4) * exp(r6);

//////////////////////////////////////////////////
/* Define PD parameters in the struct. 		*/
//////////////////////////////////////////////////

paramsPD thetaPD;

thetaPD.kmax_AM = 0.189;
thetaPD.c50_AM = 3.3;

//fix DHA values to AM values

thetaPD.kmax_L = 0.165;
thetaPD.c50_L = 125.0;

//NOTE: killing effect of metabolite DLF turned off at present
thetaPD.kmax_DLF = 0.0;
thetaPD.c50_DLF = 0.0;

////////////////////////////////////////////////////////////////////////
/* Dummy run. Max. parasite density of an untreated infection (with same 
random variables) will be used to determine pyrogenic threshold 	  */
////////////////////////////////////////////////////////////////////////

vector<double> Q(30,0.0);//vector to help assess Day0 & Day1 parasitaemia
Q[0]=0.1;//initial condition

vector<double> PC2q(30,0.0);
double PCq;
double PCsumq;
double Pvq=0;

double Scq = 1;
double Smq = 1;
double Svq = 1;
int upperq, lowerq;
double Qmax=0;//Peak parasitaemia for untreated episode
for(int k=0;k<29;k++){ 

//Innate immune response
	Scq = 1/(1+pow((Q[k]/Pcs),kaC));

//Effective var-specific adaptive immune response
	if(k>3){
		lowerq = round((k+1-4)*pow(lambda , k-4+1))-1;
		upperq = k - 4 +1;

	Pvq=0;
		for(int k2=lowerq;k2<upperq;k2++){
Pvq += Q[k2];
		}
Svq=1/(1+pow((Pvq/Pvs),kaV));
}
// General-adaptive immune response
if(Q[k]>C){
		PCq=C;
	}
	else{
		PCq=Q[k];
	}

	PC2q[k]=PCq;
		
	if(k>3){
	
	PCsumq=0;
	for(int k1=0;k1<k-4+1;k1++){
		PCsumq+=PC2q[k1];
	}
	Smq = beta + (1-beta)/(1+pow((PCsumq/Pms),kaM));
}

Q[k+1] = R[k] * Scq * Smq * Svq * Q[k];
if(Q[k+1] > Qmax){
	Qmax = Q[k+1];
}
if(Q[k+1]<pow(10,-5)){
	//cout <<"Episode ended!" <<endl;
	break;
}
}// end of Parasitaemia loop for 'dummy run'

//cout<<"Qmax is: "<< Qmax <<endl;

//////////////////////////////////////////////////
/* 			   End of dummy run 				*/
//////////////////////////////////////////////////

//////////////////////////////////////////////////
/* Fever Threshold & Waiting time for treatment */
//////////////////////////////////////////////////

//Determine pyrogenic threshold
double Ur = distribution3(generator);
double thresh = pow(10,Ur + log10(Qmax) );

cout<<"Fever Threshold:"<< thresh <<endl;

//Delay from fever breaking to 1st dose administered
double meanLNw = 0.7;
double sigmaLNw = 0.6575;

double rand1w = meanLNw + sigmaLNw * distribution(generator);
double rand2w = exp(rand1w);

while(rand2w>10.0){//truncate dist
	rand1w = meanLNw + sigmaLNw * distribution(generator);
	rand2w = exp(rand1w);
}
int rand3w = round(rand2w * 24 *(1/double(dt) ));
cout<<"Waiting time (in days): "<< rand2w <<", and in increments: "<<rand3w<<endl;

/////////////////////////////////////////////////////////////////////////////
/* Declare vector to store parasitaemia, and give initial condition at t=0 */
/////////////////////////////////////////////////////////////////////////////
vector<double> P(400,0.0);
P[0]=0.1;//initial condition

vector<double> PC2(400,0.0);
double PC;
double PCsum;
double Pv=0;

int AS = 0;
int CF = 0;
int fevk = 0;

double Sc = 1;
double Sm = 1;
double Sv = 1;
int upper, lower;

//////////////////////////////////////////////////
/* Begin Runge-Kutta routine */
//////////////////////////////////////////////////

double k1L[3], k2L[3], k3L[3], k4L[3];
double xL, yL, zL;
double k1AM[3], k2AM[3], k3AM[3], k4AM[3];
double xAM, yAM, zAM;
double killF = 0.0;

for(int k = 0; k < TL; k++){
	if(k==D1 + DELTA||k==D2 + DELTA||k==D3 + DELTA||k==D4 + DELTA||k==D5 + DELTA||k==D6 + DELTA)
	{
		DRUGL = TABLETS * 120;
		DRUGAM = TABLETS * 20;
	}
	else
	{
		DRUGL = 0.0;
		DRUGAM = 0.0;
	}

	GUTAM[k] += DRUGAM;	
	xAM = GUTAM[k];
	yAM = CENTRALAM[k];
	zAM = METABOLITEAM[k];

	k1AM[0] = dt * f1AM(xAM, yAM, zAM, &theta);
	k1AM[1] = dt * f2AM(xAM, yAM, zAM, &theta);
	k1AM[2] = dt * f3AM(xAM, yAM, zAM, &theta);

	k2AM[0] = dt * f1AM(xAM + 0.5 * k1AM[0], yAM + 0.5 * k1AM[1], zAM + 0.5 * k1AM[2], &theta);
	k2AM[1] = dt * f2AM(xAM + 0.5 * k1AM[0], yAM + 0.5 * k1AM[1], zAM + 0.5 * k1AM[2], &theta);
	k2AM[2] = dt * f3AM(xAM + 0.5 * k1AM[0], yAM + 0.5 * k1AM[1], zAM + 0.5 * k1AM[2], &theta);

	k3AM[0] = dt * f1AM(xAM + 0.5 * k2AM[0], yAM + 0.5 * k2AM[1], zAM + 0.5 * k2AM[2], &theta);
	k3AM[1] = dt * f2AM(xAM + 0.5 * k2AM[0], yAM + 0.5 * k2AM[1], zAM + 0.5 * k2AM[2], &theta);
	k3AM[2] = dt * f3AM(xAM + 0.5 * k2AM[0], yAM + 0.5 * k2AM[1], zAM + 0.5 * k2AM[2], &theta);

	k4AM[0] = dt * f1AM(xAM + k3AM[0], yAM + k3AM[1], zAM + k3AM[2], &theta);
	k4AM[1] = dt * f2AM(xAM + k3AM[0], yAM + k3AM[1], zAM + k3AM[2], &theta);
	k4AM[2] = dt * f3AM(xAM + k3AM[0], yAM + k3AM[1], zAM + k3AM[2], &theta);

	GUTAM[k+1] = GUTAM[k] +  (k1AM[0] + 2 * k2AM[0] + 2 * k3AM[0] + k4AM[0])/6;
	CENTRALAM[k+1] = CENTRALAM[k] + (k1AM[1] + 2 * k2AM[1] + 2 * k3AM[1] + k4AM[1])/6;
	METABOLITEAM[k+1] = METABOLITEAM[k] + (k1AM[2] + 2 * k2AM[2] + 2 * k3AM[2] + k4AM[2])/6;

	//Convert to the desired units!
	CENTRALAMx[k+1] = 1000 * (CENTRALAM[k+1]/theta.VCAM) ;
	METABOLITEAMx[k+1] = 1000 * (METABOLITEAM[k+1]/theta.VMAM) ;


	GUTL[k]+= DRUGL;	
	xL = GUTL[k];
	yL = CENTRALL[k];
	zL = METABOLITEL[k];

	k1L[0] = dt * f1L(xL, yL, zL, &theta);
	k1L[1] = dt * f2L(xL, yL, zL, &theta);
	k1L[2] = dt * f3L(xL, yL, zL, &theta);

	k2L[0] = dt * f1L(xL + 0.5 * k1L[0], yL + 0.5 * k1L[1], zL + 0.5 * k1L[2], &theta);
	k2L[1] = dt * f2L(xL + 0.5 * k1L[0], yL + 0.5 * k1L[1], zL + 0.5 * k1L[2], &theta);
	k2L[2] = dt * f3L(xL + 0.5 * k1L[0], yL + 0.5 * k1L[1], zL + 0.5 * k1L[2], &theta);

	k3L[0] = dt * f1L(xL + 0.5 * k2L[0], yL + 0.5 * k2L[1], zL + 0.5 * k2L[2], &theta);
	k3L[1] = dt * f2L(xL + 0.5 * k2L[0], yL + 0.5 * k2L[1], zL + 0.5 * k2L[2], &theta);
	k3L[2] = dt * f3L(xL + 0.5 * k2L[0], yL + 0.5 * k2L[1], zL + 0.5 * k2L[2], &theta);

	k4L[0] = dt * f1L(xL + k3L[0], yL + k3L[1], zL + k3L[2], &theta);
	k4L[1] = dt * f2L(xL + k3L[0], yL + k3L[1], zL + k3L[2], &theta);
	k4L[2] = dt * f3L(xL + k3L[0], yL + k3L[1], zL + k3L[2], &theta);

	GUTL[k+1] = GUTL[k] +  (k1L[0] + 2 * k2L[0] + 2 * k3L[0] + k4L[0])/6;
	CENTRALL[k+1] = CENTRALL[k] + (k1L[1] + 2 * k2L[1] + 2 * k3L[1] + k4L[1])/6;
	METABOLITEL[k+1] = METABOLITEL[k] + (k1L[2] + 2 * k2L[2] + 2 * k3L[2] + k4L[2])/6;

	//Convert to the desired units!
	CENTRALLx[k+1] = 1000 * (CENTRALL[k+1]/theta.VCL);
	METABOLITELx[k+1] = 1000 * (METABOLITEL[k+1]/theta.VML);

	//update 'kill factor'
	killF -= dt * ( thetaPD.kmax_L * (CENTRALLx[k]/(CENTRALLx[k]+thetaPD.c50_L)) /*+
	thetaPD.kmax_DLF * (METABOLITELx[k]/(METABOLITELx[k]+thetaPD.c50_DLF))*/ +
	thetaPD.kmax_AM * (CENTRALAMx[k]/(CENTRALAMx[k]+thetaPD.c50_AM)) + 
	thetaPD.kmax_AM * (METABOLITEAMx[k]/(METABOLITEAMx[k]+thetaPD.c50_AM)) );

	//Do we need to update the parasitaemia on this time step (every 48hrs)?
	if(k%(twoday)==twoday-1 && AS==0){ 
		int k1 = ((k+1) * dt /48 )-1;
		//cout<<"Check k1: "<<k1<<endl;

		//Innate immune response
		Sc = 1/(1+pow((P[k1]/Pcs),kaC));

		//Effective var-specific adaptive immune response
		if(k1>3){
			lower = round((k1+1-4)*pow(lambda , k1-4+1))-1;//Minus 1, since in Mathematica Initial Condition is at k1=1
			upper = k1 - 4 + 1;
			Pv=0;
			for(int k2=lower;k2<upper;k2++){
				Pv += P[k2];
			}
			Sv=1/(1+pow((Pv/Pvs),kaV));
		}
		// General-adaptive immune response
		if(P[k1]>C){
			PC=C;
		}
		else{
			PC=P[k1];
		}

		PC2[k1]=PC;
		
		if(k1>3){
			PCsum=0;
			for(int k2=0;k2<k1-4+1;k2++){
				PCsum+=PC2[k2];
			}
		Sm = beta + (1-beta)/(1+pow((PCsum/Pms),kaM));
		}

		cout<<"Kill factor: "<< killF <<" ,and with exp: " <<exp(killF)<<endl;

		P[k1+1] = R[k1] * Sc * Sm * Sv * P[k1] * exp(killF);

		cout<<"On Day "<< 2 * k1 <<", " << P[k1] <<" "<< Sc <<" "<< Sm <<" "<< Sv <<endl;

		if(P[k1+1]>thresh && CF==0){
			cout<<"Fever triggered!"<<endl;
			DELTA = k + rand3w;
			fevk = k;
			CF = 1;
		}

		if(P[k1+1]<pow(10,-5)){
			AS=1;
		}

		killF = 0.0; //reset
	}// end of Parasitaemia loop

	if(AS==1){
		cout<<"Episode over!"<<endl;
	break;
	}
	//finish simulation if 'Delta+28days exceeded?'
	/*if(k > (DELTA+(30 * oneday) ) ){
		cout<<"Delta + 28 days reached"<<endl;
	break;
	}*/
}// end of RK4 loop for PK ODEs


//////////////////////////////////////////////
/*					Output 					*/
//////////////////////////////////////////////

cout<<"Now output results to file"<<endl; 
 	{
	  ofstream out1(output_File1);
	  if(!out1){
	    cerr <<"Failed to open output file"<< output_File1 << endl;
	    exit(1);
	  }
	  for(int f3=0;f3<400;f3++){
	    out1 << 2*f3 <<"\t"<<P[f3]<<"\t"<<Q[f3]<<endl;
	  }
	    out1.close();
	}

	{
	  ofstream out2(output_File2);
	  if(!out2){
	    cerr <<"Failed to open output file"<< output_File2 << endl;
	    exit(1);
	  }

	  for(int f3 = 0;f3 < TL;f3+=5){
	    out2 << (1/24.0) * f3 * dt <<"\t"<<CENTRALAMx[f3]<<"\t"<<METABOLITEAMx[f3]<<"\t"<<CENTRALLx[f3]<<"\t"<<METABOLITELx[f3]<<endl;
	  
	  }
	  out2<< BW <<"\t"<< fevk <<"\t"<< dt <<"\t"<< thresh <<"\t"<< DELTA <<endl;
	    out2.close();
	}

  } //end of main
