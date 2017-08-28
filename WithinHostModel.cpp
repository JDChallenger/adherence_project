#include <iostream>
#include <cmath>
#include <stdio.h>
#include <fstream>
#include <random>
#include <cstdlib>
#include "parameters.h"
//#include <ctime>
#include <cstring>
#include <vector>
#include <chrono>

using namespace std;

int main () {

//Output file for results
const char output_File[] = "WHM_Results.txt";

//Vectors to store random growth rates
vector<double> R1(450,0.0); //These will be uncorrelated random numbers
vector<double> R ; //Random numbers correlated in time

//Construct a trivial random generator engine from a time-based seed:
unsigned seed = chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator (seed);

//Normal distribution, mean 0, standard deviation 1
normal_distribution<double> distribution (0,1);

//Uniform distribution, defined between 0 and 1
uniform_real_distribution<double> distribution2(0.0,1.0);


for(int j=0;j<450;j++){
	double number = distribution(generator);
	R1[j] = number;
}

/* Here we generate correlated random variables (used for the growth rates) from the uncorrelated random variables */
double sum = 0;
double rv = 0;
for(int j=0;j<450;j++){
	sum=0;
	for(int j1=1;j1<j+1;j1++){
		sum+=R1[j1]*pow(f,j-j1);
	}
	
	rv = pow(f,j)*R1[0]+sqrt(1-pow(f,2))*sum;

	//Check the random variable rv is within the required range- if so, add it to vector R
	if(mu + sigma * rv>1 && mu + sigma * rv<35){
  		R.push_back (mu + sigma * rv);
	}

}
//cout<<"Check length of R: "<<R.size()<<endl;//Are there enough for whole simulation?

//////////////////////////////////////////////////////////////////
/* Random numbers for the innate & general-adaptive immunities */
/////////////////////////////////////////////////////////////////
double Pms;
double Pcs;

double meanLN = 4.79*log(10);
double sigmaLN =1.2;

double rand1 = meanLN + sigmaLN * distribution(generator);
double rand2 = exp(rand1);

while(rand2 > pow(10,5.5)){
	rand1 = meanLN + sigmaLN * distribution(generator);
	rand2 = exp(rand1);
	//cout<<"Pcs Replaced"<<endl;
}

Pcs = kc * rand2;
cout<<"Pcs: " << Pcs <<endl;

//Draw from a Gompertz distribution
Pms = km * log(1-(1/gompcst2)*(log(1-distribution2(generator))))*(1/gompcst1);
cout<< "Pms: " << Pms <<endl;

/* Declare vector to store parasitaemia, and give initial condition at t=0 */
vector<double> P(400,0.0);
P[0]=0.1;//initial condition. Units are parasites per microlitre

vector<double> PC2(400,0.0);//Used to help calculate general-adaptive immune response
double PC;
double PCsum;
double Pv=0;

double Sc = 1;
double Sm = 1;
double Sv = 1;
int upper, lower;

for(int k=0;k<399;k++){ //Suggestion: start k loop from k=1, to match Mathematica, then readjust once it works?

	/* Innate immune response. Eq. 3 in article */
	Sc = 1/(1+pow((P[k]/Pcs),kaC));

	/* Effective var-specific adaptive immune response. Eq. 5 in article, with Suppl. Eq. 2a */ 
	//First determine the limits on the sum
	if(k>3){
		lower = round((k+1-4)*pow(lambda , k-4+1))-1;//Minus 1, since in Mathematica Initial Condition is at k=1
		
		upper = k - 4 +1;
		//cout<<"Lower: "<<lower << " Upper: "<<upper <<endl;

		Pv=0;
			for(int k2=lower;k2<upper;k2++){
				Pv += P[k2];
			}
		Sv=1/(1+pow((Pv/Pvs),kaV));
		//cout<<"k equals: "<<k<<" Lower: "<<lower << " Upper: "<<upper <<" Pv: "<<Pv<<endl;
			}

	/* General-adaptive immune response. Eq. 4 in article */
	if(P[k]>C){
		PC=C;
	}
	else{
		PC=P[k];
	}

	PC2[k]=PC;
		
	if(k>3){
	
		PCsum=0;
		for(int k1=0;k1<k-4+1;k1++){
			PCsum+=PC2[k1];
		}

		Sm = beta + (1-beta)/(1+pow((PCsum/Pms),kaM));
	}

	//Now combine all 3 immune responses with the growth rate to iterate model forward
	//This is Eq. 2 in the article
	P[k+1] = R[k] * Sc * Sm * Sv * P[k];
	//cout<<"On Day "<<k<<", " << P[k]<<" "<<Sc<<" "<<Sm<<" "<<Sv<<endl;
	if(P[k+1]<pow(10,-5)){
		cout <<"Episode ended!" <<endl; 
		break;
	}

}// end of Parasitaemia loop

cout<<"Now output results to file"<<endl; 
 	{
	  ofstream out(output_File);
	  if(!out){
	    cerr <<"Failed to open output file"<< output_File << endl;
	    exit(1);
	  }
	  for(int f3=0;f3<400;f3++){
	    out << 2*f3 <<"\t"<<P[f3]<<endl;
	  
	  }
	    out.close();
	}


  } //end of main
