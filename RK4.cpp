
#include "parameters.h"

//////////////////////////////////////////////////
/* Define RK4 functions */
//////////////////////////////////////////////////
double f1AM(double x,double y, double z , params* theta)//(x,y,z)=(gut,central,metabolite)
{
double temp;
temp = -(theta->kaAM) * x;
return temp;
}

double f2AM(double x,double y, double z , params* theta)//(x,y,z)=(gut,central,metabolite)
{
double temp;
temp = (theta->kaAM) * x -((theta->k23AM) + ((theta->CLAM) / (theta->VCAM)))* y;
return temp;
}

double f3AM(double x,double y, double z , params* theta)//(x,y,z)=(gut,central,metabolite)
{
double temp;
temp = (theta->k23AM) * y - ((theta->CLmetAM) / (theta->VMAM))*z;
return temp;
}

double f1L(double x,double y, double z , params* theta)//(x,y,z)=(gut,central,metabolite)
{
double temp;
temp = -(theta->kaL) * x;
return temp;
}

double f2L(double x,double y, double z , params* theta)//(x,y,z)=(gut,central,metabolite)
{
double temp;
temp = (theta->kaL) * x -((theta->k23L) + ((theta->CLL) / (theta->VCL)))* y;
return temp;
}

double f3L(double x,double y, double z, params* theta)//(x,y,z)=(gut,central,metabolite)
{
double temp;
temp = (theta->k23L) * y - ((theta->CLmetL) / (theta->VML))*z;
return temp;
}


//And now for Gametocytes. Simple structure: in and out

double gf(double x, double y, double nux, double nuy) //In / out / in / out . Could put nu in a structure
{
double temp;
temp = nux * x - nuy * y;
return temp;
}

//Need more detailed structure for gametocytes with drugs? Different drug action (if any) for different stages
//Use struct for PD parameters? Would shorten function defn.
double gfD(double x, double y, double nux, double nuy, double AMconc, double DHAconc, double LMFconc/*, double PQconc*/, paramsPD* thetaPD)
{
double temp;
temp = nux * x - nuy * y - ((thetaPD->Gkmax_AMyoung * AMconc /(AMconc + thetaPD->c50_AM)) + 
	(thetaPD->Gkmax_AMyoung * DHAconc /(DHAconc + thetaPD->c50_AM)) + (thetaPD->Gkmax_Lyoung * LMFconc /(LMFconc + thetaPD->c50_L)))*y;
return temp;
}


