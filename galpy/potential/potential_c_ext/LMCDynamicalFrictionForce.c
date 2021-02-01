#include <math.h>
#include <galpy_potentials.h>
// LMCDynamicalFrictionForce: 3 arguments: amp,ms,lnLambda 

double LMCDynamicalFrictionForceAmplitude(double R,double z, 
						    double phi,double t,
						    double r2,
						    struct potentialArg * potentialArgs,
						    double vR,double vT,
						    double vz){
  double forceAmplitude;
  double * args= potentialArgs->args;
  //Get args
  double amp= *args;
  double ms= *(args+1);
  double lnLambda= *(args+2);
  forceAmplitude= - amp*0.428*ms*lnLambda/r2;
  return forceAmplitude;
}
double LMCDynamicalFrictionForceRforce(double R,double z, double phi,
						 double t,
						 struct potentialArg * potentialArgs,
						 double vR,double vT,
						 double vz){
  double forceAmplitude;
  double r2=  R * R + z * z;
  forceAmplitude= LMCDynamicalFrictionForceAmplitude(R,z,phi,t,r2,
								 potentialArgs,
 								 vR,vT,vz);
   return forceAmplitude * vR;
}
double LMCDynamicalFrictionForcezforce(double R,double z, double phi,
						 double t,
						 struct potentialArg * potentialArgs,
						 double vR,double vT,
						 double vz){
  double forceAmplitude;
  double r2=  R * R + z * z;
  forceAmplitude= LMCDynamicalFrictionForceAmplitude(R,z,phi,t,r2,
								 potentialArgs,
								 vR,vT,vz);
  return forceAmplitude * vz;
}
double LMCDynamicalFrictionForcephiforce(double R,double z,
						   double phi,double t,
						   struct potentialArg * potentialArgs,
						   double vR,double vT,
						   double vz){
  double forceAmplitude;
  double r2=  R * R + z * z;
  forceAmplitude= LMCDynamicalFrictionForceAmplitude(R,z,phi,t,r2,
								 potentialArgs,
								 vR,vT,vz);
  return forceAmplitude * vT * R;
}

