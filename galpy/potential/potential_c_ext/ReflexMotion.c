#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <galpy_potentials.h>
// Reflex Motion potential - based on Moving Object Potential
// 3 arguments: amp, t0, tf
void constrain_range2(double * d) {
  // Constrains index to be within interpolation range
  if (*d < 0) *d = 0.0;
  if (*d > 1) *d = 1.0;
}
double ReflexMotionRforce(double R,double z, double phi,
				   double t,
				   struct potentialArg * potentialArgs){
  double amp,t0,tf,d_ind,x,y,obj_x,obj_y,obj_z,obj_R,obj_phi,RF;
  double * args= potentialArgs->args;
  //Get args
  amp= *args;
  t0= *(args+1);
  tf= *(args+2);
  d_ind= (t-t0)/(tf-t0);
  x= R*cos(phi);
  y= R*sin(phi);
  constrain_range2(&d_ind);
  // Interpolate x, y, z
  obj_x= gsl_spline_eval(*potentialArgs->spline1d,d_ind,*potentialArgs->acc1d);
  obj_y= gsl_spline_eval(*(potentialArgs->spline1d+1),d_ind,
			 *(potentialArgs->acc1d+1));
  obj_z= gsl_spline_eval(*(potentialArgs->spline1d+2),d_ind,
			 *(potentialArgs->acc1d+2));
  // Rdist= pow(pow(x-obj_x, 2)+pow(y-obj_y, 2), 0.5);
  obj_R = pow(pow(obj_x,2) + pow(obj_y,2),0.5);
  obj_phi = atan2(obj_y,obj_x);
  // Calculate R force
  RF= calcRforce(obj_R,obj_z,obj_phi,t,potentialArgs->nwrapped,
		 potentialArgs->wrappedPotentialArg);
  return -amp*RF*(cos(phi)*cos(obj_phi)+sin(phi)*sin(obj_phi));
}

double ReflexMotionzforce(double R,double z,double phi,
				      double t,
				      struct potentialArg * potentialArgs){
  double amp,t0,tf,d_ind,x,y,obj_x,obj_y,obj_z, obj_R, obj_phi;
  double * args= potentialArgs->args;
  //Get args
  amp= *args;
  t0= *(args+1);
  tf= *(args+2);
  d_ind= (t-t0)/(tf-t0);
  x= R*cos(phi);
  y= R*sin(phi);
  constrain_range2(&d_ind);
  // Interpolate x, y, z
  obj_x= gsl_spline_eval(*potentialArgs->spline1d,d_ind,*potentialArgs->acc1d);
  obj_y= gsl_spline_eval(*(potentialArgs->spline1d+1),d_ind,
			 *(potentialArgs->acc1d+1));
  obj_z= gsl_spline_eval(*(potentialArgs->spline1d+2),d_ind,
			 *(potentialArgs->acc1d+2));
  // Rdist= pow(pow(x-obj_x, 2)+pow(y-obj_y, 2), 0.5);
  obj_R = pow(pow(obj_x,2) + pow(obj_y,2),0.5);
  obj_phi = atan2(obj_y,obj_x);
  
  // Calculate z force
  return -amp * calczforce(obj_R,obj_z,obj_phi,t,potentialArgs->nwrapped,
			   potentialArgs->wrappedPotentialArg);
}

double ReflexMotionphiforce(double R,double z,double phi,
					double t,
					struct potentialArg * potentialArgs){
  double amp,t0,tf,d_ind,x,y,obj_x,obj_y,obj_z, obj_R, obj_phi, RF;
  double * args= potentialArgs->args;
  //Get args
  amp= *args;
  t0= *(args+1);
  tf= *(args+2);
  d_ind= (t-t0)/(tf-t0);
  x= R*cos(phi);
  y= R*sin(phi);
  constrain_range2(&d_ind);
  // Interpolate x, y, z
  obj_x= gsl_spline_eval(*potentialArgs->spline1d,d_ind,*potentialArgs->acc1d);
  obj_y= gsl_spline_eval(*(potentialArgs->spline1d+1),d_ind,
			 *(potentialArgs->acc1d+1));
  obj_z= gsl_spline_eval(*(potentialArgs->spline1d+2),d_ind,
			 *(potentialArgs->acc1d+2));
  // Rdist= pow(pow(x-obj_x, 2)+pow(y-obj_y, 2), 0.5);
  obj_R = pow(pow(obj_x,2) + pow(obj_y,2),0.5);
  obj_phi = atan2(obj_y,obj_x);

  // Calculate phiforce
  RF= calcRforce(obj_R,obj_z,obj_phi,t,potentialArgs->nwrapped,
		 potentialArgs->wrappedPotentialArg);
  return amp*RF*(sin(phi)*cos(obj_phi)-cos(phi)*sin(obj_phi));
}

// double MovingObjectPotentialPlanarRforce(double R, double phi,
// 				      double t,
// 				      struct potentialArg * potentialArgs){
//  double amp,t0,tf,d_ind,x,y,obj_x,obj_y,Rdist,RF;
//  double * args= potentialArgs->args;
  //Get args
//  amp= *args;
//  t0= *(args+1);
//  tf= *(args+2);
//  d_ind= (t-t0)/(tf-t0);
//  x= R*cos(phi);
//  y= R*sin(phi);
//  constrain_range(&d_ind);
  // Interpolate x, y
//  obj_x= gsl_spline_eval(*potentialArgs->spline1d,d_ind,*potentialArgs->acc1d);
//  obj_y= gsl_spline_eval(*(potentialArgs->spline1d+1),d_ind,
//			 *(potentialArgs->acc1d+1));
//  Rdist= pow(pow(x-obj_x, 2)+pow(y-obj_y, 2), 0.5);
  // Calculate R force
//  RF= calcPlanarRforce(Rdist, phi, t, potentialArgs->nwrapped,
//		       potentialArgs->wrappedPotentialArg);
//  return -amp*RF*(cos(phi)*(obj_x-x)+sin(phi)*(obj_y-y))/Rdist;
//}

// double MovingObjectPotentialPlanarphiforce(double R, double phi,
//					double t,
//					struct potentialArg * potentialArgs){
//  double amp,t0,tf,d_ind,x,y,obj_x,obj_y,Rdist,RF;
//  double * args= potentialArgs->args;
  // Get args
//  amp= *args;
//  t0= *(args+1);
//  tf= *(args+2);
//  d_ind= (t-t0)/(tf-t0);
//  x= R*cos(phi);
//  y= R*sin(phi);
//  constrain_range(&d_ind);
  // Interpolate x, y
//  obj_x= gsl_spline_eval(*potentialArgs->spline1d,d_ind,*potentialArgs->acc1d);
//  obj_y= gsl_spline_eval(*(potentialArgs->spline1d+1),d_ind,
//			 *(potentialArgs->acc1d+1));
//  Rdist= pow(pow(x-obj_x, 2)+pow(y-obj_y, 2), 0.5);
  // Calculate phiforce
//  RF= calcPlanarRforce(Rdist, phi, t, potentialArgs->nwrapped,
//		       potentialArgs->wrappedPotentialArg);
//  return -amp*RF*R*(cos(phi)*(obj_y-y)-sin(phi)*(obj_x-x))/Rdist;
//}
