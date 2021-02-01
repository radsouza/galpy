#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <galpy_potentials.h>
// MovingObjectPotential
// 3 arguments: amp, t0, tf
void constrain_range3(double * d) {
  // Constrains index to be within interpolation range
  if (*d < 0) *d = 0.0;
  if (*d > 1) *d = 1.0;
}

void cosinerule(double R1, double R2, double Rdist, double * cosphi, double *sinphi) {
// R1 is the distane of the LMC from the MW
// R2 is the distance of the Sat from the MW
// Rdist is the estimated distance of the sat from the LMC
  *cosphi = (R1*R1 - R2*R2 - Rdist*Rdist)/(2*R2*Rdist);
  *sinphi = pow(1 - (*cosphi)*(*cosphi), 0.5);
}


double MovingObjectDissipativeRforce(double R,double z, double phi,
				   double t,
				   struct potentialArg * potentialArgs,
				   double vR,double vT,
				   double vz){
  double amp,t0,tf,d_ind, Rdist, Rorb, RF, TF;
  double x,y,vx,vy, obj_x,obj_y,obj_z,obj_vx,obj_vy,obj_vz;
  double costheta, sintheta, phi2;
  double vxd, vyd, vzd, vR1, vT1;
  double * args= potentialArgs->args;
  //Get args
  amp= *args;
  t0= *(args+1);
  tf= *(args+2);
  d_ind= (t-t0)/(tf-t0);
  constrain_range3(&d_ind);

  // convert to cartesian coodinates
  x= R*cos(phi);
  y= R*sin(phi);
  vx=  vR*cos(phi) - vT*sin(phi);
  vy=  vR*sin(phi) + vT*cos(phi);

  // Interpolate x, y, z, vx, vy, vz
  obj_x= gsl_spline_eval(*potentialArgs->spline1d,d_ind,*potentialArgs->acc1d);
  obj_y= gsl_spline_eval(*(potentialArgs->spline1d+1),d_ind,
			 *(potentialArgs->acc1d+1));
  obj_z= gsl_spline_eval(*(potentialArgs->spline1d+2),d_ind,
			 *(potentialArgs->acc1d+2));
  obj_vx= gsl_spline_eval(*(potentialArgs->spline1d+3),d_ind,
			 *(potentialArgs->acc1d+3));
  obj_vy= gsl_spline_eval(*(potentialArgs->spline1d+4),d_ind,
			 *(potentialArgs->acc1d+4));
  obj_vz= gsl_spline_eval(*(potentialArgs->spline1d+5),d_ind,
			 *(potentialArgs->acc1d+5));

  // convert velocity to reference frame of LMC
  vxd = vx - obj_vx;
  vyd = vy - obj_vy;
  vzd = vz - obj_vz;

  // determine angles for transformation
  Rdist= pow(pow(x-obj_x, 2)+pow(y-obj_y, 2), 0.5);
  Rorb = pow(obj_x*obj_x + obj_y*obj_y, 0.5);
  phi2 = atan2(y-obj_y, x-obj_x);           
  cosinerule(Rorb,R,Rdist,&costheta,&sintheta);

  // convert velocity to cyclindrical coordinates in reference frame of LMC.
  vR1=    cos(phi2)*vxd    + sin(phi2)*vyd;        
  vT1=   -sin(phi2)*vxd    + cos(phi2)*vyd;

  // Calculate Radial and transverse forces
  RF= calcRforce(Rdist,(z-obj_z),phi,t,potentialArgs->nwrapped,
		 potentialArgs->wrappedPotentialArg,vR1,vT1,vzd);
  TF= calcPhiforce(Rdist,(z-obj_z),phi,t,potentialArgs->nwrapped,
		 potentialArgs->wrappedPotentialArg,vR1,vT1,vzd);

  // Transform force back into the refernce frame of the MW
  return  -amp*(RF*costheta  +TF*sintheta);
}

double MovingObjectDissipativezforce(double R,double z,double phi,
				      double t,
				      struct potentialArg * potentialArgs,
				      double vR,double vT,
				      double vz){
  double amp,t0,tf,d_ind, Rdist, Rorb; 
  double x,y,vx,vy, obj_x,obj_y,obj_z,obj_vx,obj_vy,obj_vz;
  double costheta, sintheta, phi2;
  double vxd, vyd, vzd, vR1, vT1;
  double * args= potentialArgs->args;
  //Get args
  amp= *args;
  t0= *(args+1);
  tf= *(args+2);
  d_ind= (t-t0)/(tf-t0);
  constrain_range3(&d_ind);

  // convert to cartesian coodinates
  x= R*cos(phi);
  y= R*sin(phi);
  vx=  vR*cos(phi) - vT*sin(phi);
  vy=  vR*sin(phi) + vT*cos(phi);

  // Interpolate x, y, z, vx, vy, vz
  obj_x= gsl_spline_eval(*potentialArgs->spline1d,d_ind,*potentialArgs->acc1d);
  obj_y= gsl_spline_eval(*(potentialArgs->spline1d+1),d_ind,
			 *(potentialArgs->acc1d+1));
  obj_z= gsl_spline_eval(*(potentialArgs->spline1d+2),d_ind,
			 *(potentialArgs->acc1d+2));
  obj_vx= gsl_spline_eval(*(potentialArgs->spline1d+3),d_ind,
			 *(potentialArgs->acc1d+3));
  obj_vy= gsl_spline_eval(*(potentialArgs->spline1d+4),d_ind,
			 *(potentialArgs->acc1d+4));
  obj_vz= gsl_spline_eval(*(potentialArgs->spline1d+5),d_ind,
			 *(potentialArgs->acc1d+5));

  // convert velocity to reference frame of LMC
  vxd = vx - obj_vx;
  vyd = vy - obj_vy;
  vzd = vz - obj_vz;

  // determine angles for transformation
  Rdist= pow(pow(x-obj_x, 2)+pow(y-obj_y, 2), 0.5);
  Rorb = pow(obj_x*obj_x + obj_y*obj_y, 0.5);
  phi2 = atan2(y-obj_y, x-obj_x);           
  cosinerule(Rorb,R,Rdist,&costheta,&sintheta);

  // convert velocity to cyclindrical coordinates in reference frame of LMC.
  vR1=    cos(phi2)*vxd    + sin(phi2)*vyd;        
  vT1=   -sin(phi2)*vxd    + cos(phi2)*vyd;

  // Calculate z force
  return amp * calczforce(Rdist,(z-obj_z),phi,t,potentialArgs->nwrapped,
			   potentialArgs->wrappedPotentialArg,vR1,vT1,vzd);
}

double MovingObjectDissipativephiforce(double R,double z,double phi,
					double t,
					struct potentialArg * potentialArgs,
					double vR,double vT,
					double vz){
  double amp,t0,tf,d_ind, Rdist, Rorb, RF, TF;
  double x,y,vx,vy, obj_x,obj_y,obj_z,obj_vx,obj_vy,obj_vz;
  double costheta, sintheta, phi2;
  double vxd, vyd, vzd, vR1, vT1;
  double * args= potentialArgs->args;
  //Get args
  amp= *args;
  t0= *(args+1);
  tf= *(args+2);
  d_ind= (t-t0)/(tf-t0);
  constrain_range3(&d_ind);

  // convert to cartesian coodinates
  x= R*cos(phi);
  y= R*sin(phi);
  vx=  vR*cos(phi) - vT*sin(phi);
  vy=  vR*sin(phi) + vT*cos(phi);

  // Interpolate x, y, z, vx, vy, vz
  obj_x= gsl_spline_eval(*potentialArgs->spline1d,d_ind,*potentialArgs->acc1d);
  obj_y= gsl_spline_eval(*(potentialArgs->spline1d+1),d_ind,
			 *(potentialArgs->acc1d+1));
  obj_z= gsl_spline_eval(*(potentialArgs->spline1d+2),d_ind,
			 *(potentialArgs->acc1d+2));
  obj_vx= gsl_spline_eval(*(potentialArgs->spline1d+3),d_ind,
			 *(potentialArgs->acc1d+3));
  obj_vy= gsl_spline_eval(*(potentialArgs->spline1d+4),d_ind,
			 *(potentialArgs->acc1d+4));
  obj_vz= gsl_spline_eval(*(potentialArgs->spline1d+5),d_ind,
			 *(potentialArgs->acc1d+5));

  // convert velocity to reference frame of LMC
  vxd = vx - obj_vx;
  vyd = vy - obj_vy;
  vzd = vz - obj_vz;

  // determine angles for transformation
  Rdist= pow(pow(x-obj_x, 2)+pow(y-obj_y, 2), 0.5);
  Rorb = pow(obj_x*obj_x + obj_y*obj_y, 0.5);
  phi2 = atan2(y-obj_y, x-obj_x);           
  cosinerule(Rorb,R,Rdist,&costheta,&sintheta);

  // convert velocity to cyclindrical coordinates in reference frame of LMC.
  vR1=    cos(phi2)*vxd    + sin(phi2)*vyd;        
  vT1=   -sin(phi2)*vxd    + cos(phi2)*vyd;

  // Calculate Radial and transverse forces
  RF= calcRforce(Rdist,(z-obj_z),phi,t,potentialArgs->nwrapped,
		 potentialArgs->wrappedPotentialArg,vR1,vT1,vzd);
  TF= calcPhiforce(Rdist,(z-obj_z),phi,t,potentialArgs->nwrapped,
		 potentialArgs->wrappedPotentialArg,vR1,vT1,vzd);

  // Transform force back into the refernce frame of the MW
  return  -amp*(-RF*sintheta  +TF*costheta);

}


