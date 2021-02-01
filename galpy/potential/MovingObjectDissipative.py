###############################################################################
#   MovingObjectDissipative.py: class that implements the disspative forces coming from
#                             a moving object
#
#   To implement a moving dissipative framework, we are in need of two coordinate transformations.
#   First - transform the velocities of the satellite to the coordinates system of the LMC-like object.
#   Second - transform the dissipative forces back to the coordinates system of the Galaxy.
#   The easiest way to do these transformations is as follows: 
#   The z coordinate is simple
#   The radial and the transverse directions can be interchanged between the two frames using a rotation matrix.
#   The angle of the rotation matrix can be estimated using the cosine rule.
# 
###############################################################################
import copy
import numpy
from .Potential import Potential, _isNonAxi, flatten, \
    evaluatePotentials, evaluateRforces, evaluatephiforces, evaluatezforces, evaluateDensities, _check_c
from .DissipativeForce import DissipativeForce
from .PlummerPotential import PlummerPotential
class MovingObjectDissipative(DissipativeForce):
    """
    Class that implements the dissipative force coming from a moving object by combining
    any galpy dissipative force with an integrated galpy orbit.
    """
    def __init__(self,orbit,pot=None,amp=1.0,
                 ro=None,vo=None):
        """
        NAME:

           __init__

        PURPOSE:

           initialize a MovingObjectPotential

        INPUT:

           orbit - the Orbit of the object (Orbit object)

           pot - A potential object or list of potential objects representing the potential of the moving object; should be spherical, but this is not checked 
           
           amp (=1.) another amplitude to apply to the potential

           ro=, vo= distance and velocity scales for translation into internal units (default from configuration file)

        OUTPUT:

           (none)

        HISTORY:

           2011-04-10 - Started - Bovy (NYU)

           2018-10-18 - Re-implemented to represent general object potentials using galpy potential models - James Lane (UofT)

           2021-1-13 - Changed to dissipative force

        """

        DissipativeForce.__init__(self,amp=amp,ro=ro,vo=vo,amp_units='mass')       
        # If no potential supplied use a default Plummer sphere
        if pot is None:
            raise NotImplementedError('MovingObjectDissipative should be provided a dissipative force')
        self._pot=pot
        self._orb= copy.deepcopy(orbit)
        self._orb.turn_physical_off()
        self.isNonAxi= True
        self.hasC= _check_c(self._pot)
        return None


    def _Rforce(self,R,z,phi=0.,t=0.,v=None):
        """
        NAME:
           _Rforce
        PURPOSE:
           evaluate the radial force for this potential
        INPUT:
           R - Galactocentric cylindrical radius
           z - vertical height
           phi - azimuth
           t - time
           v = current velocity in cylindrical coordinates
        OUTPUT:
           the radial force
        HISTORY:
           2011-04-10 - Written - Bovy (NYU)
           2018-10-18 - Updated for general object potential - James Lane (UofT)
           2021-1-13 - Changed to dissipative force  
        """
        # Cylindrical distance
        Rdist = _cylR(R,phi,self._orb.R(t),self._orb.phi(t))
        cosphi,sinphi=_cosinerule(self._orb.R(t),R,Rdist) 

        # Difference vector
        orbz = self._orb.z(t) if self._orb.dim() == 3 else 0
        (xd,yd,zd) = _cyldiff(R,phi,z,self._orb.R(t), self._orb.phi(t), orbz)

        # convert sat into cartesian coordinates 
        vx=v[0]*numpy.cos(phi)  - v[1]*numpy.sin(phi)
        vy=v[0]*numpy.sin(phi)  + v[1]*numpy.cos(phi)
        # convert LMC into cartesian coordinates
        vx1=   numpy.cos(self._orb.phi(t))*self._orb.vR(t)   - numpy.sin(self._orb.phi(t))*self._orb.vT(t)
        vy1=   numpy.sin(self._orb.phi(t))*self._orb.vR(t)   + numpy.cos(self._orb.phi(t))*self._orb.vT(t) 
        # In cartesian frame of the LMC
        vxd=vx-vx1
        vyd=vy-vy1
        # Convert into cylindrical coordinates of the LMC
        phi2=numpy.arctan2(yd,xd)
        vn=numpy.zeros_like(v)
        vn[0]=numpy.cos(phi2)*vxd + numpy.sin(phi2)*vyd
        vn[1]=-numpy.sin(phi2)*vxd + numpy.cos(phi2)*vyd
        vn[2]=v[2]-self._orb.vz(t)


        #Evaluate cylindrical radial and tangential force
        RF = evaluateRforces(self._pot,Rdist,zd,t=t,v=vn,use_physical=False)
        TF = evaluatephiforces(self._pot,Rdist,zd,t=t,v=vn,use_physical=False)
 
        # Return R force using appropriate rotation
        return -1.0*(RF*cosphi + TF*sinphi)

    def _zforce(self,R,z,phi=0.,t=0.,v=None):
        """
        NAME:
           _zforce
        PURPOSE:
           evaluate the vertical force for this potential
        INPUT:
           R - Galactocentric cylindrical radius
           z - vertical height
           phi - azimuth
           t - time
           v = current velocity in cylindrical coordinates
        OUTPUT:
           the vertical force
        HISTORY:
           2011-04-10 - Written - Bovy (NYU)
           2018-10-18 - Updated for general object potential - James Lane (UofT)
           2021-1-13 - Changed to dissipative force  
        """
        # Cylindrical distance
        Rdist = _cylR(R,phi,self._orb.R(t),self._orb.phi(t))
        cosphi,sinphi=_cosinerule(self._orb.R(t),R,Rdist) 

        # Difference vector
        orbz = self._orb.z(t) if self._orb.dim() == 3 else 0
        (xd,yd,zd) = _cyldiff(R,phi,z,self._orb.R(t), self._orb.phi(t), orbz)
            

        # convert sat into cartesian coordinates 
        vx=v[0]*numpy.cos(phi)  - v[1]*numpy.sin(phi)
        vy=v[0]*numpy.sin(phi)  + v[1]*numpy.cos(phi)
        # convert LMC into cartesian coordinates
        vx1=   numpy.cos(self._orb.phi(t))*self._orb.vR(t)   - numpy.sin(self._orb.phi(t))*self._orb.vT(t)
        vy1=   numpy.sin(self._orb.phi(t))*self._orb.vR(t)   + numpy.cos(self._orb.phi(t))*self._orb.vT(t) 
        # In cartesian frame of the LMC
        vxd=vx-vx1
        vyd=vy-vy1
        # Convert into cylindrical coordinates of the LMC
        phi2=numpy.arctan2(yd,xd)
        vn=numpy.zeros_like(v)
        vn[0]=numpy.cos(phi2)*vxd + numpy.sin(phi2)*vyd
        vn[1]=-numpy.sin(phi2)*vxd + numpy.cos(phi2)*vyd
        vn[2]=v[2]-self._orb.vz(t)



        #Evaluate and return z force
        return evaluatezforces(self._pot,Rdist,zd,t=t, v=vn, use_physical=False)

    def _phiforce(self,R,z,phi=0.,t=0.,v=None):
        """
        NAME:
           _phiforce
        PURPOSE:
           evaluate the azimuthal force for this potential
        INPUT:
           R - Galactocentric cylindrical radius
           z - vertical height
           phi - azimuth
           t - time
           v = current velocity in cylindrical coordinates
        OUTPUT:
           the azimuthal force
        HISTORY:
           2011-04-10 - Written - Bovy (NYU)
           2018-10-18 - Updated for general object potential - James Lane (UofT)
           2021-1-13 - Changed to dissipative force
        """
        # Cylindrical distance
        Rdist = _cylR(R,phi,self._orb.R(t),self._orb.phi(t))
        cosphi,sinphi=_cosinerule(self._orb.R(t),R,Rdist) 

        # Difference vector
        orbz = self._orb.z(t) if self._orb.dim() == 3 else 0
        (xd,yd,zd) = _cyldiff(R,phi,z,self._orb.R(t), self._orb.phi(t), orbz)

        # convert sat into cartesian coordinates 
        vx=v[0]*numpy.cos(phi)  - v[1]*numpy.sin(phi)
        vy=v[0]*numpy.sin(phi)  + v[1]*numpy.cos(phi)
        # convert LMC into cartesian coordinates
        vx1=   numpy.cos(self._orb.phi(t))*self._orb.vR(t)   - numpy.sin(self._orb.phi(t))*self._orb.vT(t)
        vy1=   numpy.sin(self._orb.phi(t))*self._orb.vR(t)   + numpy.cos(self._orb.phi(t))*self._orb.vT(t) 
        # In cartesian frame of the LMC
        vxd=vx-vx1
        vyd=vy-vy1
        # Convert into cylindrical coordinates of the LMC
        phi2=numpy.arctan2(yd,xd)
        vn=numpy.zeros_like(v)
        vn[0]=numpy.cos(phi2)*vxd + numpy.sin(phi2)*vyd
        vn[1]=-numpy.sin(phi2)*vxd + numpy.cos(phi2)*vyd
        vn[2]=v[2]-self._orb.vz(t)
        
        #Evaluate cylindrical radial and tangential force.
        RF = evaluateRforces(self._pot,Rdist,zd,t=t, v=vn, use_physical=False)
        TF = evaluatephiforces(self._pot,Rdist,zd,t=t, v=vn, use_physical=False)
        # Return phi force using appropriate rotation.
        return -1.0*(-RF*sinphi + TF*cosphi)


def _cylR(R1,phi1,R2,phi2):
    return numpy.sqrt(R1**2.+R2**2.-2.*R1*R2*numpy.cos(phi1-phi2)) # Cosine law

def _cyldiff(R1,phi1,z1,R2,phi2,z2):
    dx = R1*numpy.cos(phi1)-R2*numpy.cos(phi2)
    dy = R1*numpy.sin(phi1)-R2*numpy.sin(phi2)
    dz = z1-z2
    return (dx,dy,dz)


def _cosinerule(R1,R2,Rdist):
    """
    R1 is the distance of the LMC from the MW
    R2 is the distance of the Sat from the MW
    Rdist is the estimated distance of the Sat from the LMC.
    """

    cosphi=(R1**2. - R2**2. - Rdist**2)/(2*R2*Rdist)
    sinphi=numpy.sqrt(1-cosphi**2)

    return (cosphi, sinphi)

