################################################################
#  The reflex motion or the acceleration of the LMC on the MW 
#  can be calculated by considering the amplitude and the direction.
#
#  The internal forces balance out each other F12 = -F21 
#  This implies that M1* a12 = - M2* a21
#  where a12 is the accleration of object 2 on object 1.
#  The amplitude of the acceleration on the LMC by the MW is - G M_MW/ R\*\*2. 
#  Hence, the amplitude of the reflex accleration on the MW is   G M_LMC/R\*\*2
#
#  To calculate the direction in cylindrical coordinates:
#  First, the direction of z remains the same.
#  Second, one has to calculate the projection along the R and the phi directions. 
#  To do this, first consider the projections of the amplitude of the acceleration 
#  along the x and y directions, and then consider a rotation along the new angle (phi2) 
#  to be considered.
###############################################################################

import numpy
import copy
from .Potential import Potential, _isNonAxi, flatten, evaluatePotentials, evaluateRforces, evaluatezforces, evaluateDensities, _check_c

def _cylR(R1,phi1,R2,phi2):
    return numpy.sqrt(R1**2.+R2**2.-2.*R1*R2*numpy.cos(phi1-phi2)) # Cosine law

def _cyldiff(R1,phi1,z1,R2,phi2,z2):
    dx = R1*numpy.cos(phi1)-R2*numpy.cos(phi2)
    dy = R1*numpy.sin(phi1)-R2*numpy.sin(phi2)
    dz = z1-z2
    return (dx,dy,dz)

class ReflexMotion(Potential):
    """
    Class that implements the reflex motion from a moving object by combining
    any galpy potential with an integrated galpy orbit.
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
           
           amp (=1.) another amplitude to apply to the potential, in effect it should M_LMC/M_MW

           ro=, vo= distance and velocity scales for translation into internal units (default from configuration file)

        OUTPUT:

           (none)

        HISTORY:

           2011-04-10 - Started - Bovy (NYU)

           2018-10-18 - Re-implemented to represent general object potentials using galpy potential models - James Lane (UofT)

        """
        Potential.__init__(self,amp=amp,ro=ro,vo=vo)       
        # If no potential supplied use a default Plummer sphere
        if pot is None:
            raise NotImplementedError('Potential has to be provided')
        else:
            pot=flatten(pot)
            if _isNonAxi(pot):
                raise NotImplementedError('MovingObjectPotential for non-axisymmetric potentials is not currently supported')
            self._pot=copy.deepcopy(pot)
            self._pot._amp=self._pot._amp*amp
        self._orb= copy.deepcopy(orbit)
        self._orb.turn_physical_off()
        self.isNonAxi= True
        self.hasC= _check_c(self._pot)
        return None
    
    def _Rforce(self,R,z,phi=0.,t=0.):
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
        OUTPUT:
           the radial force
        HISTORY:
           2011-04-10 - Written - Bovy (NYU)
           2018-10-18 - Updated for general object potential - James Lane (UofT)
        """
        #Evaluate cylindrical radial force
        orbz = self._orb.z(t) if self._orb.dim() == 3 else 0
        orbphi= self._orb.phi(t)
        RF = -evaluateRforces(self._pot,self._orb.R(t),orbz,t=t,use_physical=False)
        Fx = RF*numpy.cos(orbphi)
        Fy = RF*numpy.sin(orbphi)
        # Return R force, negative of radial vector to evaluation location.
        return Fx*numpy.cos(phi)+ Fy*numpy.sin(phi)

    
    def _zforce(self,R,z,phi=0.,t=0.):
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
        OUTPUT:
           the vertical force
        HISTORY:
           2011-04-10 - Written - Bovy (NYU)
           2018-10-18 - Updated for general object potential - James Lane (UofT)
        """
        return -evaluatezforces(self._pot,self._orb.R(t),self._orb.z(t),t=t,use_physical=False)

    
    def _phiforce(self,R,z,phi=0.,t=0.):
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
        OUTPUT:
           the azimuthal force
        HISTORY:
           2011-04-10 - Written - Bovy (NYU)
           2018-10-18 - Updated for general object potential - James Lane (UofT)
        """
        #Evaluate cylindrical radial force
        orbz = self._orb.z(t) if self._orb.dim() == 3 else 0
        orbphi= self._orb.phi(t)
        RF = -evaluateRforces(self._pot,self._orb.R(t),orbz,t=t,use_physical=False)
        Fx = RF*numpy.cos(orbphi)
        Fy = RF*numpy.sin(orbphi)
        # Return R force, negative of radial vector to evaluation location.
        return -Fx*numpy.sin(phi)+ Fy*numpy.cos(phi)
