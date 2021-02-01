###############################################################################
#   LMCDynamicalFrictionForce: Class that implements the 
#                                        dynamical friction of the LMC
###############################################################################
import copy
import hashlib
import numpy
from scipy import special, interpolate
from ..util import bovy_conversion
from .DissipativeForce import DissipativeForce
from .Potential import _APY_LOADED, evaluateDensities, _check_c
from .Potential import flatten as flatten_pot
if _APY_LOADED:
    from astropy import units
_INVSQRTTWO= 1./numpy.sqrt(2.)
_INVSQRTPI= 1./numpy.sqrt(numpy.pi)
class LMCDynamicalFrictionForce(DissipativeForce):
    """Class that implements the LMC dynamical friction force

    .. math::


       \\mathbf{F}(\\mathbf{x},\\mathbf{v}) = -2\\pi\\,[G\\,M]\\,[G\\,\\rho(\\mathbf{x})]\\,\\ln[1+\\Lambda^2] \\,\\left[\\mathrm{erf}(X)-\\frac{2X}{\\sqrt{\\pi}}\\exp\\left(-X^2\\right)\\right]\\,\\frac{\\mathbf{v}}{|\\mathbf{v}|^3}\\,

    on a mass (e.g., a satellite galaxy or a black hole) :math:`M` at position :math:`\\mathbf{x}` moving at velocity :math:`\\mathbf{v}` through a background density :math:`\\rho`. The quantity :math:`X` is the usual :math:`X=|\\mathbf{v}|/[\\sqrt{2}\\sigma_r(r)`. The factor :math:`\\Lambda` that goes into the Coulomb logarithm is taken to be

    .. math::

       \\Lambda = \\frac{r/\\gamma}{\\mathrm{max}\\left(r_{\\mathrm{hm}},GM/|\\mathbf{v}|^2\\right)}\\,,

    where :math:`\\gamma` is a constant. This :math:`\\gamma` should be the absolute value of the logarithmic slope of the density :math:`\\gamma = |\\mathrm{d} \\ln \\rho / \\mathrm{d} \\ln r|`, although for :math:`\\gamma<1` it is advisable to set :math:`\\gamma=1`. Implementation here roughly follows `2016MNRAS.463..858P <http://adsabs.harvard.edu/abs/2016MNRAS.463..858P>`__ and earlier work.

    """
    def __init__(self,amp=1.,GMs=.1,gamma=0.3,
                 ro=None,vo=None):
        """
        NAME:

           __init__

        PURPOSE:

           initialize a Chandrasekhar Dynamical Friction force

        INPUT:

           amp - amplitude to be applied to the potential (default: 1)

           GMs - satellite mass; can be a Quantity with units of mass or Gxmass; can be adjusted after initialization by setting obj.GMs= where obj is your ChandrasekharDynamicalFrictionForce instance (note that the mass of the satellite can *not* be changed simply by multiplying the instance by a number, because he mass is not only used as an amplitude)

           gamma - Free-parameter in :math:`\\Lambda`



        OUTPUT:

           (none)

        HISTORY:

           2021-01-13 - Started - Richard D'Souza

        """
        DissipativeForce.__init__(self,amp=amp,ro=ro,vo=vo,
                                  amp_units='mass')
        self._gamma= gamma
        self._ms= self._amp/amp  # from handling in __init__ above, should be ms in galpy units 
        self.hasC= True
        self.isNonAxi = False
        return None


    def GMs(self,gms):
        if _APY_LOADED and isinstance(gms,units.Quantity):
            try:
                gms= gms.to(units.Msun).value\
                    /bovy_conversion.mass_in_msol(self._vo,self._ro)
            except units.UnitConversionError:
                # Try G x mass
                try:
                    gms= gms.to(units.pc*units.km**2/units.s**2)\
                        .value\
                        /bovy_conversion.mass_in_msol(self._vo,self._ro)\
                        /bovy_conversion._G
                except units.UnitConversionError:
                    raise units.UnitConversionError('GMs for %s should have units of mass or G x mass' % (type(self).__name__))
        self._ms= gms
        return None
    GMs= property(None,GMs)


    def _calc_force(self,R,phi,z,v,t):
        r= numpy.sqrt(R**2.+z**2.)
        lnLambda= self._gamma
        return -0.428*self._ms*lnLambda/r**2

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
           v= current velocity in cylindrical coordinates
        OUTPUT:
           the radial force
        HISTORY:
           2021-01-13 - Started - Richard D'Souza
        """
        return self._calc_force(R,phi,z,v,t)*v[0]
        

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
           v= current velocity in cylindrical coordinates
        OUTPUT:
           the azimuthal force
        HISTORY:
           2021-01-13 - Started - Richard D'Souza
        """
        return self._calc_force(R,phi,z,v,t)*v[1]*R
    

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
           v= current velocity in cylindrical coordinates
        OUTPUT:
           the vertical force
        HISTORY:
           2021-01-13 - Started - Richard D'Souza
        """
        return self._calc_force(R,phi,z,v,t)*v[2]

