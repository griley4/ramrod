"""
This module contains functions helping to generate and characterize particle distributions

"""
import numpy as np
from scipy.stats import truncexpon, truncnorm, uniform

def beam_distrib_norm(alpha, beta, eps, nparts, nsigma):
    """
    Generate a normally distributed particle distribution
    using the exponentially distributed actions.
    Allows cut distributions, including hollow ones.

    Parameters
    ----------
    alpha : scalar
        Twiss parameter
    beta : scalar
        Twiss parameter
    eps : scalar
        RMS emittance
    nparts : scalar
        Number of particles to be produced
    nsigma : length 2 array-like or scalar
       If scalar then we have the distribution only contains
       particles with 2J > nsigma**2*eps
       If array [n1,n2] then the distribution only contains
       oarticles with n1**2*eps > 2J > n2**2*eps
       The square is done so that the projection cut is in sigma

    Returns
    -------
    ndarray
        nparts*2 size array containing position and angle of each particle

    Examples
    --------


    Here we illustrate the function by plotting 2 distribution,
    from zero to 2 sigma on the left
    from 1 to 3 sigma on the right


    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> dist1 = pybt.tools.particles.beam_distrib_norm(0.5, 3, 3, 10000, 2)
    >>> dist2 = pybt.tools.particles.beam_distrib_norm(0.5, 3, 3, 10000, [1,3])
    >>> dist = np.array([(x,y) for x, y in dist1 if x<0] + [(x,y) for x, y in dist2 if x>0])
    >>> plt.plot(dist[:,0], dist[:,1], 'bo', markersize=0.3)
    >>> plt.xticks(np.arange(-2, 3.1, 1)*3)
    >>> plt.yticks(np.arange(-3, 3.1, 1)*np.sqrt(1.25))
    >>> plt.grid(True)
    >>> plt.show()

    """
    if np.size(nsigma) == 2:
        nsigma = np.sort(np.array(nsigma))
        cs_inv = truncexpon(b=nsigma[1]**2/2-nsigma[0]**2/2,
                            scale=2*eps).rvs(nparts)+nsigma[0]**2*eps
    elif np.size(nsigma) == 1:
        cs_inv = truncexpon(b=nsigma**2/2, scale=2*eps).rvs(nparts)
    else:
        raise ValueError('sigma is either length 2 array or scalar', nsigma)

    angles = uniform(0, 2*np.pi).rvs(nparts)

    x = np.sqrt(cs_inv*beta)*np.cos(angles)
    xp = -np.sqrt(cs_inv/beta)*(np.sin(angles)+alpha*np.cos(angles))

    return np.column_stack([x, xp])

def dpp_distrib_norm(dpp, nparts, nsigma):
    """
    Generate a 1D truncated Gaussian from a pdf with given standard deviation

    Parameters
    ----------
    dpp : scalar
        Standard deviation of the distribution
    nparts : scalar
        Number of particles to be produced
    nsigma : scalar
        sigma cut of the distribution, i.e. +- nsigma*dpp

    Returns
    -------
    ndarray
        nparts 1-D array containing the particles

    """
    return truncnorm.rvs(-nsigma, nsigma, scale=dpp, size=nparts)
