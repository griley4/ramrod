"""
This module contains simple tracking routines that use the RE first order
matrix in the twiss output of MADX
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


MAD_SPACE = {'X':1, 'PX':2, 'Y':3, 'PY':4, 'T':5, 'PT':6}


def transport_part(twiss, element, dim):
    """
    returns a function to transport a particle
    from the start of the lattice to element
    using the 6x6 matrix in twiss

    Parameters
    ----------
    twiss : DataFrame
        twiss DataFrame of the lattice
        Must contain the transfer matrix element RE

    element : string
        name of the element (or its firs occurence)

    dim : string
        dimension to be returned by the transport
        X PX Y PY T or PT

    Returns
    -------
    function
        Returns a function that transports particle coordinates


    Examples
    --------

    This example illustrates the transport of particles

    >>> import pandas as pd
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> header, twiss = pybt.read_twiss_file('../data/twiss.tfs')
    >>> distribution = pd.DataFrame(
    >>>     np.hstack((
    >>>         pybt.tools.particles.beam_distrib_norm(*twiss.iloc[0][['ALFX', 'BETX']], header['EX'], 10000, 3),
    >>>         pybt.tools.particles.beam_distrib_norm(*twiss.iloc[0][['ALFY', 'BETY']], header['EY'], 10000, 3),
    >>>         np.zeros((10000, 2))
    >>>     )),
    >>>     columns=['X', 'PX', 'Y', 'PY', 'T', 'PT'])
    >>> 
    >>> distribution['X_tcsc'] =  distribution.apply(
    >>>     pybt.tools.twiss_track.transport_part(twiss, 'TCSC.UP', 'X'),
    >>>     axis=1)
    >>> distribution['Y_tcsc'] =  distribution.apply(
    >>>     pybt.tools.twiss_track.transport_part(twiss, 'TCSC.UP', 'Y'),
    >>>     axis=1)
    >>> plt.plot(distribution['X_tcsc']*1e3, distribution['Y_tcsc']*1e3, 'ro', markersize=0.2, alpha=0.3)
    >>> plt.plot(distribution['X']*1e3, distribution['Y']*1e3, 'bo', markersize=0.2, alpha=0.3)
    >>> plt.xlabel('x (mm)')
    >>> plt.ylabel('y (mm)')
    >>> plt.axis('equal')
    >>> plt.show()


    """

    required_columns = ['RE{:d}{:d}'.format(dim, col)
                        for col in range(1, 7)
                        for dim in range(1, 7)]

    if not all([re in twiss.columns for re in required_columns]):
        raise ValueError('Some column RE are missing in the twiss file')
    try:
        twiss_elem = twiss[twiss['NAME'] == element].iloc[0]
    except IndexError:
        raise ValueError('provided element looks missing from twiss')


    def transport_function(row):
        val = (
            row['X'] *twiss_elem['RE{:d}1'.format(MAD_SPACE[dim])] +
            row['PX']*twiss_elem['RE{:d}2'.format(MAD_SPACE[dim])] +
            row['Y'] *twiss_elem['RE{:d}3'.format(MAD_SPACE[dim])] +
            row['PY']*twiss_elem['RE{:d}4'.format(MAD_SPACE[dim])] +
            row['T'] *twiss_elem['RE{:d}5'.format(MAD_SPACE[dim])] +
            row['PT']*twiss_elem['RE{:d}6'.format(MAD_SPACE[dim])] +
            twiss_elem[dim])
        return val

    return transport_function




def track_with_re(twiss, dimension, vector):
    """
    Tracks a given particle vector along the twiss lattice
    using the transfer matrix from the twiss file
    Note that it also add the central trajectory to the transported values

    Parameters
    ----------
    twiss : DataFrame
        twiss DataFrame of the lattice
        Must contain the transfer matrix element RE

    dimension : string
        dimension to be returned
        X, PX, Y, PY, T or PT

    vector : 6-long 1D-array
        input particle vector at the start of lattice

    Returns
    -------
    ndarray
        array the length of the twiss file with the trajectory asked

    Examples
    --------

    Illustrating the use of this tracker together with other functions of the package

    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> header, twiss = pybt.read_twiss_file('../data/twiss.tfs')
    >>> fig = plt.figure(figsize=(15,7))
    >>> pybt.tools.plotters.draw_optics(header, twiss, fig, ['QUADRUPOLE'])
    >>> 
    >>> particlesx = pybt.tools.particles.beam_distrib_norm(
    >>>     *twiss.iloc[0][['ALFX', 'BETX']], header['EX'], 10, [4-1e-6,4+1e-6])
    >>> 
    >>> particlesy = pybt.tools.particles.beam_distrib_norm(
    >>>     *twiss.iloc[0][['ALFY', 'BETY']], header['EY'], 10, [4-1e-6,4+1e-6])
    >>> 
    >>> _ = fig.axes[3].autoscale(False), fig.axes[4].autoscale(False)
    >>> 
    >>> for particle in np.hstack((particlesx, particlesy, np.zeros((10,2)))):
    >>>     fig.axes[3].plot(twiss['S'],
    >>>                     pybt.tools.twiss_track.track_with_re(twiss, 'X', particle),
    >>>                     'r-', alpha=0.2)
    >>>     fig.axes[4].plot(twiss['S'],
    >>>                     pybt.tools.twiss_track.track_with_re(twiss, 'Y', particle),
    >>>                     'b-', alpha=0.2)



    """


    required_columns = ['RE{:d}{:d}'.format(MAD_SPACE[dimension], col)
                        for col in range(1, 7)]

    if not all([re in twiss.columns for re in required_columns]):
        raise ValueError('Some column RE are missing in the twiss file')

    if not dimension in twiss.columns:
        raise ValueError('Some column are missing for the central trajectory')

    if not dimension in MAD_SPACE.keys():
        raise ValueError('Wrong dimension input, does not exist :', dimension)


    array = [vector[col]*twiss[required_columns[col]] for col in range(6)]
    array = np.sum(np.array(array).transpose(), axis=1)

    array += twiss[dimension]

    return array
