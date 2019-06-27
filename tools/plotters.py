"""
This module contains plotting routines
"""
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np





def draw_synoptic(ax, twiss):
    """
    Simple function to draw a synoptic as boxes
    uses the keyword column of a dataframe containing twiss

    Parameters
    ----------
    ax : matplotlib.axes
        Axes where the synoptic is to be drawn

    twiss : DataFrame
        pandas dataframe containing the lattice and in
        in particular the S, L, K1L and KEYWORD columns

    """

    for _, row in twiss.iterrows():

        if row['KEYWORD'] == 'QUADRUPOLE':
            _ = ax.add_patch(
                mpl.patches.Rectangle(
                    (row['S']-row['L'], 0), row['L'], np.sign(row['K1L']),
                    facecolor='k', edgecolor='k'))
        elif row['KEYWORD'] == 'SBEND':
            _ = ax.add_patch(
                mpl.patches.Rectangle(
                    (row['S']-row['L'], -1), row['L'], 2,
                    facecolor='None', edgecolor='k'))
        elif row['KEYWORD'] == 'SEXTUPOLE':
            if row['K2L'] == 0:
                _ = ax.add_patch(
                   mpl.patches.Rectangle(
                       (row['S']-row['L'], 0), row['L'], 0.6,
                       facecolor='grey', edgecolor='grey'))
            else:
                _ = ax.add_patch(
                   mpl.patches.Rectangle(
                       (row['S']-row['L'], 0), row['L'], 0.6*np.sign(row['K2L']),
                       facecolor='grey', edgecolor='grey'))
        elif row['KEYWORD'] == 'HKICKER':
            _ = ax.add_patch(
                mpl.patches.Rectangle(
                    (row['S']-row['L'], -1), row['L'], 2,
                    facecolor='lightgrey', edgecolor='lightgrey'))


def draw_aperture(ax, twiss, aper):
    """
    Simple function to draw apertures rectangles onto a plot

    Parameters
    ----------
    ax : matplotlib.axes
        Axes where the synoptic is to be drawn

    twiss : DataFrame
        pandas dataframe containing the lattice and in
        in particular the S, L, K1L and KEYWORD columns

    aper : string
        Either 'APER_1' for horizontal or 'APER_2' for vertical
        possibly another column
    """
    for _, row in twiss[twiss[aper] > 0].iterrows():
        _ = ax.add_patch(
            mpl.patches.Rectangle(
                (row['S']-row['L'], 1), row['L'], -1+row[aper],
                facecolor='k'))
        _ = ax.add_patch(
            mpl.patches.Rectangle(
                (row['S']-row['L'], -1), row['L'], 1-row[aper],
                facecolor='k'))



def beam_size(beta, dispersion, eps, dpp, n):
    """
    Simple calculation of beam size
    """
    beam = np.sqrt(eps*beta + (dpp*dispersion)**2)
    return beam*n


def draw_optics(header, twiss, fig,nsig=2,     
                keywords_to_label=['QUADRUPOLE', 'RBEND'],
                betlim = None, Dlim = None, Hsizelim = None, mulim = None,
                synoptic_drawer=draw_synoptic,
                aperture_drawer=draw_aperture,
                beam_sizer=beam_size):
    """
    plots some kind of optics summary using twiss file

    Parameters
    ----------

    header : dic
        Dictionary of the twiss file header, returned by pybt.tools.parsers.read_twiss_file
        inparticular uses EX, EY, SIGE

    twiss : DataFrame
        pandas dataframe containing the lattice and in particular
        BETX, BETY, DX, DY

    keywords_to_label : optional, list of strings or None
        list of KEYWORDS to be put labelled on top of the plot
        if None then nothing is written

    synoptic_drawer : optional, callable
        like pybt.tools.plotters.draw_synoptic

    aperture_drawer : optional, callable
        like pybt.tools.plotters.draw_aperture

    Examples
    --------

    Here we show the output of this function, in the case of the TT20 lattice

    >>> import matplotlib.pyplot as plt
    >>> header, twiss = pybt.read_twiss_file('../data/twiss.tfs')
    >>> fig = plt.figure(figsize=(15,7))
    >>> pybt.tools.plotters.draw_optics(header, twiss, fig, ['QUADRUPOLE'])
    >>> fig.axes[3].set_ylim(fig.axes[4].get_ylim())

    """
    gs = mpl.gridspec.GridSpec(5, 1, height_ratios=[1, 4, 4,4,4])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)
    ax3 = fig.add_subplot(gs[2], sharex=ax1)
    ax4 = fig.add_subplot(gs[3], sharex=ax1)
    ax5 = fig.add_subplot(gs[4], sharex=ax1)

    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax3.get_xticklabels(), visible=False)
    plt.setp(ax4.get_xticklabels(), visible=False)
    plt.setp(ax4.get_xticklabels(), visible=True)

    # top plot is synoptic
    ax1.axis('off')
    ax1.set_ylim(-1.2, 1)
    ax1.plot([0, twiss['S'].max()], [0, 0], 'k-')
    ax2.grid()
    ax3.grid()
    ax4.grid()
    ax5.grid()

    synoptic_drawer(ax1, twiss)

    #2nd plot is beta functions
    ax2.set_ylabel(r'$\beta$ (m)')
    ax2.plot(twiss['S'], twiss['BETX'], 'r-', linewidth = 3, label = r'$\beta_x$')
    ax2.plot(twiss['S'], twiss['BETY'], 'b-', linewidth = 3, label = r'$\beta_y$')
    ax2.legend(loc ='lower right')
    if betlim:
        ax2.set_ylim(0,betlim)

    #3rd plot is dispersion functions
    ax3.set_ylabel('D (m)')
    ax3.plot(twiss['S'], twiss['DX'], 'r-', linewidth = 3,label = r'$D_x$')
    ax3.plot(twiss['S'], twiss['DY'], 'b-', linewidth = 3, label = r'$D_y$')
    ax3.legend(loc ='lower right')
    if Dlim:
        ax3.set_ylim(0,Dlim)

    #4th plot is the horizontal beam size, with APER_1
    ax4.set_ylabel('Beam size (m)')
    ax4.plot(twiss['S'],
             beam_sizer(twiss['BETX'], twiss['DX'], header['EX'], header['SIGE'], nsig), 'r-')
    ax4.plot(twiss['S'],
             -beam_sizer(twiss['BETX'], twiss['DX'], header['EX'], header['SIGE'], nsig), 'r-', label = r'$\pm$' + str(nsig) + r'$*\sigma_{x,rms}$')
    ax4.fill_between(twiss['S'], beam_sizer(twiss['BETX'], twiss['DX'], header['EX'], header['SIGE'], nsig),-beam_sizer(twiss['BETX'], twiss['DX'], header['EX'], header['SIGE'], nsig),facecolor='red', alpha = 0.2,  interpolate=True)
    ax4.plot(twiss['S'],
             beam_sizer(twiss['BETY'], twiss['DY'], header['EY'], header['SIGE'], nsig), 'b-', label = r'$\pm$' +str(nsig) +  r'$*\sigma_{y,rms}$')
    ax4.plot(twiss['S'],
             -beam_sizer(twiss['BETY'], twiss['DY'], header['EY'], header['SIGE'], nsig), 'b-')
    ax4.fill_between(twiss['S'], beam_sizer(twiss['BETY'], twiss['DY'], header['EY'], header['SIGE'], nsig),-beam_sizer(twiss['BETY'], twiss['DY'], header['EY'], header['SIGE'], nsig),facecolor='blue', alpha = 0.2,   interpolate=True)
    ax4.legend(loc ='lower right')
    if Hsizelim:
        ax4.set_ylim(-Hsizelim,Hsizelim)

#    aperture_drawer(ax4, twiss, 'APER_1')

    #4th plot is the horizontal beam size, with APER_1
#    ax5.set_ylabel('V size (m)')#
#    ax5.plot(twiss['S'],
#             beam_sizer(twiss['BETY'], twiss['DY'], header['EY'], header['SIGE'], 4), 'b-')
#    ax5.plot(twiss['S'],
#             -beam_sizer(twiss['BETY'], twiss['DY'], header['EY'], header['SIGE'], 4), 'b-')
#    ax5.fill_between(twiss['S'], beam_sizer(twiss['BETY'], twiss['DY'], header['EY'], header['SIGE'], 4),-beam_sizer(twiss['BETY'], twiss['DY'], header['EY'], header['SIGE'], 4),facecolor='blue', alpha = 0.2,   interpolate=True)

#    aperture_drawer(ax5, twiss, 'APER_2')

    #5nd plot is phase advance
    ax5.set_ylabel(r'$\mu$ ($2*\pi$)')
    ax5.plot(twiss['S'], twiss['MUX'], 'r-', linewidth = 3, label = r'$\mu_x$')
    ax5.plot(twiss['S'], twiss['MUY'], 'b-', linewidth = 3, label = r'$\mu_y$')
    ax5.legend(loc ='lower right')
    if mulim:
        ax5.set_ylim(-mulim,mulim)

    if keywords_to_label:
        if not isinstance(keywords_to_label, list):
            raise TypeError('keywords_to_label has to be a list, if provided')

        axnames = ax1.twiny()
        axnames.spines['top'].set_visible(False)
        axnames.spines['left'].set_visible(False)
        axnames.spines['right'].set_visible(False)
        ax1._shared_x_axes.join(ax1, axnames)

        ticks, ticks_labels = list(), list()
        for keyword in keywords_to_label:
            sub_twiss = twiss[twiss['KEYWORD'] == keyword]
            ticks += list(sub_twiss['S'])
            ticks_labels += list(sub_twiss['NAME'])

        axnames.set_xticks(ticks)
        axnames.set_xticklabels(ticks_labels, rotation=90)



    ax5.set_xlabel('s (m)')
    #ax4.set_xticks(np.arange(0,np.round(np.max(twiss['S'])),2))
    #print(str(np.arange(0,np.round(np.max(twiss['S'])),2)))
#    ax4.set_xticklabels(str(np.arange(0,np.round(np.max(twiss['S'])),2)))

#    fig.tight_layout()
    fig.subplots_adjust(hspace=0)
    ax1.set_xlim(0, twiss['S'].max())
