
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tools.parsers as pa
import tools.plotters as pl

filename = 'twiss_40'
mytwiss = 'new/' +  filename + '.out'

header, tfs = pa.read_twiss_file(mytwiss)
gamma_top = header['GAMMA']
gamma_inj = 1.0074626866
header['SIGE'] = 0.0001
header['EX'] = 2E-6/np.sqrt(1-1/gamma_inj**2)/gamma_inj
header['EY'] = 2E-6/np.sqrt(1-1/gamma_inj**2)/gamma_inj
tfs['DX'] = tfs['DX']*np.sqrt(1-1/gamma_top**2)
fig = plt.figure(figsize=(12,9))
#fig.axes[3].set_ylim(fig.axes[4].get_ylim())

nsig = 2
pl.draw_optics(header, tfs, fig, nsig, ['QUADRUPOLE', 'SBEND','SEXTUPOLE','HKICKER'],
                betlim = 7., Dlim = 4.5, Hsizelim = 0.025, mulim = 2,
                )
plt.show()
fig.savefig(filename+'.png')
