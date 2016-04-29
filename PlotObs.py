#!/usr/bin/env python

import fitsio
import numpy as np
import os
import numpy.lib.recfunctions as rec
import suchyta_utils.balrog

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import suchyta_utils.plot


if __name__=='__main__':
    suchyta_utils.plot.Setup()
    simfile = os.path.join('FITS','samples','y1a1_spto_s82','sim-pseudo.fits')
    sim = fitsio.read(simfile, ext=1)

    map_fig, map_ax = plt.subplots(1,1, figsize=(13,9))
    dims = {'xoffset':0, 'yoffset':-0.5, 'cmap':cm.YlOrRd, 'dims':[150,80], 'center':[10,-35], 'ra':'alphawin_j2000_i', 'dec':'deltawin_j2000_i', 'meridians':np.arange(0,360,20), 'parallels':np.arange(-70, 0, 10)}
    suchyta_utils.plot.MapPlot(ax=map_ax, fig=map_fig, cat=sim, **dims)

    print np.sum(sim['field']==0)
    print np.sum(sim['field']==1)
    plt.tight_layout()
    plt.savefig('plots/photo-z_coverage_429.png')
    plt.show()
