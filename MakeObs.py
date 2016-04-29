#!/usr/bin/env python

import esutil
import fitsio
import numpy as np
import os
import numpy.lib.recfunctions as rec
import suchyta_utils.balrog


def Get(dir, name, file, cols, pdir='/gpfs01/astro/workarea/esuchyta/Y1-files/', pname='Y1A1NEW_COADD_SPT'):

    sim = fitsio.read(file, ext=-1, columns=cols)
    p = suchyta_utils.balrog.Y1Processing(dir=pdir, systematics=pname)
    sim = suchyta_utils.balrog.Y1Dataset(sim, processing=p)

    outfile = os.path.join(dir, '%s-pseudo.fits'%(name))
    if os.path.exists(outfile):
        os.remove(outfile)
    sim.PseudoBenchmarkGalaxies()
    if not os.path.exists(os.path.dirname(outfile)):
        os.makedirs(os.path.dirname(outfile))
    f = fitsio.FITS(outfile, 'rw')
    f.write(sim.data)

    outfile = os.path.join(dir, '%s-benchmark.fits'%(name))
    if os.path.exists(outfile):
        os.remove(outfile)
    sim.ApplyDepth(cond='>', val=22)
    sim.BenchmarkColorCuts()
    sim.ApplyCut(key='mag_auto_i', cond='>', val=17.5) 
    sim.ApplyCut(key='mag_auto_i', cond='<', val=22.0) 
    if not os.path.exists(os.path.dirname(outfile)):
        os.makedirs(os.path.dirname(outfile))
    f = fitsio.FITS(outfile, 'rw')
    f.write(sim.data)

    return sim, p


if __name__=='__main__':
    name = 'y1a1_spto_s82'
   
    bands = ['g','r','i','z','y']
    cols = ['alphawin_j2000_i', 'deltawin_j2000_i', 'alphawin_j2000_g', 'deltawin_j2000_g', 'flux_auto_g', 'fluxerr_auto_g','field']
    for band in bands:
        cols.append('mag_auto_%s'%(band))
        cols.append('flux_radius_%s'%(band))
        cols.append('flags_%s'%(band))
    cols = suchyta_utils.balrog.AddModestNeed(cols, release='y1a1')
    scols = np.append(cols, ['annz_best', 'annz_best_err', 'intrainflag'])
    dcols = np.append(cols, ['annz2_z', 'annz2_z_err', 'intrain_flag'])

    simfile = os.path.join('FITS', 'balrog', '%s.fits'%(name))
    desfile = os.path.join('FITS', 'des', '%s.fits'%(name))

    outd = os.path.join('FITS','samples', name)
    sim, sp = Get(outd, 'sim', simfile, scols, pdir='/gpfs01/astro/workarea/esuchyta/Y1-files/', pname='Y1A1NEW_COADD_SPT')
    des, dp = Get(outd, 'des', desfile, dcols, pdir='/gpfs01/astro/workarea/esuchyta/Y1-files/', pname='Y1A1NEW_COADD_SPT')
